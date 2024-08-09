function slie_calcM4LU_fromlfiw(dataset, varargin)
% Function to calculate series of monthly SLIEs:
% - min, mean, max, median, lower 10% and upper 10% 
%
% ** New Version ***
%       This does not involve stacking SLIE images to derive
%       monthly min, median, and max. Everything is calculated
%       from LFIW data, which allows us to deal with missing
%       SLIE regions from ice charts and satelite gaps
%
%       This is relatively straightforward except where LFIW is measured
%       by overlapping coastvectors, such as in Kotebue Sound. This 
%       requires an attempt to merge the two resulting SLIEs toward each
%       other and selection of the SLIE that covers the longest distance
%       with the greatest number of unique vertices along the SLIE
%
% ** Also Note: we no longer calculate modal SLIEs because they were always
%               very jagged and unreaslistic
%
%
% USAGE
%     slie_calcM5_fromlfiw
%     slie_calcM5_fromlfiw(dataset)
%     slie_calcM5_fromlfiw(__, 'krad', krad)
%     slie_calcM5_fromlfiw(__, 'zcol', zcol)
%     slie_calcM5_fromlfiw(__, 'ncol', ncol)
%     slie_calcM5_fromlfiw(__, 'lcol', lcol)
%     slie_calcM5_fromlfiw(__, 'yr0', yr0)
%     slie_calcM5_fromlfiw(__, 'yr1', yr1)
%     slie_calcM5_fromlfiw(__, 'asadir', asadir)
%     slie_calcM5_fromlfiw(__, 'outdir', outdir)
%     slie_calcM5_fromlfiw(__, 'cfgfile', cfgfile)
%
% INPUTS
%   dataset : string specifying which dataset to use.
%             See slie_get_basefolder.m for accepted values.
%             default: 'beau'
%   cfgfile : string specifying the name of the configuraion file stored
%             in the region folder containing essential information for
%             performing all analysis of SLIE data. 
%             See slie_readconfigfile.m for more details
%             default: 'slieconfig.txt'
%    asadir : string specifying folder name in region folder where
%             the data files from the all-seasons analysis are held
%             default: 'AllSeasonsAnalysis'
%    outdir : string specifying folder name inside the asadir where
%             the output files will be saved
%             default: 'MonthlySLIEs'
%      krad : scalar integer specifying radius for fattening the
%             line of pixels when drawing the mean SLIEs
%             default: 5
%      zcol : color to assign to region outside of landfast ice.
%             (i.e., background color)
%             default: [1,1,1] (white)
%      lcol : color to assign to land pixels
%             default: [0.8,0.9,0.8] (greeny light gray)
%      ncol : color to assign to no-data regions
%             (i.e. outside study region)
%             default: [0.9,0.9,0.9] (light gray)
%     shcol : color to assign to shadow regions where 
%             coast vecs don't reach
%             default: [0.9,1.0,0.9] (greeny lighter gray)
%       yr0 : the first year to include in analyis
%             default: first year in record for specified region
%       yr1 : the last year to include in analyis
%             default: last year in record for specified region
%
%
% OUTPUTS:
%    monthly images: a grayscale GeoTIFF is generated for each M4LU product
%                    for each month with SLIE data.
%    monthly stacks: RGB GeoTIFFs are generated showing the stacked MMM and 
%                    LMU images for each month
%         csv files: comma-separated ASCII files containg the coordinates
%                    defining the vertices of each monthly SLIE project
%
%
% Andy Mahoney October 2023


%  -----------------------------------------------------------------------
% Check for variable arguments specified with function call
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'krad'
            krad = varargin{1+1};
            a = a + 2;
        case 'zcol'
            zcol = varargin{a+1};
            a = a + 2;
        case 'lcol'
            lcol = varargin{a+1};
            a = a + 2;
        case 'ncol'
            ncol = varargin{a+1};
            a = a + 2;
        case 'shcol'
            shcol = varargin{a+1};
            a = a + 2;
        case 'yr0'
            yr0 = varargin{a+1};
            a = a + 2;
        case 'yr1'
            yr1 = varargin{a+1};
            a = a + 2;
        case 'asadir'
            asadir = varargin{a+1};
            a = a + 2;
        case 'outdir'
            outdir = varargin{a+1};
            a = a + 2;
        case 'cfgfile'
            cfgfile = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end    

% -------------------------------------------------------------
% Assign default paths if not specified in function call
if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023'; 
end

if ~exist('outdir', 'var')
    outdir = 'MonthlySLIEs'; 
end


% (default values for plotting-related variables are set below)


% -------------------------------------------------------------
% Figure out which dataset we're dealing with
% set the appropriate base path for the data
if ~exist('dataset', 'var')
    dataset = 'beau';
end

% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

% Assign folder for All Seasons Analysis
asapath = [sfolder filesep asadir];
if ~exist(asapath, 'file')
    mkdir(asapath);
end

% Set output folder
outpath = [asapath filesep outdir];
if ~exist(outpath, 'file')
    mkdir(outpath);
end

% Set folders for output CSV files
mapdir = [outpath filesep 'CSVfiles_mapcoords'];
if ~exist(mapdir, 'file')
    mkdir(mapdir);
end
lldir = [outpath filesep 'CSVfiles_latlons'];
if ~exist(lldir, 'file')
    mkdir(lldir);
end



% ---------------------------------------------------------------------
% Read the SLIE configuration file to get various useful things
cfg = slie_readconfigfile(sfolder, cfgfile);

% Read coast vectors file
cvecs = slie_readcvecsfile(cfg.cvecfile);

% Read land mask and get GeoTIFF information
[lmask, R] = geotiffread(cfg.lmaskfile);
gtinfo = geotiffinfo(cfg.lmaskfile);
geokeyDT = gtinfo.GeoTIFFTags.GeoKeyDirectoryTag;
% Remove any "unknown" fields from geokey
if isfield(geokeyDT, 'Unknown')
    geokeyDT = rmfield(geokeyDT, 'Unknown');
end
sz = size(lmask);

% Identify shadow objects in landmask and pixels adjacent to them
% - this info is used when filling polygons bounded by the SLIE
shadows = lmask == 64;
shadregs = regionprops(shadows, 'boundingbox', 'pixelidxlist', 'image');
Nshad = numel(shadregs);
shadadjpix = cell(Nshad);
for sh = 1:Nshad
    
    % Place feature image into large image 
    % with extra pixel all the way around
    shadimgX = zeros(size(shadregs(sh).Image)+2);
    shadimgX(2:end-1,2:end-1) = shadregs(sh).Image;
    
    % Grow object image with a 3x3 convolution
    shadconv = conv2(shadimgX, ones(3), 'same');
    
    % Identify the outer perimeter of feature
    outerperim = (shadconv > 0) & (shadimgX == 0);
    opi = find(outerperim);
    [opy,opx] = ind2sub(size(shadconv), opi);
    
    % Translate sub-image pixel index values into whole-image values
    % checking for oub-of-bounds (oob) pixels
    x0 = floor(shadregs(sh).BoundingBox(1));
    y0 = floor(shadregs(sh).BoundingBox(2));
    opxx = opx+x0-1;
    opyy = opy+y0-1;
    oob = (opxx < 1) | (opxx > sz(2)) | (opyy < 1) | opyy > sz(1);
    opxx = opxx(~oob);
    opyy = opyy(~oob);
    shadregs(sh).adjpix = sub2ind(sz, opyy, opxx);
    notalot = 0;
end

% Read coastal segment data
% - this info is used when closing polygons defined by the SLIE
csegs = slie_read_coastalsegmentsfile(cfg.coastsegmentfile);


% -----------------------------------------------------------------------
% Assign various parameters for output images

% Set land and nodata values as coded in landmask
landval = 255;
oobval = 128;
shadowval = 64;

nodataval = 111;

% Set colors and colormap
if ~exist('zcol', 'var')
    zcol = [1.0,1.0,1.0]; % No-slie color
end

if ~exist('lcol', 'var')
    lcol = [0.8,0.9,0.8]; % Land color
end

if ~exist('ncol', 'var')
    ncol = [0.9,0.9,0.9]; % Nodata color
end

if ~exist('shcol', 'var')
    shcol = [0.75,0.83,0.83]; % Coast vector shadow color
end

if ~exist('mincol', 'var')
    mincol = [0.0,0.0,1.0]; % Min SLIE color
end

if ~exist('medcol', 'var')
    medcol = [0.4,0.4,1.0]; % Median SLIE color
end

if ~exist('maxcol', 'var')
    maxcol = [0.8,0.8,1.0]; % Max SLIE color
end

if ~exist('meancol', 'var')
    meancol = [0.0,0.0,0.0]; % Mean SLIE color
end

cmap = [zcol; ...
        maxcol; ... 
        medcol; ...
        mincol; ...
        meancol; ...
        lcol; ...
        ncol; ...
        shcol];
        
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Smoothing filter
if ~exist('krad', 'var')
    krad = 5;
end
ksz = 1+(krad*2);
kx = repmat(-krad:krad,ksz,1);
ky = repmat((-krad:krad)',1,ksz);
kr = sqrt(kx.^2 + ky.^2);
skern = kr <= krad;
rval = 2*krad*cfg.pixsz;


% ========================================================================
% NOW GET TO BUSINESS

% ------------------------------------------------------------
% Read all-seasons lfiw data and sliedates
lfiwfile = [asapath filesep 'lfiw_all.csv'];
lfiw_all = read_sliedatafile(lfiwfile, '%f');
Nv = size(lfiw_all, 1);

sliedatefile = [asapath filesep 'sliedatenum_all.csv'];
lfiwdates = read_sliedatafile(sliedatefile, '%f');

% Determine month for each data point
[yr,mo,da] = datevec(lfiwdates);

if ~exist('yr0', 'var')
    yr0 = min(yr(:));
end
if ~exist('yr1', 'var')
    yr1 = max(yr(:));
end

% Exlude data outside of start/end years, if specified
% by setting values to NaN
keepi = ((yr > yr0) | ((yr == yr0) & (mo >= 10))) & ...
        ((yr < yr1) | ((yr == yr1) & (mo <= 7)));
lfiw_all(~keepi) = NaN;
lfiwdates(~keepi) = NaN;


% Pre-assign array for holding monthly min, median, mean, and max lfiws
lfiw_m4.min = zeros(Nv,12);
lfiw_m4.median = zeros(Nv,12);
lfiw_m4.mean = zeros(Nv,12);
lfiw_m4.max = zeros(Nv,12);
lfiw_m4.low = zeros(Nv,12);
lfiw_m4.upp = zeros(Nv,12);

m4 = fieldnames(lfiw_m4);

despike = 10;

wbinsz = 2; %km
wbins = 0:wbinsz:150;
wmid = (wbins(1:end-1)+wbins(2:end))/2;
for m=[1:7,10:12]
    % Create a mask identifying all lfiw values for this month
    % values for other months are assigned NaNs
    momask = 1*(mo == m);
    momask(momask == 0) = NaN;
    lfiw_mo = lfiw_all.*momask;

    if ~any(lfiw_mo(:))
        continue;
    end
    
    % Calculate min, median, mean and maximum landfast ice extents
    lfiw_m4.min(:,m) = min(lfiw_mo, [], 2, 'omitnan');
    lfiw_m4.median(:,m) = median(lfiw_mo, 2, 'omitnan');
    lfiw_m4.mean(:,m) = mean(lfiw_mo, 2, 'omitnan');
    lfiw_m4.max(:,m) = max(lfiw_mo, [], 2, 'omitnan');
        
    % Calculate upper and lower 10% landfast ice extents
    for v=1:Nv    
        fini = isfinite(lfiw_mo(v,:));
        lfiw_v = lfiw_mo(v,fini);
        lfiw_sort = sort(lfiw_v);
        Nfin = sum(fini);
        loi = max([1,floor(0.1*Nfin)]);
        upi = ceil(0.9*Nfin);
        lfiw_m4.low(v,m) = lfiw_sort(loi);
        lfiw_m4.upp(v,m) = lfiw_sort(upi);
        notalot = 0;
    end    
    
    
    % For each monthly metric:
    % - output SLIE lat/lons as csv file
    % - rasterize the landfast ice area
    mmmstack = lmask * 0;
    lmustack = lmask * 0;
    meanedge = lmask * 0;
    for i = 1:numel(m4)
        mmm = m4{i};
        
        %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Convert m5 lfiw data to a raster image
        slieimg = slie_lfiw2raster(lfiw_m4.(mmm)(:,m), cvecs, lmask, ...
                                   'pixsz', cfg.pixsz, ...
                                   'csegs', csegs, ...
                                   'shadregs', shadregs, ...
                                   'filtwid', 0);

        % Add to mmm stack if min, median, or max
        if any(strcmp({'min','median','max'}, mmm))
            mmmstack = mmmstack + uint8(slieimg > 128);
        end
        
        % Add to lmu stack if lower %, median, or upper 10%
        if any(strcmp({'low','median','upp'}, mmm))
            lmustack = lmustack + uint8(slieimg > 128);
        end


        %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Convert lfiw to pixel x,y coords
        [x,y] = slie_lfiw2xy(lfiw_m4.(mmm)(:,m), cvecs, cfg.pixsz, csegs);
        
        % Apply median filter to de-spike
        x = medianfilter(x, despike);
        y = medianfilter(y, despike);

        if strcmp('mean', mmm)
            [edgex, edgey] = jointhedots(round(x), round(y));
            edgei = sub2ind(cfg.sz, edgey, edgex);
            meanedge(edgei) = 1;
            meanedgefat = conv2(meanedge, skern, 'same')>0;
        end

        % Write data to files
        mostr = num2str(m, 'moslie_%02d');
        yrstr = num2str([yr0, mod(yr1,100)], '%04d-%02d');
        outstem = [mostr '_' mmm '_' yrstr];
        
        % Write vector lat/lon to ASCII CSV
        [mapx, mapy] = R.intrinsicToWorld(x,y);
        outfid = fopen([mapdir filesep outstem '_AkAlb.csv'], 'wt');
        fprintf(outfid, 'X,Y\n');
        fprintf(outfid, '%9.1f,%10.1f\n', [mapx, mapy]');
        fclose(outfid);
        [lat,lon] = projinv(gtinfo, mapx, mapy);
        outfid = fopen([lldir filesep outstem '_ll.csv'], 'wt');
        fprintf(outfid, 'Lat,Lon\n');
        fprintf(outfid, '%9.5f,%10.5f\n', [lat,lon]');
        fclose(outfid);
        
        % Write raster image to GeoTIFFF
        outtif = [outpath filesep outstem '.tif'];
        geotiffwrite(outtif, slieimg, R, 'GeoKeyDirectoryTag',geokeyDT);
        disp(['Written: ', outtif]);

        notalot = 0;
    end

    mmmstack(meanedgefat > 0) = 4;
    mmmstack(lmask == landval) = 5;
    mmmstack(lmask == oobval) = 6;
    mmmstack((lmask == shadowval) & (mmmstack == 0))= 7;
    outmmm = uint8(mmmstack);
    outtif = [outpath filesep mostr '_' yrstr '_MMM.tif'];
    geotiffwrite(outtif, outmmm, cmap, R, 'GeoKeyDirectoryTag',geokeyDT);
    disp(['Written: ', outtif]);
    notalot = 0;

    lmustack(meanedgefat > 0) = 4;
    lmustack(lmask == landval) = 5;
    lmustack(lmask == nodataval) = 6;
    lmustack((lmask == shadowval) & (lmustack == 0))= 7;
    outlmu = uint8(lmustack);
    outtif = [outpath filesep mostr '_' yrstr '_LMU.tif'];
    geotiffwrite(outtif, outlmu, cmap, R, 'GeoKeyDirectoryTag',geokeyDT);
    disp(['Written: ', outtif]);
    notalot = 0;

    
    
end

end


    
    
    
    
    