function slie_dailify(dataset, varargin)
% Function to interpolate slie data to a daily resolution by applying 
% nearest neighbor interpolation to the derived landfast ice width values
% A running minimum can be applied and results can be converted
% to binarized rasters
%
% USAGE:
%   slie_dailify
%   slie_dailify(dataset)
%   slie_dailify(__, 'cfgfile', cfgfile)
%   slie_dailify(__, 'asadir', asadir)
%   slie_dailify(__, 'dailydir', dailydir)
%   slie_dailify(__, 'runmin', runmin)
%   slie_dailify(__, 'skiptiff', skiptiff)
%
%
% INPUTS:
%   dataset : string specifying which dataset to use.
%             See slie_get_basefolder.m for accepted values.
%             default: 'beau'
%   cfgfile : string specifying the name of the configuraion file stored
%             in the region folder containing essential information for
%             performing all analysis of SLIE data. 
%             See slie_readconfigfile.m for more details
%             default: 'slieconfig.txt'
%    asadir : string specifying folder name in region folder where
%             the data files from the all-seasons analysis will be saved
%             default: 'AllSeasonsAnalysis2023'
%  dailydir : string specifying folder name in where the daily-resolution
%             data files will be saved
%             default: A folder with the same name as that where data for 
%                      input are stored, but with '_daily' appended
%    runmin : integer specifying length of window in days for applying
%             a running minimum
%             default: 1 (no running minimum applied)
%  skiptiff : a scalar flag specifying whether or not to skip creation of
%             binarized raster GeoTIFFs 
%             Default: 0 (do not skip tiff creation)
%
% Andy Mahoney October 2023
%              January 2024: amended handling of NaNs to avoid filling
%                            in trailing NaN values
% ----------------------------------------------------------------


% - - - - - - - - - - - - - - - - - -
% Check for variable arguments
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'cfgfile'
            cfgfile = varargin{a+1};
            a = a + 2;
        case 'asadir'
            asadir = varargin{a+1};
            a = a + 2;
        case 'runmin'
            runmin = varargin{a+1};
            a = a + 2;
        case 'skiptiff'
            skiptiff = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end   

% ----------------------------------------------------------------
% Set default dataset if not specified
if ~exist('dataset', 'var')
    dataset = 'beau';
end

% Identify the location where data are stored for specfied dataset
sfolder = slie_get_basefolder(dataset);

% Set other default values for parameters not specified on command line
if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end
asapath = [sfolder filesep asadir];

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

if ~exist('runmin', 'var')
    runmin = 1;
end

if runmin > 1
    runminsuff = num2str(runmin, '_rm%0.0f');
else
    runminsuff = '';
end

if ~exist('dailydir', 'var')
    dailydir = [sfolder '_daily' runminsuff];
end


if ~exist('skiptiff', 'var')
    skiptiff = 0;
end


% ----------------------------------------------------------------
% Make folder for Daily SLIEs if necessary
dailydir = [sfolder '_Daily'];
if ~exist(dailydir, 'file')
    mkdir(dailydir);
end



% ----------------------------------------------------------------
% Read config file and other data files specified therein
cfg = slie_readconfigfile(sfolder, cfgfile);

% Read coast vector file
[cvecs, Nv] = slie_readcvecsfile(cfg.cvecfile);

% Read coast vector labels
[cvlabs, labs_i, Nlabs] = slie_readcveclabsfile(cfg.cveclabelsfile);
cveclabs = cell(Nv,1);
cveclabs(labs_i) = cvlabs;

% Read coastal segment data
csegs = slie_read_coastalsegmentsfile(cfg.coastsegmentfile);

       
% Read land mask and get GeoTIFF information
[lmask, R] = geotiffread(cfg.lmaskfile);
gtinfo = geotiffinfo(cfg.lmaskfile);
geokeyDT = gtinfo.GeoTIFFTags.GeoKeyDirectoryTag;

% Remove any "unknown" fields from geokey
if isfield(geokeyDT, 'Unknown')
    geokeyDT = rmfield(geokeyDT, 'Unknown');
end
sz = size(lmask);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Identify individual shadow regions in land mask objects
% and indentify the non-shadow pixels adjacent to them
if skiptiff == 0
    shadows = lmask == 64;
    shad = regionprops(shadows, 'boundingbox', 'pixelidxlist', 'image');
    Nshad = numel(shad);
    %shadadjpix = cell(Nshad);
    for sh = 1:Nshad
        
        % Place feature image into large image 
        % with extra pixel all the way around
        shadimgX = zeros(size(shad(sh).Image)+2);
        shadimgX(2:end-1,2:end-1) = shad(sh).Image;
        
        % Grow object image with a 3x3 convolution
        shadconv = conv2(shadimgX, ones(3), 'same');
        
        % Identify the outer perimeter of feature
        outerperim = (shadconv > 0) & (shadimgX == 0);
        opi = find(outerperim);
        [opy,opx] = ind2sub(size(shadconv), opi);
        
        % Translate sub-image pixel index values into whole-image values
        % checking for oub-of-bounds (oob) pixels
        x0 = floor(shad(sh).BoundingBox(1));
        y0 = floor(shad(sh).BoundingBox(2));
        opxx = opx+x0-1;
        opyy = opy+y0-1;
        oob = (opxx < 1) | (opxx > sz(2)) | (opyy < 1) | opyy > sz(1);
        opxx = opxx(~oob);
        opyy = opyy(~oob);
        shad(sh).adjpix = sub2ind(sz, opyy, opxx);
        notalot = 0;
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set days of year defining beginning and end of each season
doy0_mmdd = [10, 01];
doy1_mmdd = [7, 31];

% Count number of seasons
Nseas = numel(cfg.seasons);
for ss = 1:Nseas
    
    % Create season folder in Daily SLIE directory if necessary
    sddir = [dailydir filesep cfg.seasons{ss}];
    if ~exist(sddir, 'file')
        mkdir(sddir);
    end
    
    % Read LFIW and water depth data for this season
    lfiwfile = [sfolder filesep cfg.seasons{ss} filesep 'lfiw.csv'];
    [lfiw, slienames] = read_sliedatafile(lfiwfile, '%f');
    [sdates, spref] = slie_name2date(slienames);
    wdepfile = [sfolder filesep cfg.seasons{ss} filesep 'wdep.csv'];
    wdep = read_sliedatafile(wdepfile, '%f');
    
    % Fill-in any NaNs using the nearest non-NaN observation in time
    lfiw_nonan = fillmissing(lfiw, 'nearest', 2, ...
                                   'samplepoints', sdates', ...
                                   'maxgap', 1, ...
                                   'endvalues', 'none');
    wdep_nonan = fillmissing(wdep, 'nearest', 2, ...
                                   'samplepoints', sdates', ...
                                   'maxgap', 1, ...
                                   'endvalues', 'none');

    % Count SLIEs this season
    Nslie = size(lfiw, 2);
    
    % Create daily dates for this season
    yr0 = str2double(cfg.seasons{ss}(1:4));

    % Clip data to specified beginning / end doys if necessary
    d0 = max(datenum([yr0,doy0_mmdd]), min(sdates));
    d1 = min(datenum([yr0+1,doy1_mmdd]), max(sdates));
    sd = (d0:d1)';
    Nd = numel(sd);
    
    % Determine which SLIE date is nearest to each daily date
    sdi = interp1(sdates, 1:Nslie, sd, 'nearest');
    
    % Create SLIE filenames for each day
    dailynames = cell(Nd,1);
    for d=1:Nd
        dailynames{d} = [spref{sdi(d)} datestr(sd(d), 'yyyymmdd'), ...
                        '_dailyslie' runminsuff '.tif'];
    end

    % Interpolate lfiw data, apply running minimum and write to file
    lfiwDaily = interp1(sdates, lfiw_nonan', sd, 'nearest')';
    lfiwDaily_rmin = runningminimum(lfiwDaily, runmin);
    lfiwfile = [sddir filesep 'lfiw.csv'];
    slie_writesliedata(lfiwfile, lfiwDaily_rmin, dailynames, ...
                       'datafmt', '%4.1f', 'veclabels', cveclabs);

    % Interpolate wdep data, apply running minimum and write to file
    wdepDaily = interp1(sdates, wdep_nonan', sd, 'nearest')';
    wdepDaily_rmin = runningminimum(wdepDaily, runmin);
    wdepfile = [sddir filesep 'wdep.csv'];
    sdstr = mat2cell(datestr(sd, 'yyyy-mm-dd'), ones(Nd,1));
    slie_writesliedata(wdepfile, wdepDaily_rmin, dailynames, ...
                       'datafmt', '%0.0f', 'veclabels', cveclabs);


    % Write daily slie dates to file
    sdnumfile = [sddir filesep 'sliedatenum.csv'];
    sd2 = repmat(sd', Nv, 1);
    slie_writesliedata(sdnumfile, sd2, dailynames, ...
                       'datafmt', '%0.0f', 'veclabels', cveclabs);
    sdstrfile = [sddir filesep 'sliedatestr.csv'];
    sdstr2 = repmat(sdstr', Nv, 1);
    slie_writesliedata(sdstrfile, sdstr2, dailynames, ...
                       'datafmt', '%s', 'veclabels', cveclabs);
    

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    % Skip the follow GeoTIFF creation steps if specified in function call
    if skiptiff == 1
        continue;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Create dailified GeoTIFFs 
    for d=1:Nd

        slieimg = slie_lfiw2raster(lfiwDaily_rmin(:,d), cvecs, lmask, ...
                                'pixsz', cfg.pixsz, ...
                                'csegs', csegs, ...
                                'shadregs', shad);

        geotiffwrite([sddir filesep dailynames{d}], slieimg, R, ...
                     'GeoKeyDirectoryTag',geokeyDT);
        disp(['Written: ', dailynames{d}]);
        notalot = 0;

    end    
        
    notalot = 0;

end

end
       


function out = runningminimum(in, wid)
% Function to calculate a runnning minimum 
% along the 2nd dimension of a 2-D array

% USAGE
%  out = runningminimum(in, wid)
%
% INPUTS
%       in : a m-by-n array of values
%      wid : an odd integer specifying width in the n-dimension over which
%            to calculate the running minimum. If specified value is even
%            then a value of wid-1 is used
%
% OUTPUT
%      out : a m-by-n array of values with the running minimum applied
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Check dimensions of input data
sz = size(in);

% Ensure running window width is odd
wid = wid - 1*(mod(wid,2) == 0);
halfwid = floor(wid/2);

% Check to see if half-width is > 0
if halfwid > 0

    % Pre-assign output array
    out = zeros(sz);
    
    % Calculate running minimum for every column of input array
    for x = 1:sz(2)
        win0 = max(1,x-halfwid);
        win1 = min(sz(2),x+halfwid);
        win = win0:win1;
    
        out(:,x) = min(in(:,win), [], 2);
    end
else
    out = in;
end

end





       