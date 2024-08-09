function slie_catalogbreakouts(dataset, varargin)
% A function to identify breakouts affecting broad sections of coast
% according to quantitative criteria specified below and write a list
% of these breakouts to a comma-separated ASCII file
%
%
% USAGE
%        slie_catalogbreakouts
%        slie_catalogbreakouts(dataset)
%        slie_catalogbreakouts(__, 'cfigfile', cfgfile)
%        slie_catalogbreakouts(__, 'asadir', asadir)
%        slie_catalogbreakouts(__, 'outdir', outdir)
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
%    outdir : string specifying folder name inside the asadir where
%             all the output files will be saved
%             default: 'BreakOuts'
%
% OUTPUT:
%     
%   A breakout catalog
%      An ASCII file containing a list of each breakout specifying where
%      and when it occurred and how extensive it was
%
%
% Andy Mahoney - October 2023


% =================== DEFINE BREAKOUT CRITERIA ========================
% It's only a breakout if:
% 1. the remaining ice width is less than a certain fraction of the 
%    modal ice width at a given location, specified by rem_width 
% 2. the ice retreats by a distance greater than a certain fraction 
%    of the median ice width at a given location, specified by min_retr
% 3. it extends for more than a minimum number of consecutive coast 
%    vectors specified by min_extent
% 4. it occurs during specified months
rem_width = 0.5; 
min_retr = 0.33; 
min_extent = 25; 
bo_mon = [1,2,3,4];
% ===================  ======================== ========================
        


% ----------------------------------------------------------------
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
        case 'outdir'
            outdir = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end   


% Set default values for parameters not specified on command line
if ~exist('dataset', 'var')
    dataset = 'beau';
end

if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

if ~exist('outdir', 'var')
    outdir = 'BreakOuts';
end

% -----------------------------------------------------------------------
% Identify the location where data are stored for specfied dataset
sfolder = slie_get_basefolder(dataset);


% Read config file and other data files specified therein
cfg = slie_readconfigfile(sfolder, cfgfile);


% Read coast vector file
[cvecs, Nvecs] = slie_readcvecsfile(cfg.cvecfile);

% Read landmask
[lmask, R] = geotiffread(cfg.lmaskfile);
gtinfo = geotiffinfo(cfg.lmaskfile);
mstruct = geotiff2mstruct(gtinfo);

% Read all-season lfiw data
[lfiw, slienames] = slie_readsliedata([sfolder filesep asadir filesep ...
                                        'lfiw_all.csv']);
Nslie = numel(slienames);                                   
                                   

% Calcuate all-season median lfiw for each vector
lfiw_med = median(lfiw, 2, 'omitnan');


% Read SLIE dates from file
sdates2 = slie_readsliedata([sfolder filesep asadir filesep ...
                            'sliedatenum_all.csv']);
sdates = sdates2(1,:);

% Initialize breakout counter
Nbo = 0;

% Check to see if output folder exists
outpath = [sfolder filesep asadir filesep outdir];
if ~exist(outpath, 'file')
    mkdir(outpath);
end

% Open output file for writing data to text file
outfile = [outpath filesep 'breakout_catalog_' dataset '.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, ['BrkOutID,Season,SLIE,SLIEdate,vec0,vec1,' ...
                 'meanwid_km,coastlen_km,area_km2,Latitude,Longitude\n']);
outfstr = '%4.0f,%s,%s,%s,%4.0f,%4.0f,%5.1f,%6.1f,%8.1f,%10.6f,%11.6f\n';



% Go through LFIW data by season
Nseasons = numel(cfg.seasons);
for ss = 1:Nseasons

    % Identify lfiw data that falls within this season
    yr0 = str2double(cfg.seasons{ss}(1:4));
    yr1 = yr0+1;
    sep01 = datenum(yr0,9,1);
    aug31 = datenum(yr1,8,31);
    ssi = (sdates >= sep01) & (sdates <= aug31);
    Nslie = sum(ssi);

    % Extract lfiw data for this season
    lfiw_s = lfiw(:,ssi);
    sd = sdates(ssi);
    sn = slienames(ssi);

    % Calculate retreat between consective SLIEs, as a positive number
    % (cannot do this for first SLIE)
    lfiw_r = lfiw_s*0;
    lfiw_r(:,2:end,:) = -diff(lfiw_s, 1, 2);

    % Set retreat to zero for any advancing SLIEs
    lfiw_r = lfiw_r .* (lfiw_r > 0);

    % Normalize both width and retreat by median lfiw
    lfiw_m2 = repmat(lfiw_med,[1,Nslie]);
    lfiw_n = lfiw_s./lfiw_m2;
    lfiw_nr  = lfiw_r./lfiw_m2;


    % Calculate normalized ice retreat
    % ** impossible to identify breakout in first SLIE of a season
    
    % Threshold ice retreat and width remaining to identify breakouts
    lfiw_bo = 1*((lfiw_nr > min_retr) & (lfiw_n < rem_width));

    % Now go through each SLIE, identifying any breakouts
    for s=1:Nslie

        % Figure out which month we're in
        [~,mo] = datevec(sd(s));
        
        % Skip breakouts not within specified months
        if ~any(bo_mon == mo)
            continue;
        end

        % Extract breakout state for this SLIE
        lfiw_bo_s = squeeze(lfiw_bo(:,s));

        % Take along-coast differential to find locations
        % where breakout state changes between retreat and not-retreat
        bodiff = diff(lfiw_bo_s);
        i0 = find([lfiw_bo_s(1); bodiff] == 1);
        i1 = find([bodiff; -lfiw_bo_s(end)] == -1);
        nboseg = numel(i0);
        if nboseg ~= numel(i1)
            wtf = 99;
        end


        % Now go through each segment and expand on either side to include
        % all vectors where ice retreated but didn't necessarily meet
        % the breakout criteria
        boseg = zeros(nboseg,2);
        for b=1:nboseg
            seg0 = find(lfiw_r(1:i0(b),s) == 0, 1, 'last')+1;
            seg1 = i1(b) + find(lfiw_r(i1(b):end,s) == 0, 1, 'first')-2;
            if ~any(seg0), seg0 = 1; end
            if ~any(seg1), seg1 = Nvecs; end            
            boseg(b,:) = [seg0, seg1];
        end

        % Merge any that are just separated by a few vectors
        % - simplest way to do this is to make more duplicates
        %   which are removed when we select unique breakouts
        distbetween = boseg(2:end,1) -  boseg(1:end-1,2);
        closeenough = find(distbetween < min_extent/5);
        Nce = numel(closeenough);
        for c=1:Nce
            cc = closeenough(c);
            boseg(cc,2) = boseg(cc+1,2);
            boseg(cc+1,1) = boseg(cc,1);
        end

        % Remove any duplicates create when adjacent segments are merged
        boseg = unique(boseg, 'rows');

        
        % Exclude any segments that are too short
        % (i.e. occur across too few coast vectors)
        tooshort = diff(boseg, 1, 2) <=  min_extent;
        boseg = boseg(~tooshort, :);

        % See how many breakouts we're left with
        Nboseg = size(boseg, 1);

        % Calculate properties of each breakout 
        % and write info to catalog
        for b=1:Nboseg
            v0 = boseg(b,1);
            v1 = boseg(b,2);

            % Calculate length, mean width, and area of breakout
            boW = mean(lfiw_r(v0:v1, s), 'omitnan');
            vxdif = (cvecs(v0+1:v1,1)-cvecs(v0:v1-1,1));
            vydif = (cvecs(v0+1:v1,2)-cvecs(v0:v1-1,2));
            boL = sum(sqrt(vxdif.^2 + vydif.^2))*cfg.pixsz;

            bow0 = lfiw_s(v0:v1, s-1);
            [box0, boy0] = w2xy(bow0, cvecs(v0:v1, :), cfg.pixsz);
            bow1 = lfiw_s(v0:v1, s);
            [box1, boy1] = w2xy(bow1, cvecs(v0:v1, :), cfg.pixsz);
            bopolyx = [box0; box1(end:-1:1)];
            bopolyy = [boy0; boy1(end:-1:1)];
            boA = polyarea(bopolyx, bopolyy)*(cfg.pixsz^2);
            bopoly = polyshape(bopolyx, bopolyy);
            [bocx, bocy] = centroid(bopoly);
            [boCx, boCy] = R.intrinsicToWorld(bocx,bocy);
            [bolat,bolon] = projinv(mstruct, boCx, boCy);

            % Write info for this breakout to file
            Nbo = Nbo + 1;
            fprintf(outfid, outfstr, ...
                            Nbo, ...
                            cfg.seasons{ss}, ...
                            sn{s}, ...
                            datestr(sd(s), 'yyyy-mm-dd'), ...
                            v0, v1, boW, boL, boA, bolat, bolon);
           
            notalot = 0;
        end
    end
end

% Close catalog file
fclose(outfid);

fprintf('%0.0f breakouts found.', Nbo);
disp(['Written: ' outfile]);


end



function [x,y] = w2xy(w,cvecs,pixsz)
% Short single-purpose function to find the pixel x-y coordinates
% corresponding to width measurements along particular coast vecs
%
% Kinda like slie_lfiw2xy, but without all the coast segementation stuff

dx = cvecs(:,3) - cvecs(:,1);
dy = cvecs(:,4) - cvecs(:,2);
az = atan2(dy,dx);   % Counterclockwise from x-axis

pw = w/pixsz;
x = cvecs(:,1) + pw.*cos(az);
y = cvecs(:,2) + pw.*sin(az);

% Diagnostics
%fig = figure;
%ax = axes('parent', fig);
%line(cvecs(:,1), cvecs(:,2), 'parent', ax, 'color', 'black');
%for v=1:numel(w)
%    line(cvecs(v,[1,3]), cvecs(v,[2,4]), 'parent', ax, ...
%                                         'color', [0.8,0.8,0.8]);
%end
%line(x,y, 'parent', ax, 'color', 'blue');
%notalot = 0;

end










