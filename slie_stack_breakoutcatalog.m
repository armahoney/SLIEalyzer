function slie_stack_breakoutcatalog(dataset, varargin)
% Function to create stacked image of breakouts 
%
% USAGE
%   slie_stack_breakoutcatalog
%   slie_stack_breakoutcatalog(dataset)
%   slie_stack_breakoutcatalog(__, 'cfgfile', cfgfile)
%   slie_stack_breakoutcatalog(__, 'asadir', asadir)
%   slie_stack_breakoutcatalog(__, 'bodir', outdir)
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
%     bodir : string specifying folder name inside the asadir where
%             all the breakout data files are stored
%             default: 'BreakOuts'
%
% OUTPUT:
%    A geotiff stored within the output folder representing the 
%    cumulative number of breakouts affecting each pixel in the study
%    area of the specified data set
%
%
% NOTE: 
%   This stacking process is largely the same as what is already performed 
%   during the creation of the breakout catalog by slie_find_breakouts.m.
%   However, this function recalculates the stacked image from the catalog
%   without having to find each breakout again and may have small
%   differences such as division by the number of seasons to calculuate
%   breakout occurence as a mean annual count
%
% Andy Mahoney - October 2023
%
% ------------------------------------------------------------------------



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
        case 'bodir'
            bodir = varargin{a+1};
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

if ~exist('bodir', 'var')
    bodir = 'BreakOuts';
end

% -----------------------------------------------------------------------
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

% Read config file and other data files specified therein
cfg = slie_readconfigfile(sfolder, cfgfile);

% Count number of seasons in this dataset
Nseasons = numel(cfg.seasons);

% Read coast vector file
[cvecs, Nvecs] = slie_readcvecsfile(cfg.cvecfile);

% Read landmask
lmask = geotiffread(cfg.lmaskfile);

% Read breakout catalog for this dataset
bopath = [sfolder filesep asadir filesep bodir];
bofile = [bopath filesep 'breakout_catalog_' dataset '.csv'];
bo = slie_readbreakoutcatalog(bofile);

% Go through each breakout listed in catalog
Nbo = numel(bo.brkoutid);
for b=1:Nbo

    % List all SLIE files in season this breakout occured in
    sspath = [sfolder filesep bo.season{b}];
    sslist = filedigger(sspath, '*slie.tif');
    
    % Sort list by SLIE date
    % (accounts for instances in 2008 when Envisat-derived SLIEs
    %  appear out of chronological order when sorted alphabetically)
    [~, sorti] = sort(slie_name2date(sslist));
    sslist = sslist(sorti);
    
    % Find the SLIE image associated with the breakout
    % and the one preceding it
    boi = find(strcmp(sslist, bo.slie{b}));
    preboname = sslist{boi-1};

    % Print some stuff for the user to see
    fprintf('%0.0f/%0.0f: %s\n', b, Nbo, bo.slie{b});

    % Read breakout and pre-breakout GeoTIFFs
    gtinfo = geotiffinfo([sspath filesep bo.slie{b}]);
    [botif, R] = geotiffread([sspath filesep bo.slie{b}]);
    pbotif = geotiffread([sspath filesep preboname]);
    sz = size(botif);

    % If this if first breakout in catalog, 
    % initialize array for counting breakout occurrence 
    if b == 1
        bostack = zeros(sz);
    end

    % Identify coast vectors bounding breakout from catalog
    % - and use them to construct polygon around breakout
    v0 = bo.vec0(b);
    v1 = bo.vec1(b);
    bopolyx = [cvecs(v0,1); ...
               cvecs(v0,3); ...
               cvecs(v1,3); ...
               cvecs(v1,1); ...
               cvecs(v1:-1:v0,1)];
    bopolyy = [cvecs(v0,2); ...
               cvecs(v0,4); ...
               cvecs(v1,4); ...
               cvecs(v1,2); ...
               cvecs(v1:-1:v0,2)];

    % Subtract tiffs to identify pixels where landfast ice retreated
    bopix = (pbotif > botif);

    % Identify discrete features within this breakout image
    bofeatures = regionprops(bopix, 'Centroid', 'PixelIdxList');
    boCxy = vertcat(bofeatures(:).Centroid);
    boCx = boCxy(:,1);
    boCy = boCxy(:,2);

    % Figure out which breakout feature(s) lies within the bounded region
    inbo = find(inpolygon(boCx, boCy, bopolyx, bopolyy));
    
    % Increment corresponding breakout stack pixels
    for f=1:numel(inbo)
        boi = bofeatures(inbo(f)).PixelIdxList;
        bostack(boi) = bostack(boi) + 1;
    end

    % Some diagnostics
    %diagnimg = double(pbotif>0) + double(botif>0) + ...
    %           double(inbo) + 3*double(bopix);
    %imagesc(diagnimg);
    notalot = 0;

end

% Divide by number of seasons to derive mean annual count
bostack = bostack / Nseasons;

% Pull GeoTIFF information
geokeyDT = gtinfo.GeoTIFFTags.GeoKeyDirectoryTag;
% Remove any "unknown" fields from geokey
if isfield(geokeyDT, 'Unknown')
    geokeyDT = rmfield(geokeyDT, 'Unknown');
end

% Write breakout stack image to GeoTIFF
outtif = [bopath filesep 'BreakOut_stack_' dataset '.tif'];
geotiffwrite(outtif, bostack, R, 'GeoKeyDirectoryTag',geokeyDT);

disp(['Written: ', outtif]);

end
