function slie_visualizebreakouts(dataset, varargin)
% A function to visualize each breakout listed the catalog
% produced by slie_findbreakouts.m
%
% USAGE
%      slie_visualizebreakouts
%      slie_visualizebreakouts(dataset)
%      slie_visualizebreakouts(__, 'cfgfile', cfgfile)
%      slie_visualizebreakouts(__, 'asadir', asadir)
%      slie_visualizebreakouts(__, 'bodir', outdir)
%      slie_visualizebreakouts(__, 'plotvis', plotvis)
%      slie_visualizebreakouts(__, 'overwrite', overwrite)
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
%     bodir : string specifying folder name inside the asadir where
%             all the breakout data files are stored
%             default: 'BreakOuts'
%   plotvis : string specifying visibility of figure created to visualize
%             timing and extent of breakout. Turning visibility off will
%             speed up the process. Valid values are 'on' or 'off'
%             default: 'off'
%
% OUTPUT:
%    
%
%
% NOTE: 
%   
% Andy Mahoney - October 2023
%
% ------------------------------------------------------------------------

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
        case 'plotvis'
            plotvis = varargin{a+1};
            a = a + 2;
        case 'overwrite'
            overwrite = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end   


% ------------------------------------------------------------------------
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

if ~exist('plotvis', 'var')
    plotvis = 'off';
end

if ~exist('overwrite', 'var')
    overwrite = 0;
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
sz = size(lmask);


% -----------------------------------------------------------------------
% Set up some things for plotting and output

% Create color table
bo_ctab = [[0.8,0.8,0.8]; ...   % No landfast ice
           [  1,  0,  0]; ...   % Retreat area 
           [  0,  0,  1]; ...   % Advance area
           [  1,  1,  1]; ...   % Landfast in both SLIEs
           [0.6,0.8,0.6]; ...   % Land
           [0.6,0.6,0.6]];      % Out of bounds
    

% Create output directory structure if necessary
outpath = [sfolder filesep asadir filesep bodir];
if ~exist(outpath, 'file')
    mkdir(outdir)
end
if ~exist([outpath filesep 'PNGs'], 'file')
    mkdir([outpath filesep 'PNGs']);
end
if ~exist([outpath filesep 'FIGs'], 'file')
    mkdir([outpath filesep 'FIGs']);
end


% -----------------------------------------------------------------------
% Start plotting figures for each breakout

% Read breakout catalog for this dataset
bopath = [sfolder filesep asadir filesep bodir];
bofile = [bopath filesep 'breakout_catalog_' dataset '.csv'];
bo = slie_readbreakoutcatalog(bofile);

% Set a dummy start season
currseason = 'yyyy-yy';

% Set limits for breakout map axes
xlim = [0,sz(2)]+0.5;
ylim = [0,sz(1)]+0.5;

% Create figure
% (make invisible to speed things up)
fig=figure ('units', 'inches', 'position', [4,4,12,8], ...
                               'paperposition', [0,0,12,8], ...
                               'visible', plotvis);

% Create axes for LFIW timeseries
ax1 = axes('parent', fig, 'pos', [0.05,0.35,0.25,0.3]);
ax1_th = title('', 'parent', ax1, 'interpreter', 'none');   % title handle

% Create axes for plotting breakout image
ax2 = axes('parent', fig, 'pos', [0.35,0.1,0.55,0.8], ...
           'xlim', xlim, 'ylim', ylim, 'ydir', 'reverse', ...
           'dataaspectratio', [1,1,1], ...
           'visible', 'off');
colormap(ax2, bo_ctab);
ax2_imh = image('parent', ax2, 'cdata', 1:6, ...
                'xdata', [1,sz(2)], 'ydata', [1,sz(1)]);
ax2_v1h = line(0,0, 'parent', ax2, 'color', 'cyan', 'linewidth', 2);
ax2_v2h = line(0,0, 'parent', ax2, 'color', 'magenta', 'linewidth', 2);


% Go through each breakout listed in catalog
Nbo = numel(bo.brkoutid);
boname = 'blank';
for b=1:Nbo
    
    % Count how many breakouts occurred in this SLIE so that we can
    % assign each breakout a unique filename
    sameslie = strcmp(bo.slie{b}, bo.slie);
    Nss = sum(sameslie);
    if Nss > 1
        nsuff = find(bo.brkoutid(sameslie) == bo.brkoutid(b));
        suffix = num2str(nsuff, '_%02d');
    else
        suffix = '';
    end

    % Now check to see if output files already exist for this breakout 
    outstem = ['BreakOut_' bo.slie{b}(1:13) suffix];
    outpng =  ['PNGs' filesep outstem '.png'];
    outfig =  ['FIGs' filesep outstem '.fig'];
    if exist([outpath filesep outfig], 'file') && (overwrite == 0)
        disp(['Outfile already exists for breakout ' outstem]);
        continue;
    end


    % If this breakout is from a different season than the last one: 
    % - read landfast ice width data for the season
    % - figures out dates of SLIEs 
    % - set x-limits and ticks for timeseries axes
    if ~strcmp(bo.season{b}, currseason)
        [lfiw, slienames] = slie_readsliedata([sfolder filesep ...
                                                 bo.season{b} filesep ...
                                                'lfiw.csv']);
        sdates = slie_name2date(slienames);
        [yr0, ~, ~] = datevec(sdates(1));
        [yr1, ~, ~] = datevec(sdates(end));
        tlim = [datenum(yr0,11,1),datenum(yr1,8,1)];
        tt = datenum([zeros(10,1)+yr0,(11:20)',ones(10,1),zeros(10,3)]);
        ttl = datestr(tt(1:end-1), 'mmm');
        ttm = (tt(1:end-1)+tt(2:end))/2;
        
        set(ax1, 'xlim', tlim, 'xtick', tt, 'xticklabels', []);
        text(ttm, zeros(1,9), ttl, 'parent', ax1, ...
                                   'fontsize', 14, ...
                                   'horizontalal', 'center', ...
                                   'verticalal', 'top');

        currseason = bo.season{b};
    end


    % Find the SLIE image associated with the breakout
    % and the one preceding it
    bosi = find(strcmp(slienames, bo.slie{b}));

    % Print some stuff for the user to see
    fprintf('%0.0f/%0.0f\n', b, Nbo);

    % Check to see if we need to load new SLIE image data

    if ~strcmp(boname, bo.slie{b})
        boname = bo.slie{b};
        preboname = slienames{bosi-1};

        % Read breakout and pre-breakout GeoTIFFs
        sspath = [sfolder filesep bo.season{b}];
        botif = geotiffread([sspath filesep boname]);
        pbotif = geotiffread([sspath filesep preboname]);

        fprintf('Read GeoTIFFs: %s, %s\n', boname, preboname);
    end

    % --------------------------------------------------------------------           
    % Plot LFIW timeseries for each vector in breakout
    v0 = bo.vec0(b);
    v1 = bo.vec1(b);
    bovi = v0:v1;
    Nbov = numel(bovi);
    lcols = generatecolorramp({'cya', 'mag'}, Nbov);
    ax1_lh1 = zeros(Nbov, 1);
    for v = 1:Nbov
        ax1_lh1(v) = line(sdates, lfiw(bovi(v), :), ...
                          'parent', ax1, ...
                          'color', lcols(v,:));
    end

    % Draw vertical red line indicating time of breakout
    lfiw_max = max(lfiw(v0:v1, :), [], 'all');
    sd = sdates(bosi);
    ax1_lh2 = line([sd,sd],[0,lfiw_max],'parent', ax1, ...
                                        'color', 'red', ...
                                        'linewidth', 2);

    % Set title with SLIE / date info
    tstr = [bo.slie{b} ': v=',num2str(v0,'%0d'), ...
                        ' - ',num2str(v1,'%0d')];
    set(ax1_th, 'string', tstr);


    % --------------------------------------------------- 
    % Plot SLIE images showing breakup region

    % Stack before and after images
    % to create visualization of breakout
    boimg = 1*(pbotif > 0) + 2*(botif > 0);
    boimg = (lmask==0).*boimg + 4*(lmask==255) + 5*(lmask==128);


    % Display combined before- and after-breakout image
    set(ax2_imh, 'cdata', boimg+1);

    % Add bounding coast vectors
    set(ax2_v1h, 'xdata', [cvecs(v0,1), cvecs(v0,3)], ...
                 'ydata', [cvecs(v0,2), cvecs(v0,4)]);
    set(ax2_v2h, 'xdata', [cvecs(v1,1), cvecs(v1,3)], ...
                 'ydata', [cvecs(v1,2), cvecs(v1,4)]);

    % --------------------------------------------------- 
    % Save figure as image and fig file
    print(fig, '-dpng', '-r300', [outpath filesep outpng]);
    savefig(fig, [outpath filesep outfig])
    
    notalot = 0;

    % Delete lines in timeseries plot
    delete(ax1_lh1);
    delete(ax1_lh2);
        
    
end


end