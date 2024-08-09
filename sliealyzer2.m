function sliealyzer2(dataset, varargin)
% Function to analyze SLIE geotiffs created from NIC ice charts
%
%
% USAGE
%   sliealyzer2
%   sliealyzer2(dataset)
%   sliealyzer2(__, 'cfgfile', cfgfile)
%   sliealyzer2(__, 'asadir', asadir)
%   sliealyzer2(__, 'skiplfiw', skiplfiw)
%   sliealyzer2(__, 'timeshift', timeshift)
%   sliealyzer2(__, 'doyspan', doyspan)
%
% INPUTS:
%   dataset : string specifyingplot( which dataset to use.
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
%  skiplfiw : a scalar flag indicating whether or not to skip calculation
%             of landfast ice width (lfiw) from the SLIE images
%             This saves time if recalculating derived data products
%             like key event dates, etc.
%             Default: 0 (don't skip lfiw calculation)
%   doyspan : a 2-element vector specifying the days of year (doy) between
%             which to run the sliealyzer. If the first value is less than
%             the second, it is assumed that the desired range wraps around
%             the end of the caledar year
%             Default: [274, 212] (Oct 1 - Jul 31)
% timeshift : a scalar flag indicating whether or not  to apply
%             "timeshifting" when calculating the date associated with
%             a particlar lfiw measurement. See slie_calcsliedates.m
%             for more details.
%             Default: 0 (do NOT apply timeshifting)
%             ** NOTE ** Timeshfting was applied in the analysis described
%                        by Mahoney et al. (2007; 2014), but innaccuracies
%                        have since been discovered and it is no longer
%                        recommended
%             
%
% OUTPUTS:
%   This code generates many comma-separated ASCII files, which are saved
%   either in the data folder for each season, or in the All-Seasons
%   Analyis folder (asadir).
%   These data files include:
%   - lfiw : landfast ice width measured along coast vectors
%            (see slie_calclfiw.m)
%   - wdep : water depth at intersection point with SLIE and coast vectors
%            (see slie_calclfiw.m)
%   - key events: dates of key events within landfast ice annual cycle
%                 (see slie_calcslieevents.m)
%
% NOTES:
%    - All the necessary information for finding the necessary input
%      datafiles and processing parameters are contained in the 
%      specifically-formated configuration file designated by configfile
%
%
% Andy Mahoney October 2023
%



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
        case 'skiplfiw'
            skiplfiw = varargin{a+1};
            a = a + 2;
        case 'doyrange'
            doyrange = varargin{a+1};
            a = a + 2;
        case 'timeshift'
            timeshift = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end   

% Set default values for parameters not specified on command line
if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

if ~exist('skiplfiw', 'var')
    skiplfiw = 0;
end

if ~exist('doyrange', 'var')
    doyrange = [274, 212];
end

if ~exist('timeshift', 'var')
    timeshift = 0;
end


% ----------------------------------------------------------------
% Set default dataset if not specified
if ~exist('dataset', 'var')
    dataset = 'beau';
end

% Identify the location where data are stored for specfied dataset
sfolder = slie_get_basefolder(dataset);

% Make directory for All Seasons Analysis if necessary
asapath = [sfolder filesep asadir];
if ~exist(asapath, 'file')
    mkdir(asapath);
end


% ----------------------------------------------------------------
% Read config file and other data files specified therein
cfg = slie_readconfigfile(sfolder, cfgfile);

% Read coast vector file
cvecs = slie_readcvecsfile(cfg.cvecfile);
Nv = numel(cvecs)/4;

% Read coastal distance file
cdist = geotiffread(cfg.cdistfile);

% Read bathymetry file
bathy = geotiffread(cfg.bathyfile);

% Read coast vector labels
[cvlabs, labs_i, Nlabs] = slie_readcveclabsfile(cfg.cveclabelsfile);
cveclabs = cell(Nv,1);
cveclabs(labs_i) = cvlabs;

Nseas = numel(cfg.seasons);

% Define regular expression for extracting dates from SLIE filenames
sdrexp = ['[' cfg.slieprefix '](?<yr>[1-2][0-9][0-9][0-9])' ...
          '_(?<d0>[0-9][0-9][0-9])-' ...
           '(?<d1>[0-9][0-9][0-9]).*'];


% ----------------------------------------------------------------
% Calculate landfast ice width (lfiw) and water depth at SLIE (wdep)
% for all SLIEs in each season
Nslie_all = 0;

for ss = 1:Nseas
    seasondir = [sfolder filesep cfg.seasons{ss}];


    % Find SLIE tif files for this folder
    % - taking account of the possibility for multiple file prefices
    %   specified in config file
    sp = cfg.slieprefix(1);
    sliefiles = dir([seasondir filesep sp '*slie.tif']);
    for p=2:numel(cfg.slieprefix)
        sp = cfg.slieprefix(p);
        sliefiles = [sliefiles; dir([seasondir filesep sp '*slie*.tif'])];
    end
    Nslie_in = numel(sliefiles);

    % Make a cell array of filenames and extract dates
    flist = cell(Nslie_in,1);
    mosdates = zeros(Nslie_in,3);
    for s=1:Nslie_in
        flist{s} = sliefiles(s).name;
        sd = regexp(flist{s}, sdrexp, 'names');
        %fprintf('%s: %s, %s, %s\n', flist{s}, sd.yr,sd.d0,sd.d1);
        syr = str2double(sd.yr);
        sd0 = str2double(sd.d0);
        sd1 = str2double(sd.d1);
        sd1 = sd1 + (sd1<sd0)*(365+1*(mod(syr,4)==0));
        mosdates(s,1) = datenum(syr,1,1) + sd0 -1;
        mosdates(s,3) = datenum(syr,1,1) + sd1 -1;
        mosdates(s,2) = mean(mosdates(s,[1,3]));
    end

    % Determine which files correspond to dates inside doyrange
    mosdoy = dayofyear(mosdates(:,2));
    if doyrange(1) < doyrange(2)
        keepi = (mosdoy >= doyrange(1)) & (mosdoy <= doyrange(2));
    else
        keepi = (mosdoy <= doyrange(2)) | (mosdoy >= doyrange(1));
    end
    sliefiles = sliefiles(keepi);
    flist = flist(keepi);
    Nslie = sum(keepi);    
    fprintf('%d SLIE found. Keeping %d between %s and %s.\n', ...
            Nslie_in, Nslie, datestr(doyrange(1), 'mmm dd'), ...
                             datestr(doyrange(2), 'mmm dd'));


    if ~skiplfiw             
        % Pre-assign arrays for holding data for this season
        lfiw_seas = zeros(Nv,Nslie);
        wdep_seas = zeros(Nv,Nslie);


        % Calculate landfast ice width and water depth at SLIE
        for s=1:Nslie
            disp(sliefiles(s).name);
            slieimg = geotiffread([seasondir filesep sliefiles(s).name]);

            [lfiw, wdep] = slie_calclfiw(slieimg, cvecs, cdist, 'bathy', bathy);

            %plot(lfiw);
            %line(1:numel(lfiw), wdep);

            lfiw_seas(:,s) = lfiw;
            wdep_seas(:,s) = wdep;
        end
        
        % Check for missing landfast ice polygons 
        % based on spurious, shortlived reductions in landast ice width
        mp = slie_findmissingpolygondata(lfiw_seas, cvecs, cfg.pixsz);
        lfiw_seas(mp==1) = NaN;
        wdep_seas(mp==1) = NaN;
        
        slie_writemisssingpolygondata(mp, seasondir, sliefiles);
        
        outname = [seasondir filesep 'lfiw.csv'];
        slie_writesliedata(outname, lfiw_seas, flist, ...
                           'veclabels', cveclabs, ...
                           'datafmt', '%5.1f');

        outname = [seasondir filesep 'wdep.csv'];
        slie_writesliedata(outname, wdep_seas, flist, ...
                           'veclabels', cveclabs, ...
                           'datafmt', '%4.0f');
                                      
        % Time-shift the SLIE dates if specified in function call
        if timeshift
            % Determine best-guess dates for SLIEs based on
            %  i) ice chart date ranges
            % ii) lfiw trends
            sliedates = slie_calcsliedates(lfiw_seas, mosdates(keepi,:));
        else
            % Just use mid-point of SLIE period
            sdates = mosdates(keepi,2);
            sliedates = repmat(sdates', Nv, 1);
        end

        % Write slie dates as date numbers
        outname = [seasondir filesep 'sliedatenum.csv'];
        slie_writesliedata(outname, sliedates, flist, ...
                           'veclabels', cveclabs, ...
                           'datafmt', '%0.0f');

        % Write slie dates as date strings
        sliedatesstr = cell(Nv,Nslie);
        for s=1:Nslie
            datechararray = datestr(sliedates(:,s), 'yyyy/mm/dd');
            sliedatesstr(:,s) = mat2cell(datechararray, ones(Nv,1), 10);
        end
        outname = [seasondir filesep 'sliedatestr.csv'];
        slie_writesliedata(outname, sliedatesstr, flist, ...
                           'veclabels', cveclabs, ...
                           'datafmt', '%s');
                              
                                      
    else
        % Otherwise, read lfiw, wdep and slie dates from previously generated file
        inname = [seasondir filesep 'lfiw.csv'];
        [lfiw_seas, flist] = slie_readsliedata(inname, '%f');
        inname = [seasondir filesep 'wdep.csv'];
        wdep_seas = slie_readsliedata(inname, '%f');
        inname = [seasondir filesep 'sliedatenum.csv'];
        sliedates = slie_readsliedata(inname, '%f');
        inname = [seasondir filesep 'sliedatestr.csv'];
        sliedatesstr = slie_readsliedata(inname, '%10s');
        Nslie = numel(flist);
        notalot = 0;
    end
    

    % Calculate dates of key fast ice events
    ev = slie_calcslieevents(lfiw_seas, wdep_seas, sliedates, ...
                             'min_ice', 0.5, ...
                             'stab_dep', 15);

    notalot = 0;

    if ss == 1
        flist_all = flist;
        lfiw_all = lfiw_seas;
        wdep_all = wdep_seas;
        sliedates_all = sliedates;
        sliedatesstr_all = sliedatesstr;
        ev_all = ev;
    else
        flist_all = [flist_all; flist];
        lfiw_all = [lfiw_all, lfiw_seas];
        wdep_all = [wdep_all, wdep_seas];
        sliedates_all = [sliedates_all, sliedates];
        sliedatesstr_all = [sliedatesstr_all,sliedatesstr];
        evf = fieldnames(ev);
        for e=1:numel(evf)
            ev_all.(evf{e}) = [ev_all.(evf{e}), ev.(evf{e})];
        end
    end

    notalot = 0;

    Nslie_all = Nslie_all + Nslie;
end

% Write concatentated data for all seasons to ASCII files
slie_writesliedata([asapath filesep 'lfiw_all.csv'], ...
               lfiw_all, flist_all, ...
               'veclabels', cveclabs, ...
               'datafmt', '%5.1f');
slie_writesliedata([asapath filesep 'wdep_all.csv'], ...
               wdep_all, flist_all, ...
               'veclabels', cveclabs, ...
               'datafmt', '%4.0f');
slie_writesliedata([asapath filesep 'sliedatenum_all.csv'], ...
               sliedates_all, flist_all, ...
               'veclabels', cveclabs, ...
               'datafmt', '%0.0f');
slie_writesliedata([asapath filesep 'sliedatestr_all.csv'], ...
               sliedatesstr_all, flist_all, ...
               'veclabels', cveclabs, ...
               'datafmt', '%s');
               
% Do same for each SLIE event
evdir = [asapath filesep 'SLIE_events'];
if ~exist(evdir, 'file')
    mkdir(evdir);
end
evf = fieldnames(ev);
for e=1:numel(evf)
    estr = num2str(e, 'event_%02d_');
    slie_writesliedata([evdir filesep estr evf{e} '_datenum.csv'], ...
                   ev_all.(evf{e}), cfg.seasons, ...
                   'veclabels', cveclabs, ...
                   'datafmt', '%0.0f');
        
    evdatestr =  cell(Nv,Nseas);
    for ss=1:Nseas
        evdate_ss = ev_all.(evf{e})(:,ss);
        fini = isfinite(evdate_ss);
        evdate_fin = evdate_ss(fini);
        Nfin = sum(fini);
        datecell = mat2cell(datestr(evdate_fin, 'mm/dd'), ones(Nfin,1), 5);
        evdatestr(fini,ss) = datecell;
    end
    slie_writesliedata([evdir filesep estr evf{e} '_datestr.csv'], ...
                   evdatestr, cfg.seasons, ...
                   'veclabels', cveclabs, ...
                   'datafmt', '%5s');

end          
           


notalot = 0;
end