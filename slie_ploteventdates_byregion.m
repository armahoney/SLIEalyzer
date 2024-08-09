function slie_ploteventdates_byregion(varargin)
% Function to plot the mean dates of key events 
% in the landfast ice cycle for each sub-region
%
% USAGE
%      slie_ploteventdates_byregion
%      slie_ploteventdates_byregion('dataset', dataset)
%      slie_ploteventdates_byregion('asadir', asadir)
%      slie_ploteventdates_byregion('cfgfile', cfgfile)
%      slie_ploteventdates_byregion('evfolder', evfolder)
%      slie_ploteventdates_byregion('frzthaw', frzthaw)
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
%  evfolder : string specifying folder name inside the asadir where
%             the output files will be saved
%             default: 'SLIE_Events'
%   frzthaw : a binary flag specifying whether or not to plot dates of
%             freezing an thawing onset, based on air temperature data.
%             ** NOTE **
%                  This requires the existence of a Matlab data file
%                  named 'cp_freezethaw.mat' inside the folder specified
%                  by asadir. This datafile is contains freeze/thaw data
%                  for every coast point (cp) and is created by the 
%                  the function slie_calc_sfctemp_at_cps.m base on NCEP
%                  surface air temperature data
% OUTPUTS
%       fig : A figure handle for the resulting plot 
%
%
% Andy Mahoney - October 2023

% ------------------------------------------------------------------------
% Handle any variable arguments specified with function call
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'dataset'
            dataset = varargin{a+1};
            a = a + 2;
        case 'cfgfile'
            cfgfile = varargin{a+1};
            a = a + 2;
        case 'asadir'
            asadir = varargin{a+1};
            a = a + 2;
        case 'evfolder'
            dataset = varargin{a+1};
            a = a + 2;
        case 'frzthaw'
            frzthaw = 1;
            a = a + 1;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end    
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set default values for parameters not specified in function call
if ~exist('dataset', 'var')
    dataset = 'beau';
end
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

% Plot freeze-thaw dates unless told not to do so
if ~exist('frzthaw', 'var'), frzthaw = 1; end

% Set default events data folder if not spec'd in function call
if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end
if ~exist('evfolder')
    evfolder = 'SLIE_Events';
end
evfolder = [sfolder filesep asadir filesep evfolder];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Read the SLIE configuration file to get subregion and vector label info
cfg = slie_readconfigfile(sfolder, cfgfile);



% read the coast vector subregions file
[cvecregs, regnames, regbounds, Nregs] = ...
                   slie_readcvecregsfile(cfg.cvecregionsfile);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Read in the different event data files
slieevent_strs = {'firstice', ...
                  'stableice', ...
                  'maxice', ...
                  'unstableice', ...
                  'breakup', ...
                  'icefree'};
Nevents = numel(slieevent_strs);

slieevent_labs = {{'First','Ice'}, ...
                  {'Stable','Ice'}, ...
                  {'Max','Ice'}, ...
                  {'Unstable','Ice'}, ...
                  {'Break','-up'}, ...
                  {'Ice','Free'}};


for e = 1:Nevents
    evfname = ['event_',num2str(e,'%02d_'), ...
               slieevent_strs{e},'_datenum.csv'];
    evfile = [evfolder  filesep  evfname];
              
    disp(['Reading ' evfname ' ...']);
    [evdata, Nvecs, Nseasons] = slie_readeventdates(evfile);
    
    % Use size data to initialize data array for first event
    if e == 1, evdates = zeros(Nvecs, Nseasons, Nevents); end
    
    evdates(:,:,e) = evdata;
end
    


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Calculate regional events dates for each subregion
% - mean
% - upper and lower quartiles
% ** Standard deviation across all years and vectors within a region
%    is too large to useful plot 
evdates_regmean = zeros(Nregs,Nevents) + NaN;     % Initialize arrays
evdates_reglow = zeros(Nregs,Nevents) + NaN;    % filled with NaNs
evdates_regupp = zeros(Nregs,Nevents) + NaN;    % filled with NaNs
for r = 1:Nregs
    ri = cvecregs == r;
    evdates_r = evdates(ri,:,:);
    for e = 1:Nevents
        evdates_e = evdates_r(:,:,e);
        evdates_sorted = sort(evdates_e(:));
        fin_i = isfinite(evdates_e);
        if any(fin_i(:))
            evdates_regmean(r,e) = mean(evdates_e(fin_i));
            lowi = floor(0.25*sum(fin_i(:)));
            uppi = floor(0.75*sum(fin_i(:)));
            evdates_reglow(r,e) = evdates_sorted(lowi);
            evdates_regupp(r,e) = evdates_sorted(uppi);            
        else
            notmuch=0;
        end
    end
end


% ========================= FIGURE PLOTTING ===============================                          
% Plot region-wide mean events dates for each region
% - x-axis represents day of year
% - sub-regions are spaced discretely along y-axis 

fig_wd = 17.5;
fig_ht = 9.5;
fig2 = figure('units', 'inches', ...
              'position', [1,5,fig_wd,fig_ht], ...
              'paperposition', [0,0,fig_wd,fig_ht]);

% Create x-tickmarks for first day of each month
x_yrs = [ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1];
x_mos = [ 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8];
x_das = zeros(1,12)+1;
xtick = datenum(x_yrs, x_mos, x_das);
xticklab = datestr(xtick, 'mmm dd');

% Set x-axis limits
xlim = [datenum(0,10,1),datenum(1,7,31)];


% Set y-axis limts and ticks, etc
ph = 0.4;  % The height of the patch for each region (in y-units)
ylim = [1-1.2*ph,Nregs+2*ph];

% Create regions labels by swap the hyphen for "to"
% and inserting newlines before and after it
reglabs = regexprep(regnames, '([a-zA-Z\s]) - ([a-zA-Z\s])', '$1\nto\n$2');


% Set colors for patches
growcol = [0.5,0.8,0.8];    % Muted cyan
retrcol = [0.5,0.6,0.8];    % Gray-blue
brupcol = [0.8,0.5,0.8];    % Muted magenta
noicecol = [0.9,1.0,1.0];   % Very pale blue
evvarcol = [0.7,0.7,0.7];   % Light gray

% Create axes for holding plot
lefmar = 0.2;
rigmar = 0.05;
botmar = 0.1;
topmar = 0.05;
pwd = 1 - (lefmar + rigmar);
pht = 1 - (botmar + topmar);
ax2 = axes('parent', fig2, 'pos', [lefmar,botmar,pwd,pht], ...
          'xlim', xlim, 'xtick', xtick, 'xticklabel', xticklab, ...
          'xgrid', 'on', 'fontsize', 18, ...
          'ytick', 1:Nregs, 'yticklabel', [], ...
          'ylim', ylim);
        

% Add y-tick labels as text objects
% (this allows more flexibility than trying to assign labels
%  via the axes definition)
ytl_x = zeros(Nregs,1) + xlim(1) - 0.02*diff(xlim);
ytl_y = 1:Nregs;
text(ytl_x, ytl_y, reglabs, 'parent', ax2, ...
                            'horizontalal', 'right', ...
                            'verticalal', 'middle', ...
                            'fontsize', 18, 'fontweight', 'bold');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
% Create patches for periods between key events
all_py = zeros(4, Nregs);
pre_ice_px = zeros(4, Nregs);
icegrow_px = zeros(4, Nregs);
iceretr_px = zeros(4, Nregs);
icebrup_px = zeros(4, Nregs);
postice_px = zeros(4, Nregs);
for r=1:Nregs
    
    % Set y-coords for all patches
    all_py(:,r) = r+[-ph, ph, ph, -ph];
    
    % Set x-coords specific for each patch
    pre_ice_px(1:2,r) = xlim(1);              % left-side of plot
    pre_ice_px(3:4,r) = evdates_regmean(r,1); % First ice

    icegrow_px(1:2,r) = evdates_regmean(r,1); % First ice
    icegrow_px(3:4,r) = evdates_regmean(r,3); % Max ice
    
    iceretr_px(1:2,r) = evdates_regmean(r,3); % Max ice
    iceretr_px(3:4,r) = evdates_regmean(r,5); % Break up

    icebrup_px(1:2,r) = evdates_regmean(r,5); % Break up
    icebrup_px(3:4,r) = evdates_regmean(r,6); % Ice free

    postice_px(1:2,r) = evdates_regmean(r,6); % Break up
    postice_px(3:4,r) = xlim(2);              % right-side of plot

end

patch(pre_ice_px, all_py, noicecol, 'parent', ax2, ...
                                    'edgecolor', 'none');
patch(icegrow_px, all_py, growcol, 'parent', ax2, ...
                                   'edgecolor', 'none');
patch(iceretr_px, all_py, retrcol, 'parent', ax2, ...
                                   'edgecolor', 'none');
patch(icebrup_px, all_py, brupcol, 'parent', ax2, ...
                                   'edgecolor', 'none');
patch(postice_px, all_py, noicecol, 'parent', ax2, ...
                                    'edgecolor', 'none');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
% Create patches for variability of each key event
all_py = zeros(4, Nregs);
firstice_px = zeros(4, Nregs);
max_ice_px = zeros(4, Nregs);
breakup_px = zeros(4, Nregs);
lastice_px = zeros(4, Nregs);
postice_px = zeros(4, Nregs);
for r=1:Nregs
    
    % Set y-coords for all patches
    all_py(:,r) = r+[-ph, ph, ph, -ph];
    
    % Set x-coords specific for each patch
    firstice_px(1:2,r) = evdates_reglow(r,1);
    firstice_px(3:4,r) = evdates_regupp(r,1); 

    max_ice_px(1:2,r) = evdates_reglow(r,3);
    max_ice_px(3:4,r) = evdates_regupp(r,3);
    
    breakup_px(1:2,r) = evdates_reglow(r,5);
    breakup_px(3:4,r) = evdates_regupp(r,5);

    lastice_px(1:2,r) = evdates_reglow(r,6);
    lastice_px(3:4,r) = evdates_regupp(r,6);

end

patch(firstice_px, all_py, evvarcol, 'parent', ax2, ...
                                     'edgecolor', 'none', ...
                                     'facealpha', 0.5);
patch(max_ice_px, all_py, evvarcol, 'parent', ax2, ...
                                    'edgecolor', 'none', ...
                                    'facealpha', 0.5);
patch(breakup_px, all_py, evvarcol, 'parent', ax2, ...
                                    'edgecolor', 'none', ...
                                    'facealpha', 0.5);
patch(lastice_px, all_py, evvarcol, 'parent', ax2, ...
                                    'edgecolor', 'none', ...
                                    'facealpha', 0.5);


% Plot event dates over the top -------------------

% Initialize arrays to hold x and y values
% - we can can plot all regions using one line for each region
%   if we concatenate the x and y values with NaNs between each region
all_y3 = zeros(Nregs*3,1);
firstice_x = zeros(Nregs*3,1);
maxice_x = zeros(Nregs*3,1);
breakup_x = zeros(Nregs*3,1);
icefree_x = zeros(Nregs*3,1);
all_y6 = zeros(Nregs*6,1);
firstice_var_x = zeros(Nregs*6,1);
maxice_var_x = zeros(Nregs*6,1);
breakup_var_x = zeros(Nregs*6,1);
icefree_var_x = zeros(Nregs*6,1);


for r = 1:Nregs
    r0 = 3*(r-1);
    r00 = 6*(r-1);

    % Set y values for all events
    all_y3(r0+(1:3)) = r+[-ph, ph, NaN];
    all_y6(r00+(1:6)) = r+[-ph, ph, NaN,-ph, ph, NaN];
    
    % Set x values for each event indivivually
    firstice_x(r0+(1:2)) = evdates_regmean(r,1);
    firstice_x(r0+3) = NaN;

    maxice_x(r0+(1:2)) = evdates_regmean(r,3);
    maxice_x(r0+3) = NaN;

    breakup_x(r0+(1:2)) = evdates_regmean(r,5);
    breakup_x(r0+3) = NaN;

    icefree_x(r0+(1:2)) = evdates_regmean(r,6);
    icefree_x(r0+3) = NaN;

    firstice_var_x(r00+(1:2)) = evdates_reglow(r,1);
    firstice_var_x(r00+3) = NaN;
    firstice_var_x(r00+(4:5)) = evdates_regupp(r,1);
    firstice_var_x(r00+6) = NaN;

    maxice_var_x(r00+(1:2)) = evdates_reglow(r,3);
    maxice_var_x(r00+3) = NaN;
    maxice_var_x(r00+(4:5)) = evdates_regupp(r,3);
    maxice_var_x(r00+6) = NaN;

    breakup_var_x(r00+(1:2)) = evdates_reglow(r,5);
    breakup_var_x(r00+3) = NaN;
    breakup_var_x(r00+(4:5)) = evdates_regupp(r,5);
    breakup_var_x(r00+6) = NaN;

    icefree_var_x(r00+(1:2)) = evdates_reglow(r,6);
    icefree_var_x(r00+3) = NaN;
    icefree_var_x(r00+(4:5)) = evdates_regupp(r,6);
    icefree_var_x(r00+6) = NaN;

end

% Plot lines for mean dates of each event
line(firstice_x, all_y3, 'parent', ax2, 'color', 'black', ...
                                       'linewidth', 2, ...
                                       'linestyle', '--');
                                 
line(maxice_x, all_y3, 'parent', ax2, 'color', 'black', ...
                                     'linewidth', 2, ...
                                     'linestyle', '--');
                                 
line(breakup_x, all_y3, 'parent', ax2, 'color', 'black', ...
                                      'linewidth', 2, ...
                                      'linestyle', '--');

line(icefree_x, all_y3, 'parent', ax2, 'color', 'black', ...
                                      'linewidth', 2, ...
                                      'linestyle', '--');

% Plot lines for mean +/- standard deviations for each event
line(firstice_var_x, all_y6, 'parent', ax2, 'color', 'black', ...
                                           'linewidth', 1, ...
                                           'linestyle', ':');
                                 
line(maxice_var_x, all_y6, 'parent', ax2, 'color', 'black', ...
                                         'linewidth', 1, ...
                                         'linestyle', ':');
                                 
line(breakup_var_x, all_y6, 'parent', ax2, 'color', 'black', ...
                                          'linewidth', 1, ...
                                          'linestyle', ':');

line(icefree_var_x, all_y6, 'parent', ax2, 'color', 'black', ...
                                          'linewidth', 1, ...
                                          'linestyle', ':');

% Add legend labels above top plot
legax = axes('parent', fig2, ...
             'pos', [lefmar,1-(topmar-0.025),pwd,0.5*topmar], ...
             'xlim', [0,1], 'ylim', [0,1], 'visible', 'off');


text_x = evdates_regmean(Nregs, [1,3,5,6]);
text_y = zeros(1,4)+Nregs+1.1*ph;
text_h = text(text_x, text_y, slieevent_labs([1,3,5,6]), ...
                              'parent', ax2, ...
                              'fontsize', 18, ...
                              'fontweight', 'bold', ...
                              'horizontalal', 'center', ...
                              'verticalal', 'bottom');



end