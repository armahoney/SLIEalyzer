function slie_ploteventdates_alongcoast(varargin)
% Function to plot the dates of key events within the landfast ice as
% functions of position along the coast with interannual variability
% indicated by the vertical extent of colored polygons around the mean
%
% USAGE
%      slie_ploteventdates_alongcoast
%      slie_ploteventdates_alongcoast('dataset', dataset)
%      slie_ploteventdates_alongcoast('asadir', asadir)
%      slie_ploteventdates_alongcoast('cfgfile', cfgfile)
%      slie_ploteventdates_alongcoast('evfolder', evfolder)
%      slie_ploteventdates_alongcoast('frzthaw', frzthaw)
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
%   ev2plot : array specifying the events to be included in plot
%             default: [1,2,5,6] (see below for event numbers)
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
% NOTES:
%    The events are numbered accordingly:
%           1: 'First ice' ('firstice')
%           2: 'Stable ice' ('stableice')
%           3: 'Max ice' ('maxice')
%           4: 'Unstable ice' ('unstableice')
%           5: 'Breakup' ('breakup')
%           6: 'Ice free' ('icefree')
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
            evfolder = varargin{a+1};
            a = a + 2;
        case 'ev2plot'
            ev2plot = varargin{a+1};
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

% Specify which events we're plotting
if ~exist('ev2plot')
    ev2plot = [1,2,5,6];
end
Nevp = numel(ev2plot);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Read the SLIE configuration file to get subregion and vector label info
cfg = slie_readconfigfile(sfolder, cfgfile);



% read the coast vector subregions file
[cvecregs, regnames, regbounds, Nregs] = ...
                   slie_readcvecregsfile(cfg.cvecregionsfile);

% Read the coast vector labels file
[cveclabs, labs_i, Nlabs] = slie_readcveclabsfile(cfg.cveclabelsfile);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Read in the different event data files
slieevent_strs = {'firstice', ...
                  'stableice', ...
                  'maxice', ...
                  'unstableice', ...
                  'breakup', ...
                  'icefree'};
Nevents = numel(slieevent_strs);

slieevent_labs = {'First ice', ...
                   'Stable ice', ...
                   'Max ice', ...
                   'Unstable ice', ...
                   'Breakup', ...
                   'Ice free'};

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
% Calculate means and stdevs for each vector using array maths
% Generate an array with NaNs swapped for zeros
fin_i = isfinite(evdates);
evdates_fin = evdates;
evdates_fin(~fin_i) = 0;

% Count the number of finite values across seasons
Nfin = sum(fin_i, 2);

% Sum the event dates across seaons
evdates_sum = sum(evdates_fin, 2);

% Calculate mean
evdates_mean = evdates_sum./Nfin;

% Calculate standard deviation
evdates_mean_rep = fin_i.*repmat(evdates_mean, [1,Nseasons, 1]);
evdates_stdev = sqrt( sum((evdates_fin - evdates_mean_rep).^2,2) ./ Nfin);
lessthantwo = Nfin < 2;
evdates_stdev(lessthantwo) = NaN;

% Remove unnecessary dimensions
evdates_mean = squeeze(evdates_mean);
evdates_stdev = squeeze(evdates_stdev);

% Calculate median and upper and lower quantiles
evdates_med = zeros(Nvecs, Nevents);
evdates_low = zeros(Nvecs, Nevents);
evdates_upp = zeros(Nvecs, Nevents);
q = 0.20;
qstr = sprintf('%02.0f%%', q*100);
for ev = 1:Nevents
    for v = 1:Nvecs
        evd = evdates_fin(v,:,ev);
        evd = evd(evd > 0);
        if numel(evd) > 2
            evsort = sort(evd);
            Nevd = numel(evd);
            midi = round(Nevd/2);
            evdates_med(v,ev) = evsort(midi);
            lowi = ceil(q*Nevd);
            evdates_low(v,ev) = evsort(lowi);
            uppi = floor((1-q)*Nevd);
            evdates_upp(v,ev) = evsort(uppi);
            notalot = 0;
        else
            evdates_low(v,ev) = NaN;
            evdates_upp(v,ev) = NaN;
        end
        notalot = 0;
    end
end

% Read in surface air temperature freeze/thaw data if frzthaw is specified
if frzthaw    
    % This requies that this datafile exists in this location!
    ft = load([sfolder filesep asadir filesep ...
               'FrzThaw_analysis/cp_freezethaw.mat'], ...
               'yr0', 'yr1', 'date0', 'date1', 'thaw_onset', 'frz_onset');
    ftsz = size(ft.thaw_onset);
    Nfrzthawyrs = ftsz(2);
           
    % Calculate mean and stdev thaw dates
    fini = isfinite(ft.thaw_onset);
    thaw_fin = dayofyear(ft.thaw_onset) + 365; % Add 365 to move into second calendar year of season
    thaw_fin(~fini) = 0;
    Nfin = sum(fini, 2);
    thaw_mean = sum(thaw_fin, 2)./Nfin;
    thaw_mean_rep = fini.*repmat(thaw_mean, [1,Nfrzthawyrs]);
    thaw_stdev = sqrt( sum((thaw_fin - thaw_mean_rep).^2, 2) ./ Nfin);
    lessthantwo = Nfin < 2;
    thaw_stdev(lessthantwo) = NaN;

    thaw_sort = sort(dayofyear(ft.thaw_onset), 2) + 365;
    midi = round(Nfin/2);
    lowi = ceil(q*Nfin);
    uppi = floor((1-q)*Nfin);
    midii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', midi);
    lowii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', lowi);
    uppii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', uppi);
    thaw_med = thaw_sort(midii);
    thaw_low = thaw_sort(lowii);
    thaw_upp = thaw_sort(uppii);

    fini = isfinite(ft.frz_onset);
    frz_fin = dayofyear(ft.frz_onset);
    frz_fin(~fini) = 0;
    Nfin = sum(fini, 2);
    frz_mean = sum(frz_fin, 2)./Nfin;
    frz_mean_rep = fini.*repmat(frz_mean, [1,Nfrzthawyrs]);
    frz_stdev = sqrt( sum((frz_fin - frz_mean_rep).^2, 2) ./ Nfin);
    lessthantwo = Nfin < 2;
    frz_stdev(lessthantwo) = NaN;
    
    frz_sort = sort(dayofyear(ft.frz_onset), 2);
    midi = round(Nfin/2);
    lowi = ceil(q*Nfin);
    uppi = floor((1-q)*Nfin);
    midii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', midi);
    lowii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', lowi);
    uppii = sub2ind([Nvecs, Nseasons], (1:Nvecs)', uppi);
    frz_med = frz_sort(midii);
    frz_low = frz_sort(lowii);
    frz_upp = frz_sort(uppii);
    
    disp('Loaded freeze / thaw analysis data');
end

% Report failure rates for events 
disp('Events failure rates:');
for ev = 1:Nevents
    Nfail = sum(~isfinite(evdates(:,:,ev)), 'all');
    fr = 100*Nfail / (Nvecs*Nseasons);
    fprintf('%s: %4.1f\n', slieevent_labs{ev}, fr);
    notalot = 0;
end


% ========================= FIGURE PLOTTING ===============================                          
% Plot mean event dates with interannual variability for all coast vectors
% - X-axis represents location along coast, divided up by regions 
%   (this is achieved by creating separate axes for each region)
% - Y-axis represents day of year

fig_wd = 17.5;
fig_ht = 9.5;
fig = figure('units', 'inches', ...
              'position', [1,5,fig_wd,fig_ht], ...
              'paperposition', [0,0,fig_wd,fig_ht]);
lefmar = 0.06;
rigmar = 0.12;
topmar = 0.1;
botmar = 0.125;
hgap = 0.0025;
hspace = 1-(lefmar+rigmar+((Nregs-1)*hgap));
vgap = 0;
vspace = 1-(topmar+botmar);

% Big blank axes to hold labels outside of other axes
bbax = axes('parent', fig, 'pos', [0,0,1,1], 'color', 'white', ...
                                              'xtick', [], 'ytick', [], ...
                                              'xcolor', 'white', ...
                                              'ycolor', 'white');

% Set up y-axis label stuff
x_yrs = [ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1];
x_mos = [ 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8];
x_das = zeros(1,12)+1;
xtick = datenum(x_yrs, x_mos, x_das);
xticklab = datestr(xtick, 'mmm dd');
xlim = [datenum(0,9,1),datenum(1,8,31)];

ax = zeros(Nregs,1);

% Set colors for each event
slieevent_cols = {'magenta', ...
                  'blue', ...
                  'black', ...
                  'cyan' ...
                  'green', ...
                  'red'};

% Set colors for freeze / thaw dates
thaw_col = [1.0, 0.5, 0.0]; % Orange
thaw_style = '--';
frz_col = [0.5, 0.0, 1.0]; % Purply colour;
frz_style = '--';
              
% Initialize location of first set of axes for first region
ax_x0 = lefmar-hgap;
pwd = 0;

for r=1:Nregs
    % Figure out which coast vectors are in this region
    ri = cvecregs == r;
    rv = find(ri)';
    rx0 = min(rv);
    rx1 = max(rv);    
    rNv = sum(ri);
    rx = shiftdim(1:rNv,1);
        
    % Set up axes for this region
    % - allow a little gap between axes
    ax_x0 = ax_x0 + pwd + hgap;
    pwd = rNv*hspace/Nvecs;
    ax_y0 = botmar;
    pht = vspace;
    
    % Only use yaxis tick labels for leftmost axes
    % (set other ytick labels to empty array)
    if r > 1
        xticklab = {};
    end
    
    % Figure out x-axis ticks and labels
    rlabs_i = (labs_i >= rx0) & (labs_i <= rx1);
    rlabs = cveclabs(rlabs_i);
    rxtick = 1+ labs_i(rlabs_i) - rx0;

    % Create set of axes for this region
    ax(r) = axes('parent', fig, 'pos', [ax_x0, ax_y0, pwd, pht], ...
                 'xtick', rxtick, 'xticklabel', {}, ...
                 'xlim', [1,rNv], 'xgrid', 'on', ...
                 'ylim', xlim, 'ytick', xtick, 'yticklabel', xticklab, ...
                 'ygrid', 'on', ...
                 'fontsize', 14, 'fontweight', 'bold', ...
                 'box', 'on', 'linewidth', 1);

    % Add x-tick labels as text objects
    % (this allows more flexibility than trying to assign labels
    %  via the axes definition)
    tx_norm = ax_x0 + rxtick*pwd/rNv;
    ty_norm = zeros(numel(rlabs),1) + ax_y0 - 0.01;
    text(tx_norm, ty_norm, rlabs, 'parent', bbax, ...
                                  'horizontalal', 'right', ...
                                  'verticalal', 'middle', ...
                                  'rotation', 45, ...
                                  'units', 'norm', ...
                                  'fontsize', 12, 'fontweight', 'bold');

    % Insert newlines to split labels onto multiple lines
    if any(regexp(regnames{r}, '-'))
        rexp0 = '([a-zA-Z\.])\s-\s([a-zA-Z\.])';
        rexp1 = '$1\nto\n$2';
    else
        rexp0 = '([a-zA-Z])\s([a-zA-Z])';
        rexp1 = '$1\n$2';
    end
    reglab = regexprep(regnames{r}, rexp0, rexp1);
    
    % Use region label as plot title
    title(ax(r), reglab, 'fontsize', 14, 'fontweight', 'bold');
                               
    if frzthaw
        % Draw patch and lines for thaw onset
        %uppbound = thaw_mean(rv) + thaw_stdev(rv);
        %lowbound = thaw_mean(rv) - thaw_stdev(rv);
        uppbound = thaw_upp(rv);
        lowbound = thaw_low(rv);
        
        px = [rx; rx(end:-1:1)];
        py = [uppbound; lowbound(end:-1:1)];
        patch(px,py, thaw_col, 'parent', ax(r), ...
                    'edgecolor', 'none', ...
                    'facealpha', 0.3);
                
        line(rx, uppbound, 'parent', ax(r), ...
                            'color', thaw_col, ...
                            'linestyle', thaw_style, ...
                            'linewidth' ,1);
        line(rx, lowbound, 'parent', ax(r), ...
                             'color', thaw_col, ...
                             'linestyle', thaw_style, ...
                             'linewidth' ,1);
        
        line(rx, thaw_med(rv), 'parent', ax(r), ...
                                'displayname', 'Thaw onset', ...
                                'color', thaw_col, ...
                                'linestyle', thaw_style, ...
                                'linewidth', 2);
    

        % Draw patch and lines for freeze onset
        %uppbound = frz_mean(rv) + frz_stdev(rv);
        %lowbound = frz_mean(rv) - frz_stdev(rv);
        uppbound = frz_upp(rv);
        lowbound = frz_low(rv);
        
        px = [rx; rx(end:-1:1)];
        py = [uppbound; lowbound(end:-1:1)];
        patch(px,py, frz_col, 'parent', ax(r), ...
                    'edgecolor', 'none', ...
                    'facealpha', 0.3);
                
        line(rx, uppbound, 'parent', ax(r), ...
                            'color', frz_col, ...
                            'linestyle', frz_style, ...
                            'linewidth' ,1);
        line(rx, lowbound, 'parent', ax(r), ...
                             'color', frz_col, ...
                             'linestyle', frz_style, ...
                             'linewidth' ,1);
        
        line(rx, frz_med(rv), 'parent', ax(r), ...
                               'displayname', 'Freeze onset', ...
                               'color', frz_col, ...
                               'linestyle', frz_style, ...
                               'linewidth', 2);

    
    end
    
    % Now plot the events we're told to by ev2plot variable
    for ee = 1:Nevp
        % Assign short variable name to idenitfy the event
        % being plotted in this loop iteration
        e = ev2plot(ee);
                
        % Smooth dates along the coast
        smwid = 101;
       
        % Calculate mean event date 
        % and standard deviation 
        % across all years
        ev_mean = evdates_mean(ri,e);
        ev_stdev = evdates_stdev(ri,e);
        ev_med = evdates_med(ri,e);
        ev_low = evdates_low(ri,e);
        ev_upp = evdates_upp(ri,e);

        
        % Apply a running mean along the coast
        % - using a technique that handles and 
        %   preserves any NaNs in the data
        ev_mean_sm = rmean_nan(ev_mean, smwid);
        ev_stdev_sm = rmean_nan(ev_stdev, smwid);
        ev_med_sm = rmean_nan(ev_med, smwid);
        ev_low_sm = rmean_nan(ev_low, smwid);
        ev_upp_sm = rmean_nan(ev_upp, smwid);
        
        % Calculate the mean plus and minus the standard deviation
        %uppbound = ev_mean_sm + ev_stdev_sm;
        %lowbound = ev_mean_sm - ev_stdev_sm;
        uppbound = ev_upp_sm;
        lowbound = ev_low_sm;
        

        % Find NaNs and finite data points
        dev_nan = find(isnan(uppbound));     
        dev_fin = find(isfinite(uppbound));

        % Draw patches for region between +/- 1 standard dev
        % - identify continuous segments of data separated by NaNs 
        if any(dev_fin)
            seg_start = dev_fin(1);
            seg_end = -1;
            while seg_end < dev_fin(end)
                next_nan = find(dev_nan > seg_start, 1, 'first');
                if any(next_nan)
                    seg_end = dev_nan(next_nan)-1;
                else
                    seg_end = dev_fin(end);
                end
            
                segi = seg_start:seg_end;
                segrev = seg_end:-1:seg_start;
    
               
                px = [rx(segi); rx(segrev)];
                py = [uppbound(segi); lowbound(segrev)];
                patch(px,py, slieevent_cols{e}, 'parent', ax(r), ...
                    'edgecolor', 'none', ...
                    'facealpha', 0.3);
                
                next_fin = find(dev_fin > seg_end, 1, 'first');
                seg_start = dev_fin(next_fin);
            end
        end
        
        % Add lines to bound the upper and lower edges of the patches
        line(rx, uppbound, 'parent', ax(r), ...
                            'color', slieevent_cols{e}, ...
                            'linestyle', '-', ...
                            'linewidth' ,1);
        line(rx, lowbound, 'parent', ax(r), ...
                            'color', slieevent_cols{e}, ...
                            'linestyle', '-', ...
                            'linewidth' ,1);
        
        % Plot mean event date over the top
        line(rx, ev_med_sm, 'parent', ax(r), ...
            'displayname', slieevent_strs{e}, ...
            'color', slieevent_cols{e}, ...
            'linestyle', '-', ...
            'linewidth', 2);

        
    end       


end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Now add a legend ------------

% Initialize legend axes
leg_x0 = 1-(0.9*rigmar);
leg_y0 = botmar;
leg_wd = 0.8*rigmar;
leg_ht = vspace;
legax = axes('parent', fig, 'pos', [leg_x0,leg_y0,leg_wd,leg_ht], ...
                             'xlim', [0,1], 'ylim', [0,1], ...
                             'visible', 'off');

% Assign sizing and spacing parameters for legend entries
lgap = 0.1;
Nleg = Nevp + (2*frzthaw);
lht = (1 - ((Nleg-1)*lgap))/Nleg;

lx0 = 0;
lx1 = 0.3;
llx = lx1+0.1;

% Create legend entries for each event
for ee=1:Nevp
    e = ev2plot(ee);
    
    ly0 = (ee+frzthaw-1)*(lht+lgap);
    ly1 = ly0+lht;
    
    lpx = [lx0,lx1,lx1,lx0];
    lpy = [ly0,ly0,ly1,ly1];
    
    % Draw legend entry patch
    patch(lpx,lpy, slieevent_cols{e}, 'parent', legax, ...
                                      'edgecolor', 'none', ...
                                      'facealpha', 0.3);
    
    % +/- 1 s.d lines
    line([lx0,lx1],[ly1,ly1], 'parent', legax, ...
                              'color', slieevent_cols{e}, ...
                              'linestyle', '-', ...
                              'linewidth' ,1);
    line([lx0,lx1],[ly0,ly0], 'parent', legax, ...
                              'color', slieevent_cols{e}, ...
                              'linestyle', '-', ...
                              'linewidth' ,1);
    
    % Mean line
    lym = (ly0+ly1)/2;
    line([lx0,lx1], [lym,lym], 'parent', legax, ...
                               'displayname', slieevent_strs{e}, ...
                               'color', slieevent_cols{e}, ...
                               'linestyle', '-', ...
                               'linewidth', 2);
    
    
    % Labels
    text(lx1,ly1, ['upp ' qstr], 'parent', legax, ...
                             'fontsize', 10, 'fontweight', 'bold', ...
                             'horizontalal', 'left', ...
                             'verticalal', 'middle');
    text(lx1,lym, slieevent_labs{e}, 'parent', legax, ...
                                     'fontsize', 12, ...
                                     'fontweight', 'bold', ...
                                     'horizontalal', 'left', ...
                                     'verticalal', 'middle');
    text(lx1,ly0, ['low ' qstr], 'parent', legax, ...
                             'fontsize', 10, 'fontweight', 'bold', ...
                             'horizontalal', 'left', ...
                             'verticalal', 'middle');
    
end

% Add legend entries for freeze/thaw dates if they are plotted
% - these are included above and below the event legend entries
if frzthaw
    % Add legend item for Thaw Onset
    ly0 = (Nevp+1)*(lht+lgap);
    ly1 = ly0+lht;
    
    lpx = [lx0,lx1,lx1,lx0];
    lpy = [ly0,ly0,ly1,ly1];
    
    % Draw legend entry patch
    patch(lpx,lpy, thaw_col, 'parent', legax, ...
                             'edgecolor', 'none', ...
                             'facealpha', 0.3);
    
    % +/- 1 s.d lines
    line([lx0,lx1],[ly1,ly1], 'parent', legax, ...
                              'color', thaw_col, ...
                              'linestyle', thaw_style, ...
                              'linewidth' ,1);
    line([lx0,lx1],[ly0,ly0], 'parent', legax, ...
                              'color', thaw_col, ...
                              'linestyle', thaw_style, ...
                              'linewidth' ,1);
    
    % Mean line
    lym = (ly0+ly1)/2;
    line([lx0,lx1], [lym,lym], 'parent', legax, ...
                               'displayname', 'Thaw onset', ...
                               'color', thaw_col, ...
                               'linestyle', thaw_style, ...
                               'linewidth', 2);
    
    
    % Labels
    text(lx1,ly1, ['upp ' qstr], 'parent', legax, ...
                             'fontsize', 10, 'fontweight', 'bold', ...
                             'horizontalal', 'left', ...
                             'verticalal', 'middle');
    text(lx1,lym, 'Thaw onset', 'parent', legax, ...
                                'fontsize', 12, ...
                                'fontweight', 'bold', ...
                                'horizontalal', 'left', ...
                                'verticalal', 'middle');
    text(lx1,ly0, ['low ' qstr], 'parent', legax, ...
                             'fontsize', 10, 'fontweight', 'bold', ...
                             'horizontalal', 'left', ...
                             'verticalal', 'middle');
        

                         
                         
    % Add legend item for Freeze Onset
    ly0 = 0;
    ly1 = ly0+lht;
    
    lpx = [lx0,lx1,lx1,lx0];
    lpy = [ly0,ly0,ly1,ly1];
    
    % Draw legend entry patch
    patch(lpx,lpy, frz_col, 'parent', legax, ...
                            'edgecolor', 'none', ...
                            'facealpha', 0.3);
    
    % +/- 1 s.d lines
    line([lx0,lx1],[ly1,ly1], 'parent', legax, ...
                              'color', frz_col, ...
                              'linestyle', frz_style, ...
                              'linewidth' ,1);
    line([lx0,lx1],[ly0,ly0], 'parent', legax, ...
                              'color', frz_col, ...
                              'linestyle', frz_style, ...
                              'linewidth' ,1);
    
    % Mean line
    lym = (ly0+ly1)/2;
    line([lx0,lx1], [lym,lym], 'parent', legax, ...
                               'displayname', 'Freeze onset', ...
                               'color', frz_col, ...
                               'linestyle', frz_style, ...
                               'linewidth', 2);
    
    
    % Labels
    text(lx1,ly1, ['upp ' qstr], 'parent', legax, ...
                                 'fontsize', 10, 'fontweight', 'bold', ...
                                 'horizontalal', 'left', ...
                                 'verticalal', 'middle');
    text(lx1,lym, 'Freeze onset', 'parent', legax, ...
                                'fontsize', 12, ...
                                'fontweight', 'bold', ...
                                'horizontalal', 'left', ...
                                'verticalal', 'middle');
    text(lx1,ly0, ['low ' qstr], 'parent', legax, ...
                                 'fontsize', 10, 'fontweight', 'bold', ...
                                 'horizontalal', 'left', ...
                                 'verticalal', 'middle');
end


end

