function [fig, fig2] = slie_ploteventdates_byyear(varargin)
% A function to plot timeseries of sub-regional mean dates of key SLIE events


% USAGE:
%    [fig, fig2] = slie_ploteventdates_byyear()
%    [fig, fig2] = slie_ploteventdates_byyear('dataset', dataset)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'cfgfile', cfgfile)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'asadir', asadir)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'ev2plot', ev2plot)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'shadebounds', shadebounds)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'norm', norm)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'cveci', cveci)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'regioncol', rcol)
%    [fig, fig2] = slie_ploteventdates_byyear(__, 'hatching', hatching)
%
% INPUTS:
%     dataset : Character string specifying dataset to plot
%               See slie_get_basefolder.m for accepted values.
%               Default: "beau"
%     cfgfile : string specifying the name of the configuraion file stored
%               in the region folder containing essential information for
%               performing all analysis of SLIE data. 
%               See slie_readconfigfile.m for more details
%               default: 'slieconfig.txt'
%      asadir : string specifying folder name in region folder where
%               the data files from the all-seasons analysis are saved
%               default: 'AllSeasonsAnalysis2023'
%      evplot : integer array specifying which events to plot
%               Default: [1,2,5,6]
%               (i.e., FreezeUp, Stable Ice, Break Up, Last Ice)
% shadebounds : character string specifying what measures to use for upper
%               and lower bounds of shaded regions. Can be either:
%               - 'stdev' for +/- 1 standard deviation
%               - 'qtile' for upper and lower 10% quantiles
%                Default: 'qtile'
%        norm : flag specifying whether or not to subtract mean dates
%               Default: 0 (Don't normalize)
%        rcol : Color specification for single-event plots
%               Default: 'blue'
%    hatching : Angle to use for hatch-fill of polygon patches
%               Default: 'none'


% Deal with variable arguments
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
        case 'ev2plot'
            ev2plot = varargin{a+1};
            a = a + 2;
        case 'shadebounds'
            shadebounds = varargin{a+1};
            a = a + 2;
        case 'norm'
            norm = 1;
            a = a + 1;
        case 'cveci'
            cveci = varargin{a+1};
            a = a + 2;
        case 'regioncol'
            rcol = varargin{a+1};
            a = a + 2;
        case 'hatching'
            hatching = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end    
    
        
% ----------------------------------------------------------------
% Assign default values if not specified
if ~exist('dataset', 'var')
    dataset = 'beau';
end

if ~exist('cfgfile', 'var'), cfgfile = 'sliealyzer_config.txt'; end

if ~exist('asadir', 'var'), asadir = 'AllSeasonsAnalysis2023'; end

if ~exist('ev2plot', 'var'), ev2plot = [1,2,5,6]; end

if ~exist('shadebounds', 'var'), shadebounds = 'qtile'; end

if ~exist('norm', 'var'), norm = 0; end

if ~exist('rcol', 'var'), rcol = 'blue'; end

if ~exist('hatching', 'var'), hatching = 'none'; end


% ----------------------------------------------------------------
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);


% Read the SLIE configuration file to get various useful things
cfg = slie_readconfigfile(sfolder, cfgfile);



% Read in the regional mean event data files
evfolder = [sfolder filesep asadir filesep 'SLIE_Events'];


% Assign event names and labels
slieevent_strs_all = {'01_firstice', ...
                      '02_stableice', ...
                      '03_maxice', ...
                      '04_unstableice', ...
                      '05_breakup', ...
                      '06_icefree', ...
                      '07_firstwid', ...
                      '08_lastwid'};
slieevent_labs_all = {'First ice', ...
                      'Stable ice', ...
                      'Maximum width', ...
                      'End of stable period', ...
                      'Breakup', ...
                      'Ice free', ...
                      'Wide ice', ...
                      'End of wide ice'};
evcol_all = {'magenta', ...
             'blue', ...
             'black', ...
             'cyan', ...
             'green', ...
             'red', ...
             [0.5,0.5,1.0], ...
             [1.0,0.5,0.5]};


slieevent_strs = slieevent_strs_all(ev2plot);
slieevent_labs = slieevent_labs_all(ev2plot);
evcol = evcol_all(ev2plot);
Nevents = numel(slieevent_strs);

              
% Initialize min, max, mean arrays
Nseasons = numel(cfg.seasons);
evmean = zeros(Nseasons, Nevents) + NaN;
evstdev = zeros(Nseasons, Nevents) + NaN;
evmin = zeros(Nseasons, Nevents) + NaN;
evmax = zeros(Nseasons, Nevents) + NaN;
evmedian = zeros(Nseasons, Nevents) + NaN;
evlow = zeros(Nseasons, Nevents) + NaN;
evupp = zeros(Nseasons, Nevents) + NaN;

% Set quantile bounds
q = 0.2;
qlow = q;
qupp = 1-q;
qlowstr = sprintf('Lower %02d%%', qlow*100);
quppstr = sprintf('Upper %02d%%', qlow*100);

for e = 1:Nevents
    evfile = [evfolder, filesep, ...
              'event_',slieevent_strs{e},'_datenum.csv'];

    [evdates, Nvecs, ~, ~] = slie_readeventdates(evfile);
    
    if exist('cveci', 'var')
        evdates_cv = evdates(cveci,:);
    else
        evdates_cv = evdates;
    end

    % Normalize by subtracting mean if spec'd in function call
    if norm
        fini = isfinite(evdates_cv);
        evdates_fin = evdates_cv;
        evdates_fin(~fini) = 0;
        normmean = sum(evdates_fin, 2)./sum(1*fini,2);
        evdates_cv = evdates_cv - repmat(normmean,1,Nseasons);
    end
    
    % Calculate seasonal mean for each event
    for s=1:Nseasons
        fini = isfinite(evdates_cv(:,s));
        if any(fini)
            evmean(s,e) = mean(evdates_cv(fini,s));
            evstdev(s,e) = std(evdates_cv(fini,s));
            evmin(s,e) = min(evdates_cv(fini,s));
            evmax(s,e) = max(evdates_cv(fini,s));
            evmedian(s,e) = median(evdates_cv(fini,s));
            evsort = sort(evdates_cv(fini,s));
            lowi = max(1, floor(qlow*sum(fini)));
            evlow(s,e) = evsort(lowi);
            uppi = ceil(qupp*sum(fini));
            evupp(s,e) = evsort(uppi);
        end
    end

    notalot= 0;
end
    

% Get ready to plot things ---------------

yr0 = zeros(Nseasons,1);
for s=1:Nseasons
    yr0(s) = str2double(cfg.seasons{s}(1:4));
end


% fig = zeros(Nevents,1);
% ax = zeros(Nevents,1);
% 
% for e = 1:Nevents
%     fig(e) = figure;
%     
%     ax(e) = axes('parent', fig(e));
%     
%     title(slieevent_labs{e});
%     for r=1:Nregs
%         line(yr0, evdates(r,1:Nseasons,e), 'parent', ax(e), ...
%                                   'displayname', subregs{r});
%     end
% end


% =========================================================================
% Plot each event on separate axes
% % switch lower(region)
% %     case 'chuk'
% %         rcol = 'red';
% %         hatchangle = -45;
% %     case 'beau'
% %         rcol = 'blue';
% %         hatchangle = 45;
% %     case 'chuknic'
% %         rcol = 'red';
% %         hatchangle = -45;
% %     case 'beaunic'
% %         rcol = 'blue';
% %         hatchangle = 45;
% %     case 'chukasip'
% %         rcol = 'red';
% %         hatchangle = -45;
% %     case 'beauasip'
% %         rcol = 'blue';
% %         hatchangle = 45;
% % end


fig = figure;
lefmar = 0.1;
rigmar = 0.1;
topmar = 0.3;
botmar = 0.3;
hspace = 1-(lefmar+rigmar);
vspace = 1-(topmar+botmar);
hgap = 0.05;
vgap = 0.04;
Ncols = Nevents;
Nrows = 1;
pwd = (hspace - (Ncols-1)*hgap)/Ncols;
pht = (vspace - (Nrows-1)*vgap)/Nrows;

bbax = axes('parent', fig, 'pos', [0,0,1,1], 'visible', 'off');
axg = zeros(Nevents,1);
ax = zeros(Nevents,1);

for e = 1:Nevents

    % Set position of axes
    c = ceil(e/Nrows);
    r = mod(e,Nrows);
    r = r + Nrows*(r == 0);
    %disp(['r=',num2str(r),', c=',num2str(c)]);
    x0 = lefmar+((c-1)*(pwd+hgap));
    y0 = 1-(topmar+(r*pht)+((r-1)*vgap));
    
    % Set x-axis labels, etc
    xt = yr0;
    %Nxt = numel(xt);
    %xtl = cell(Nxt,1);
    %for t=1:Nxt,
    %    xtl{t} = [{num2str(xt(t),'%4d'),},{['- ',num2str(xt(t)+1, '%4d')]}];
    %end
    xtl = cfg.seasons;
    xlim = [min(yr0)-0.5, max(yr0)+0.5];
    
    % Create blank grid axes (with gray gridlines)
    axg(e) = axes('parent', fig, 'pos', [x0,y0,pwd,pht], ...
                 'xtick', xt, 'xticklabel', [], 'xlim', xlim, ...
                 'xgrid', 'on', 'xcolor', [0.5,0.5,0.5], ...
                 'ygrid', 'on', 'ycolor', [0.5,0.5,0.5]);
    
    
    ax(e) = axes('parent', fig, 'pos', [x0,y0,pwd,pht], ...
                 'xtick', xt, 'xticklabel', [], 'xlim', xlim, ...
                 'fontsize', 14, 'fontweight', 'normal', ...
                 'color', 'none');
             
    title(slieevent_labs{e}, 'fontsize', 16, 'fontweight', 'bold');

    % Add text objects for xtick labels on bottom row       
    if r == Nrows
        for t= 1:numel(xt)
            xtN = x0 + pwd*(xt(t)-xlim(1))/(xlim(2)-xlim(1));
            text(xtN, y0-0.01, xtl{t}, 'parent', bbax, ...
                                       'units', 'normal', ...
                                       'fontsize', 14, ...
                                       'verticalal', 'middle', ...
                                       'horizontalal', 'right', ...
                                       'rotation', 90);
        end
    end
    
    
    % Plot standard deviation as a patch with top and bottom lines
    switch shadebounds
        case 'stdev'
            centerline = evmean(:,e);
            lowbound = evmean(:,e)-evstdev(:,e);
            uppbound = evmean(:,e)+evstdev(:,e);
        case 'qtile'
            centerline = evmedian(:,e);
            lowbound = evlow(:,e);
            uppbound = evupp(:,e);
        otherwise
            disp(['Unrecognized shadebounds value: ' shadebounds]);
            disp('Using upper/lower quantiles');
            centerline = evmedian(:,e);
            lowbound = evlow(:,e);
            uppbound = evupp(:,e);
    end

    fini = find(isfinite(lowbound) & isfinite(uppbound));
    revi = numel(fini):-1:1;
    px = [yr0(fini);yr0(fini(revi))];
    py = [lowbound(fini); uppbound(fini(revi))];
    shadepatch = patch(px, py, rcol, 'parent', ax(e), ...
                                     'edgecol', 'none', ...
                                     'facealpha', 0.25);

    if isnumeric(hatching)
        set(shadepatch, 'facecol', 'none');
        hatch(shadepatch, hatching, rcol, '-', 3,1);
    end
    
    lowlh = line(yr0, lowbound, 'parent', ax(e), ...
                                      'color', rcol, ...
                                      'linestyle', '-', ...
                                      'linewidth', 1);
    upplh = line(yr0, evupp(:,e), 'parent', ax(e), ...
                                     'color', rcol, ...
                                     'linestyle', '-', ...
                                     'linewidth', 1);
                                
    % Plot max a min lines as dashed lines
    minlh = line(yr0, evmin(:,e), 'parent', ax(e), ...
                                    'color', rcol, ...
                                    'linestyle', '--', ...
                                    'linewidth', 1);
    maxlh = line(yr0, evmax(:,e), 'parent', ax(e), ...
                                    'color', rcol, ...
                                    'linestyle', '--', ...
                                    'linewidth', 1);
                                     
    
    % Plot mean/median line in black
    centerlh = line(yr0, centerline, 'parent', ax(e), ...
                                     'color', rcol, ...
                                     'linestyle', '-', ...
                                     'linewidth', 2, ...
                                     'marker', '.', 'markersize', 8);
   
    ylim = get(ax(e), 'ylim');
    yt = ylim(1)+10:20:ylim(2);
    ytl = datestr(yt, 'mmm dd');
    set(ax(e), 'ytick', yt, 'yticklabel', ytl, 'ygrid', 'on');
    set(axg(e), 'ytick', yt, 'yticklabel', []);
    
end




% =========================================================================
% Plot all events on single axes
fig2 = figure;
lefmar = 0.10;
rigmar = 0.15;
topmar = 0.03;
botmar = 0.17;
hspace = 1-(lefmar+rigmar);
vspace = 1-(topmar+botmar);
hgap = 0.05;
vgap = 0.04;
Ncols = 1;
Nrows = 1;
pwd = (hspace - (Ncols-1)*hgap)/Ncols;
pht = (vspace - (Nrows-1)*vgap)/Nrows;


% Create blank grid axes (with gray gridlines)
xt = yr0;
%Nxt = numel(xt);
%xtl = cell(Nxt,1);
%for t=1:Nxt,
%    xtl{t} = [{num2str(xt(t),'%4d'),},{['- ',num2str(xt(t)+1, '%4d')]}];
%end
xtl = cfg.seasons;
xlim = [min(yr0)-0.5, max(yr0)+0.5];
ylim = [datenum(2,9,1),datenum(3,8,31)]-datenum(1,12,31);
ytv = [ 0, 9, 1,0,0,0; ...
        0,10, 1,0,0,0; ...
        0,11, 1,0,0,0; ...
        0,12, 1,0,0,0; ...
        1, 1, 1,0,0,0; ...
        1, 2, 1,0,0,0; ...
        1, 3, 1,0,0,0; ...
        1, 4, 1,0,0,0; ...
        1, 5, 1,0,0,0; ...
        1, 6, 1,0,0,0; ...
        1, 7, 1,0,0,0; ...
        1, 8, 1,0,0,0; ...
        1, 9, 1,0,0,0];
yt = datenum(ytv);   
ytl = datestr(yt, 'mmm dd');
axpos = [lefmar,botmar, hspace, vspace];
axg = axes('parent', fig2, 'pos', axpos, ...
           'xtick', xt, 'xticklabel', [], 'xlim', xlim, ...
           'xgrid', 'on', 'xcolor', [0.5,0.5,0.5], ...
           'ytick', yt, 'yticklabel', [], 'ylim', ylim, ...
           'ygrid', 'on', 'ycolor', [0.5,0.5,0.5]);
    
  
% Create separate axes for plotting data       
ax = axes('parent', fig2, 'pos', axpos, ...
          'xtick', xt, 'xticklabel', [], 'xlim', xlim, ...
          'ytick', yt, 'yticklabel', ytl, 'ylim', ylim, ...
          'fontsize', 14, 'fontweight', 'normal', ...
          'color', 'none');

% Initialize array to hold text for linear regression
regrstr = cell(Nevents,1);

for e = 1:Nevents
        
    % Plot standard deviation as a patch with top and bottom lines
    switch shadebounds
        case 'stdev'
            centerline = evmean(:,e);
            lowbound = evmean(:,e)-evstdev(:,e);
            uppbound = evmean(:,e)+evstdev(:,e);
        case 'qtile'
            centerline = evmean(:,e);
            lowbound = evlow(:,e);
            uppbound = evupp(:,e);
        otherwise
            disp(['Unrecognized shadebounds value: ' shadebounds]);
            disp('Using upper/lower quantiles');
            centerline = evmean(:,e);
            lowbound = evlow(:,e);
            uppbound = evupp(:,e);
    end


    fini = find(isfinite(lowbound) & isfinite(uppbound));
    revi = numel(fini):-1:1;
    px = [yr0(fini);yr0(fini(revi))];
    py = [lowbound(fini); uppbound(fini(revi))];
    shadepatch = patch(px, py, evcol{e}, 'parent', ax, ...
                                         'edgecol', 'none', ...
                                         'facealpha', 0.25);

    if isnumeric(hatching)
        set(shadepatch, 'facecolor', 'none');
        hatch(shadepatch, hatching, evcol{e}, '-', 3,1);
    end
    
    lowlh = line(yr0, lowbound, 'parent', ax, ...
                                      'color', evcol{e}, ...
                                      'linestyle', '-', ...
                                      'linewidth', 1);
    upplh = line(yr0, uppbound, 'parent', ax, ...
                                     'color', evcol{e}, ...
                                     'linestyle', '-', ...
                                     'linewidth', 1);

    
    % Plot mean line in black
    centerlh = line(yr0, centerline, 'parent', ax, ...
                                     'color', evcol{e}, ...
                                     'linestyle', '-', ...
                                     'linewidth', 2, ...
                                     'marker', '.', 'markersize', 8);

    % Calculate linear regression on center line
    
    rstats = regstats(centerline, yr0, 'linear');
    if rstats.fstat.pval < 0.05
        regrstr{e} = {sprintf('%4.1f d/yr', rstats.beta(2)), ...
                     sprintf('R^2 = %4.2f', rstats.rsquare)};
        line(yr0, rstats.yhat, 'parent', ax, ...
                               'color', evcol{e}, ...
                               'linestyle', '--', ...
                               'linewidth', 2, ...
                               'marker', 'none');
    else
        regrstr{e} = {'no significant', 'trend'};
    end
    fprintf('%s: %s, %s\n', slieevent_labs{e}, regrstr{e}{:});

    notalot = 0;
end

% Add text objects for xtick labels
for t= 1:numel(xt)
    %xtN = x0 + pwd*(xt(t)-xlim(1))/(xlim(2)-xlim(1));
    text(xt(t), ylim(1)-10, xtl{t}, 'parent', ax, ...
                                    'fontsize', 14, ...
                                    'verticalal', 'middle', ...
                                    'horizontalal', 'right', ...
                                    'rotation', 90);
end

% Now add a legend ------------

% Initialize legend axes
leg_x0 = 1-(0.9*rigmar);
leg_y0 = botmar;
leg_wd = 0.8*rigmar;
leg_ht = vspace;
legax = axes('parent', fig2, 'pos', [leg_x0,leg_y0,leg_wd,leg_ht], ...
                             'xlim', [0,1], 'ylim', [0,1], ...
                             'visible', 'off');

lgap = 0.1;
lht = (1 - ((Nevents-1)*lgap))/Nevents;

lx0 = 0;
lx1 = 0.3;
llx = lx1+0.1;
for e=1:Nevents
    
    ly0 = (e-1)*(lht+lgap);
    ly1 = ly0+lht;
    
    lpx = [lx0,lx1,lx1,lx0];
    lpy = [ly0,ly0,ly1,ly1];
    
    % Draw legend entry patch
    patch(lpx,lpy, evcol{e}, 'parent', legax, ...
                             'edgecolor', 'none', ...
                             'facealpha', 0.3);
    
    % upper and lower bounds lines
    line([lx0,lx1],[ly1,ly1], 'parent', legax, ...
                              'color', evcol{e}, ...
                              'linestyle', '-', ...
                              'linewidth' ,1);
    line([lx0,lx1],[ly0,ly0], 'parent', legax, ...
                              'color', evcol{e}, ...
                              'linestyle', '-', ...
                              'linewidth' ,1);
    
    % Center line
    lym = (ly0+ly1)/2;
    line([lx0,lx1], [lym,lym], 'parent', legax, ...
                               'displayname', slieevent_strs{e}, ...
                               'color', evcol{e}, ...
                               'linestyle', '-', ...
                               'linewidth', 2);
    
    
    % Labels
    switch shadebounds
        case 'stdev'
            upplab = '+1 s.d.';
            lowlab = '-1 s.d.';
        case 'qtile'
            upplab = quppstr;
            lowlab = qlowstr;
        otherwise
            upplab = quppstr;
            lowlab = qlowstr;
    end
    text(lx1,ly1, upplab, 'parent', legax, ...
                          'fontsize', 10, 'fontweight', 'bold', ...
                          'horizontalal', 'left', ...
                          'verticalal', 'middle');
    text(lx1,lym, slieevent_labs{e}, 'parent', legax, ...
                                     'fontsize', 12, ...
                                     'fontweight', 'bold', ...
                                     'horizontalal', 'left', ...
                                     'verticalal', 'middle');
    text(lx1,ly0, lowlab, 'parent', legax, ...
                          'fontsize', 10, 'fontweight', 'bold', ...
                          'horizontalal', 'left', ...
                          'verticalal', 'middle');
    
end


end






