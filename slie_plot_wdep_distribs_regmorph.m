function ax = slie_plot_wdep_distribs_regmorph(dataset, varargin)
% Function to plot the the water-depth-at-SLIE distributions according to:
% - sub-regions
% - months
% - coastal morphologies (headlands, bays, open coasts)
%
% USAGE
%     ax = slie_plot_wdep_distribs_regmorph(dataset)
%     ax = slie_plot_wdep_distribs_regmorph(dataset, asadir)
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
%             the data files from the all-seasons analysis will be saved
%             default: 'AllSeasonsAnalysis2023'
%
% OUTPUTS
%        ax : axes handles for output plots
%
% Andy Mahoney - October 2023

% ------------------------------------------------------------------------
% Check for variable arguments in function call
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
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end   

% Set defaults for anything not specified
if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

% Set default dataset if not specified
if ~exist('dataset', 'var')
    dataset = 'beau';
end

% -----------------------------------------------------------------------

% Identify the location where data are stored for specfied dataset
sfolder = slie_get_basefolder(dataset);


% Read the SLIE configuration file to get region and vector label info
cfg = slie_readconfigfile(sfolder, cfgfile);

% read the coast vector regions file
[cvecregs, regnames, regbounds, Nregs] = ...
                          slie_readcvecregsfile(cfg.cvecregionsfile);


% Set up headlands and bays info
[hbo_regs, hbonames] = slie_group_bays_headlands('region', dataset);


% Exclude lagoon regions for this analysis
% - these are identified by a decimal region value of 0.4;
% - (see slie_group_bays_headlands.m)
Nhbo = 3;
hbo_regs = floor(hbo_regs);
hbonames = hbonames(1:Nhbo); 


% Read in the water depth data
dfolder = [sfolder filesep asadir];
    
wdep = read_sliedatafile([dfolder filesep 'wdep_all.csv'], '%f');
%wdep = read_sliedatafile([dfolder filesep 'lfiw_all.csv'], '%f');


% Read in SLIE date data
sdates = read_sliedatafile([dfolder filesep 'sliedatenum_all.csv'], '%f');

% Figure out year, month and day for each data point
[syr, smo, sda] = datevec(sdates);

% =======================================================================
% Now get to business

% Bin by month and region and calculate distributions
mobins = [10,11,12,1,2,3,4,5,6,7; ...
          10,11,12,1,2,3,4,5,6,7];
Nmos = numel(mobins)/2;

dep0 = 0;
dep1 = 200;
binsz = 1;
depbins = [dep0:binsz:dep1, inf];
Nbins = numel(depbins);
regmo_dist = zeros(Nmos, Nregs, Nhbo, Nbins);


for r=1:Nregs
    regi = find(cvecregs == r);
    for hb = 1:Nhbo
        hboi = find(hbo_regs(regi) == hb);
        rmo = smo(hboi,:);
        rdep = wdep(hboi,:);
        for m=1:Nmos
            moi = find((rmo >= mobins(1,m)) & (rmo <= mobins(2,m)));
            distrib = histc(rdep(moi), depbins);
            regmo_dist(m,r,hb,:) = distrib;
            wtf = 0;
        end
    end
end
    
% Normalize distributions for each month and region
for r = 1:Nregs
    for m = 1:Nmos
        distsum = sum(sum(regmo_dist(m,r,:,:),4),3);
        regmo_dist(m,r,:,:) = regmo_dist(m,r,:,:)/distsum;
    end
end



% cols = [  0,  0,  0; ... % Oct - black
%           0,  0,  0; ... % Nov - black
%         255,  0,  0; ... % Dec - red
%         255,  0,128; ... % Jan - red/blue
%         128,  0,255; ... % Feb - blue/red
%           0,  0,255; ... % Mar - blue
%           0,128,255; ... % Apr - blue/green
%           0,255,255; ... % May - cyan
%           0,255,128; ... % Jun - green/blue
%           0,255,  0];    % Jul - green
% cols = cols/255;
      
fig = figure;
ax = zeros(Nmos,Nregs);
bbax = axes('parent', fig, 'pos', [0,0,1,1], ...
            'color', 'none', 'visible', 'off');


topmar = 0.1;
botmar = 0.05;
lefmar = 0.1;
rigmar = 0.05;
vspace = 1-(topmar+botmar);
hspace = 1-(lefmar+rigmar);
vgap = 0.02;
hgap = 0.01;
pht = (vspace-((Nmos-1)*vgap))/Nmos;
pwd = (hspace-((Nregs-1)*hgap))/Nregs;


dlines = zeros(Nmos, Nregs);

xlim = [0,50];
ylim = [0,0.3];


hbocols = [1  ,0.1,0.1; ... 
           0.1,0.9,0.1; ...
           0  ,0  ,  1];

for m=1:Nmos
    % Determine y postion of bottom edge of axes
    y0 = (1-topmar)-((m*pht)+(m-1)*(vgap));
    
    % Add month labels in left margin
    if mobins(1,m) == mobins(2,m)    
        monstr = datestr([1999,mobins(1,m),1,0,0,0], 'mmm');
    else
        monstr = {datestr([1999,mobins(1,m),1,0,0,0], 'mmm -'), ...
                  datestr([1999,mobins(2,m),1,0,0,0], 'mmm')};
    end
    text(lefmar/2, y0+pht/2, monstr, 'parent', bbax, ...
                                     'fontsize', 14, ...
                                     'fontweight', 'bold');
                  
    for r=1:Nregs
        % Determine x position of left edge of axes
        x0 = lefmar+(r-1)*(pwd+hgap);
        
        xtick = xlim(1):10:xlim(2);
        ytick = ylim(1):0.1:ylim(2);
        
        % Sort out x-tick labels for bottom row
        if (m == Nmos)
            xtickstr = num2str(shiftdim(xtick,1), '%3d');
        else
            xtickstr = '';
        end
        
        % Sort out y-tick labels for leftmost column
        if (r == 1)
            ytickstr = num2str(shiftdim(ytick,1), '%3.1f');
        else
            ytickstr = '';
        end
        
        % Draw axes
        ax(m,r) = axes('parent', fig, 'pos', [x0,y0,pwd,pht], ...
                                      'xlim', xlim, 'xtick', xtick, ...
                                      'xticklabel', xtickstr, 'xgrid', 'on', ...
                                      'ylim', ylim, 'ytick', ytick, ...
                                      'yticklabel', ytickstr, 'ygrid', 'on');
        
                        
        % Draw distribution                          
        d = squeeze(regmo_dist(m,r,:,1:end-1))';
%         lh1 = line(depbins, d, 'parent', ax(m,r), ...
%                          'displayname', monstr, ...
%                          'color', 'blue', ...
%                          'linewidth', 1);
        db = depbins(1:end-1) + (binsz/2);
        bh = bar(db, d, 1, 'stacked', 'parent', ax(m,r));
        
        % Set bar colors and displaynames
        for h = 1:Nhbo
            set(bh(h), 'displayname', hbonames{h}, ...
                       'facecolor', hbocols(h,:));
        end
                 
        % Add a horizonal line as a y-scale guide
        roundval = 0.05;
        dstack = sum(d, 2);
        ymax = roundval*ceil(max(dstack(2:end)/roundval));
        if isnan(ymax), ymax = 1; end
        if ymax == 0, ymax = roundval; end
        
        hlh = line(xlim, [ymax,ymax], 'parent', ax(m,r), ...
                                          'color', 'black', ...
                                          'linewidth', 1, ...
                                          'linestyle', '--');
                                                                            
        th1 = text(xlim(2),ymax,[num2str(100*ymax,'%2.0f'),'%'], ...
                                     'horizontalal', 'right', ...
                                     'verticalal', 'top');

        % Draw the segment from 0-1m in grey
%         lh3 = line(depbins(1:2), d(1:2), 'parent', ax(m,r), ...
%                                         'displayname', monstr, ...
%                                         'color', 'red', ...
%                                         'linewidth', 1);
        px = [db(1:2),db(2:-1:1)]-(binsz/2);
        py = [0,0,dstack(1),dstack(1)];
        ph = patch(px, py, [0.5,0.5,0.5], 'parent', ax(m,r), ...
                                          'displayname', monstr, ...
                                          'edgecolor', 'black');
                                 

                              
        % Now reset a bunch of stuff that was screwed by the bar plot call
        % Add title for top row
        if m==1
            title(ax(m,r), regnames{r}, 'fontsize', 14, ...
                                        'fontweight', 'bold');
        end
              
        set(ax(m,r), 'xlim', xlim, 'xtick', xtick, ...
                     'xticklabel', xtickstr, 'xgrid', 'on', ...
                     'ylim', [0,ymax], 'ytick', [0,ymax], ...
                     'yticklabel', '', 'ygrid', 'on', ...
                     'box', 'off');
        set(bh, 'edgecolor', 'black');
                                    
end

end

end





