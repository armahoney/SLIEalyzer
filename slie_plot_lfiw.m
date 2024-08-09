function fig = slie_plot_lfiw(dataset, varargin)

% Function to plot a timeseries of landfast ice width for a specified
% set of coast vectors
%
% USAGE
%     fig = slie_plot_lfiwatcvec(dataset)
%     fig = slie_plot_lfiwatcvec(__, 'cvrange', cvrange)
%     fig = slie_plot_lfiwatcvec(__, 'yr0', yr0)
%     fig = slie_plot_lfiwatcvec(__, 'yr1', yr1)
%     fig = slie_plot_lfiwatcvec(__, 'foldbyseason', foldbyseason)
%     fig = slie_plot_lfiwatcvec(__, 'quantile', quantile)
%     fig = slie_plot_lfiwatcvec(__, 'asadir', asadir)
%     fig = slie_plot_lfiwatcvec(__, 'cfgfile', cfgfile)
%
% INPUTS
%      dataset : string specifying which dataset to use.
%                See slie_get_basefolder.m for accepted values.
%                default: 'beau'
%    subregion : a string corresponding to the name of a subregion
%                specified in the cvecregs file
%                Default: '' (blank)
%                default: 1:10 (first 10 coast vectors)
%      cvrange : Either: 
%                -  a Nv-by-1 vector specifying indices for coast vectors
%                   to be plotted; or
%                -  a character string specifying either "all" or the name
%                   of a subregion specified in subregfile
%                default: 1:10 (first 10 coast vectors)
%          yr0 : Scalar integer specifying starting year to plot
%                Default: first year in record
%          yr1 : Scalar integer specifying ending year to plot
%                Default: last year in record
%     quantile : a scalar value between 0 and 0.5 expressing the upper and
%                lower bounds to use for plotting range of lfiw values
%                default: 0.1 (upper and lower 10%)
% foldbyseason : a scalar flag specifying whether to plot values as one
%                long timeseries or to fold the data by season
%                default: 0 (do not fold by season)
%      cfgfile : string specifying the name of the configuraion file stored
%                in the region folder containing essential information for
%                performing all analysis of SLIE data. 
%                See slie_readconfigfile.m for more details
%                default: 'slieconfig.txt'
%       asadir : string specifying folder name in region folder where
%                the data files from the all-seasons analysis will be saved
%                default: 'AllSeasonsAnalysis2023'
%
% OUTPUTS
%        ax : axes handles for output plots
%
% Andy Mahoney - January 2024
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
        case 'subregion'
            subregion = varargin{a+1};
            a = a + 2;
        case 'cvrange'
            cvrange = varargin{a+1};
            a = a + 2;
        case 'yr0'
            yr0 = varargin{a+1};
            a = a + 2;
        case 'yr1'
            yr1 = varargin{a+1};
            a = a + 2;
        case 'quantile'
            quantile = varargin{a+1};
            a = a + 2;
        case 'foldbyseason'
            foldbyseason = varargin{a+1};
            a = a + 2;
        case 'asadir'
            asadir = varargin{a+1};
            a = a + 2;
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
   

if ~exist('cvrange', 'var')
    cvrange = 1:10;
end

if ~exist('quantile', 'var')
    quantile = 0.1;
end

if ~exist('foldbyseason', 'var')
    foldbyseason = 0;
end

if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end



% ------------------------------------------------------------------------
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

% Read config info
cfg = slie_readconfigfile(sfolder, cfgfile);

% Set upper and lower bounds
qlow = quantile;
qupp = 1-quantile;

% ------------------------------------------------------------------------
% Read data from All Seasons Analysis folder
lfiwfile = [sfolder filesep asadir filesep 'lfiw_all.csv'];
[lfiw, slienames] = slie_readsliedata(lfiwfile);
sdatefile = [sfolder filesep asadir filesep 'sliedatenum_all.csv'];
sdates = slie_readsliedata(sdatefile);
[syr, smo, sda] = datevec(sdates);
sdoy = dayofyear(sdates);

% Assign start and end years if not specified in function call
if ~exist('yr0', 'var')
    yr0 = min(syr(:));
end
if ~exist('yr1', 'var')
    yr1 = max(syr(:));
end

% Discard data outside year bounds
% (by setting to NaNs)
oct1 = dayofyear(datenum('1999/09/01'));
ltyr0 = (syr < yr0) | ((syr == yr0) & sdoy < oct1);
gtyr1 = (syr > yr1) | ((syr == yr1) & sdoy >= oct1);
lfiw(ltyr0 | gtyr1) = NaN;
sdates(ltyr0 | gtyr1) = NaN;

% Select coast vector range
if ischar(cvrange)
    if strcmpi(cvrange, 'all')
        cvrange = 1:size(lfiw, 1);
    else
        [~, subregnames, subregbounds] = ...
                          slie_readcvecregsfile(cfg.cvecregionsfile);
        regi = strcmp(subregnames, cvrange);
        if any(regi)
            cv1 = subregbounds(regi, 1);
            cv2 = subregbounds(regi, 2);
            cvrange = cv1:cv2;
        else
            disp(['cvrange value not recognized: ' cvrange]);
            return
        end
    end
end
lfiw_cv = lfiw(cvrange, :);
sdates_cv = sdates(cvrange, :);

% Fold by season if specified
if foldbyseason
    oct1 = dayofyear(datenum('1999/10/01'));
    sdoy2 = dayof2year(sdates_cv, oct1);
    Nday = 303;
    dd = oct1 + (1:Nday)' -1;
    lfiw_mean = zeros(Nday,1);
    lfiw_low = zeros(Nday,1);
    lfiw_upp = zeros(Nday,1);
    for d=1:Nday
        di = sdoy2 == (oct1 + d - 1);
        if any(di(:))
            lfiwd = lfiw_cv(di);
            fini = ~isnan(lfiwd);
            lfiw_fin = lfiwd(fini);
            lfiw_mean(d) = median(lfiw_fin);
            Nfin = sum(fini);
            lfiw_sort = sort(lfiw_fin);
            lowi = max(1, floor(qlow*Nfin));
            lfiw_low(d) = lfiw_sort(lowi);
            uppi = ceil(qupp*Nfin);
            lfiw_upp(d) = lfiw_sort(uppi);
        end
        notalot = 0;
    end
    notalot = 0;
else
    sdates_fin = sdates_cv(isfinite(sdates_cv));
    dd = sort(unique(sdates_fin));
    Nday = numel(dd);
    lfiw_mean = zeros(Nday,1);
    lfiw_low = zeros(Nday,1);
    lfiw_upp = zeros(Nday,1);
    for d=1:Nday
        di = sdates_cv == dd(d);
        lfiwd = lfiw_cv(di);
        fini = ~isnan(lfiwd);
        lfiw_fin = lfiwd(fini);
        lfiw_mean(d) = median(lfiw_fin);
        Nfin = sum(fini);
        if Nfin > 0
            lfiw_sort = sort(lfiw_fin);
            lowi = max(1, floor(qlow*Nfin));
            lfiw_low(d) = lfiw_sort(lowi);
            uppi = ceil(qupp*Nfin);
            lfiw_upp(d) = lfiw_sort(uppi);
        else
            lfiw_low(d) = NaN;
            lfiw_upp(d) = NaN;
        end
        disp([d, Nday]);
        notalot = 0;
    end
end

% ------------------------------------------------------------------------
% Initialize figure
fig = figure('units', 'inches', 'position', [0,0,24,5], ...
                                'papersize', [24,5]);
if foldbyseason
    xlim = [oct1, oct1 + Nday - 1];
    xt = datenum([zeros(11,1), (1:11)'+9, ones(11,1), zeros(11,3)]);
    xtmid = (xt(1:end-1) + xt(2:end))/2;
    xtl = datestr(xtmid, 'mmm');
else
    [yr0,~,~] = datevec(dd(1));
    [yr1,~,~] = datevec(dd(end));
    xlim = datenum([[yr0;yr1+1],[1;1],[1;1],zeros(2,3)]);
    Nyrs = 1+yr1-yr0;
    xt = datenum([(yr0:(yr1+1))', ...
                  ones(Nyrs+1,1), ...
                  ones(Nyrs+1,1), ...
                  zeros(Nyrs+1,3)]);
    xtmid = (xt(1:end-1) + xt(2:end))/2;
    xtl = datestr(xtmid, 'yyyy');
end

ax = axes('parent', fig, 'pos', [0.05,0.1,0.93,0.85], ...
                         'xlim', xlim, 'xtick', xt, 'xticklabel', [], ...
                         'ylim', [0,130], 'ytick', 0:20:120, ...
                         'fontsize', 30);
ylabel(ax, 'km', 'fontsize', 30);


text(xtmid, zeros(numel(xt)-1,1), xtl, 'parent', ax, ...
                                     'fontsize', 20, ...
                                     'horizontalal', 'center', ...
                                     'verticalal', 'top', ...
                                     'units', 'data');


px = [dd; dd(end:-1:1)];
py = [lfiw_low; lfiw_upp(end:-1:1)];
pcol = [0.8,0.8,0.95]; % Pale blue ish

patch(px, py, pcol, 'parent', ax, 'edgecolor', 'none', ...
                                  'displayname', 'range');
line(dd, lfiw_low, 'parent', ax, 'color', 'blue', ...
                                 'displayname', 'low');
line(dd, lfiw_upp, 'parent', ax, 'color', 'blue', ...
                                 'displayname', 'upp');
line(dd, lfiw_mean, 'parent', ax, 'color', 'black', ...
                                  'linewidth', 2, ...
                                  'displayname', 'median');

end










