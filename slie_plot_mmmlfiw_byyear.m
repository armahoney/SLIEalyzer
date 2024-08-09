function fig = slie_plot_mmmlfiw_byyear(dataset, varargin)

% Function to plot a timeseries of landfast ice width for a specified
% set of coast vectors
%
% USAGE
%     fig = slie_plot_lfiwatcvec(dataset)
%     fig = slie_plot_lfiwatcvec(__, 'cvrange', cvrange)
%     fig = slie_plot_lfiwatcvec(__, 'doy0', doy0)
%     fig = slie_plot_lfiwatcvec(__, 'doy1', doy1)
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
%         doy0 : Scalar integer specifying starting day of year to plot
%                Default: 274 (October 1)
%         doy1 : Scalar integer specifying ending day of year to plot
%                Default: 212 (July 31)
%     quantile : a scalar value between 0 and 0.5 expressing the upper and
%                lower bounds to use for plotting range of lfiw values
%                default: 0.1 (upper and lower 10%)
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
%        fig : figure handle for output plots
%
% Andy Mahoney - June 2024
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
        case 'cvrange'
            cvrange = varargin{a+1};
            a = a + 2;
        case 'doy0'
            doy0 = varargin{a+1};
            a = a + 2;
        case 'doy1'
            doy1 = varargin{a+1};
            a = a + 2;
        case 'quantile'
            quantile = varargin{a+1};
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

if ~exist('doy0', 'var')
    doy0 = 274;
end

if ~exist('doy1', 'var')
    doy1 = 212;
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


% Discard data outside doy bounds
% (by setting to NaNs)
[syr, smo, sda] = datevec(sdates_cv);
sdoy = dayofyear(sdates_cv);
if doy0 < doy1
    keepi = (sdoy >= doy0) & (sdoy <= doy1);
else
    keepi = (sdoy >= doy0) | (sdoy <= doy1);
end
lfiw_cv(~keepi) = NaN;
sdates_cv(~keepi) = NaN;


% Plot lfiw by season / year
Nseas = numel(cfg.seasons);
lfiw_median = zeros(Nseas,1);
lfiw_low = zeros(Nseas,1);
lfiw_upp = zeros(Nseas,1);
seasdoy0 = dayofyear(datenum('2001/10/01'));
seasdoy1 = dayofyear(datenum('2001/07/31'));
for s=1:Nseas
    seasyr0 = str2double(cfg.seasons{s}(1:4));
    seasyr1 = seasyr0 + 1;
    si = ((syr == seasyr0) & (sdoy >= seasdoy0)) | ...
         ((syr == seasyr1) & (sdoy <= seasdoy1));
    lfiw_seas = lfiw_cv(si);

    fini = ~isnan(lfiw_seas);
    lfiw_fin = lfiw_seas(fini);
    lfiw_median(s) = median(lfiw_fin);
    Nfin = sum(fini);
    lfiw_sort = sort(lfiw_fin);
    lowi = max(1, floor(qlow*Nfin));
    lfiw_low(s) = lfiw_sort(lowi);
    uppi = ceil(qupp*Nfin);
    lfiw_upp(s) = lfiw_sort(uppi);
end


% ------------------------------------------------------------------------
% Initialize figure
fig = figure('units', 'inches', 'position', [0,0,24,5], ...
                                'papersize', [24,5]);
xt = (1:Nseas)';
xlim = [0,Nseas]+0.5;
xtl = cfg.seasons;


ax = axes('parent', fig, 'pos', [0.05,0.1,0.93,0.85], ...
                         'xlim', xlim, 'xtick', xt, 'xticklabel', xtl, ...
                         'ylim', [0,130], 'ytick', 0:20:120, ...
                         'fontsize', 30);
ylabel(ax, 'km', 'fontsize', 30);



px = [xt; xt(end:-1:1)];
py = [lfiw_low; lfiw_upp(end:-1:1)];
pcol = [0.8,0.8,0.95]; % Pale blue ish

patch(px, py, pcol, 'parent', ax, 'edgecolor', 'none', ...
                                  'displayname', 'range');
line(xt, lfiw_low, 'parent', ax, 'color', 'blue', ...
                                 'displayname', 'low');
line(xt, lfiw_upp, 'parent', ax, 'color', 'blue', ...
                                 'displayname', 'upp');
line(xt, lfiw_median, 'parent', ax, 'color', 'black', ...
                                  'linewidth', 2, ...
                                  'displayname', 'median');

end










