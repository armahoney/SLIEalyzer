function [headbays, hbonames] = group_bays_headlands(varargin)
% Function to assign coast points into headland or bays based on
% how these features are encoded in the headbaymask file
%
% USAGE
%    headbays = group_bays_headlands
%    headbays = group_bays_headlands(__, 'dataset', dataset)
%    headbays = group_bays_headlands(__, 'cvecs', cvecs)
%    headbays = group_bays_headlands(__, 'headbaymask', headbaymask)
%    headbays = group_bays_headlands(__, 'lagoonmask', lagoonmask)
%    headbays = group_bays_headlands(__, 'display', display)
%   [headbays, hbonames] = group_bays_headlands(__)
%
% INPUTS
%     dataset : string specifying which region to use.
%               See slie_get_regionfolder.m for accepted values.
%               default: 'beau'
%       cvecs : cvecs : Nvec-by-4 array of coast vector coordinates
%               Default: read-in from file specified in config file
% headbaymask : full path to non-geoTIFF with headlands/bays encoded
%               Default: read-in from file specified in config file
%  lagoonmask : full path to non-geoTIFF with lagoons encoded
%               Default: read-in from file specified in config file
%     display : Scalar flag (0 or 1) specifying wither to display vectors 
%               assigned to each coastal type for diagnostics and 
%               use in figures
%               Default: 0 (no display)
%
% OUTPUTS
%    headbays : a 1-by-Nvecs array of values encoding the morphology type
%               for each coast vector, using the following values:
%               -   headlands : 1           
%               - open coasts : 2 (or 2.4 if lagoon also present)
%               -  embayments : 3 (or 3.4 if lagoon also present)
%               -     lagoons : 0.4 added to other morphology code
%    hbonames : Cell array containing names assigned below to 
%               each type of morphology:
%                {'Headlands', ['Open', 'coastlines'], ...
%                 'Embayments', 'Lagoons'};
%
% Andy Mahoney - October 2023
%

% ------------------------------------------------------------------------
% Assign values that are use to code features 
% in annotated bay/lagoon/headland images
headval = 18;
bayval = 54;
openval = 100;
lagval = 182;
landval = 255;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Deal with arguments specified in function call
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'cvecs'
            cvecs = varargin{a+1};
            a = a + 2;
        case 'headbaymask'
            headbaymask = varargin{a+1};
            a = a + 2;
        case 'lagoonmask'
            lagoonmask = varargin{a+1};
            a = a + 2;
        case 'region'
            region = varargin{a+1};
            a = a + 2;
        case 'display'
            displayvecs = 1;
            a = a + 1;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Figure out which dataset we're dealing with
% set the appropriate base path for the data
if ~exist('dataset', 'var')
    dataset = 'beau';
end

% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);


% Read the config file
cfg = read_SLIEconfigfile(sfolder);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Assign default values for arguments not specified in function call
% - read data from files specified in config file if necessary
if ~exist('displayvecs', 'var'), displayvecs = 0; end

if ~exist('cvecs', 'var')
    % Read cvecs file
    [cvecs, Nvecs]  = read_cvecsfile(cfg.cvecfile);
end

if ~exist('headbaymask', 'var')
    % Read headland / bay mask
    headbaymask  = imread(cfg.headbaymaskfile, 'tiff');
end
sz = size(headbaymask);

if ~exist('lagoonmask', 'var')
    % Read lagoon mask
    lagoonmask  = imread(cfg.lagoonmaskfile, 'tiff');
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Assign points along coast from cvecs array
cpx = cvecs(:,1); % Add 1 to convert from IDL zero-based indices
cpy = cvecs(:,2);
cpi = (sz(1)*(cpx-1)) + cpy;

% Get values at these points from headland/bay mask
cphb = headbaymask(cpi);
head_i = cphb == headval;
bay_i = cphb == bayval;
open_i = cphb == openval;

% And do the same from lagoon mask
cpl = lagoonmask(cpi);
lag_i= cpl == lagval;

% Create an array to code the type of coastline at each coast point
headbays = zeros(1,Nvecs);
headbays(head_i) = 1;           
headbays(open_i) = 2;
headbays(bay_i) = 3;
headbays(lag_i) = headbays(lag_i) + 0.4;

% Assign names to morphology-based reions
hbonames = {'Headlands', ['Open', 'coastlines'], 'Embayments', 'Lagoons'};
Nhbo = numel(hbonames);


% Display vectors assigned to each coastal type for 
% diagnostics and use in figures
if displayvecs
    regveccols = {'red','green','blue','magenta','cyan','yellow'};
    
    fig = figure;
    ax = axes('parent', fig);
    
    imagesc([1,sz(2)]-0.05,[1,sz(1)]-0.5, headbaymask, 'parent', ax);
    colormap(ax, 'gray');

    for r = 1:Nhbo
        if  ~strcmp(hbonames{r}, 'Lagoons')
            ri = find(floor(headbays) == r);
        else
            ri = find(round(10*(headbays-floor(headbays))) == r);
        end
        
        Nv = numel(ri);

        rx = zeros(3*Nv, 1);
        ry = zeros(3*Nv, 1);
        for v=1:Nv
            li = 3*(v-1)+(1:3);
            rx(li) = [cvecs(ri(v),1), cvecs(ri(v),3), NaN];
            ry(li) = [cvecs(ri(v),2), cvecs(ri(v),4), NaN];
        end
    
        line(rx, ry, 'parent', ax, 'color', regveccols{r}, ...
                     'displayname', hbonames{r});
    end
end


end
