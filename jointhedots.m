function [jx,jy] = jointhedots(x,y,varargin)
% A function to connect a series of pixel coordinates with the coordinates
% of the contiguous pixels in between them.
%
% NOTE: For these purposes pixel coordinates are vertices 
%       specified by integer x,y coordinates

% USAGE
%    [jx,jy] = jointhedots(x,y)
%    [jx,jy] = jointhedots(__, 'maxsep', maxsep)
%    [jx,jy] = jointhedots(__, 'minseg', minseg)
%
% INPUTS
%       x : m-by-1 vector of x coordinates
%       y : m-by-1 vector of y coordinates
%  maxsep : maximum separation between adjacent points to be considered
%           part of a continuous segment. Gaps will be left between
%           points separated by more than maxsep pixels
%           Default: Infinity (all points will be joined)
%  minseg : minimum segment length considered for joining. Segments defined
%           by fewer points will be exluded.
%           Default: 0 (all points will be joined)
%
% OUTPUTS
%     jx : M-by-1 vector of x-coordinate for joined-up points
%     jy : M-by-1 vector of y-coordinate for joined-up points

% -------------------------------------------------------------------------
% Deal with variable arguments
a = 1;
while a < numel(varargin)
    switch varargin{a}
        case 'maxsep'
            maxsep = varargin{a+1};
            a = a + 2;
        case 'minseg'
            minseg = varargin{a+1};
            a = a + 2;
        otherwise
            disp('jointhedots: Argument not recognized:');
            disp(varargin{a});
            a = a + 1;
    end
end

% Set infinite maximum separation if not specified
if ~exist('maxsep', 'var')
    maxsep = Inf;
end

% Set zero minimum segment length if not specified
if ~exist('minseg', 'var')
    minseg = 0;
end


% -------------------------------------------------------------------------
% Check number of vertices specified in function call
Nx = numel(x);
if numel(y) ~= Nx
    disp('jointhedots: x and y must be equal size.')
    yx = [];
    yx = [];
    return
end

% -------------------------------------------------------------------------
% Calculate distance between adjacent pairs of points
% This allows us to:
% - identify separate segments
% - identify and exclude duplicate points
% - estimate required number of joined-up points
dx = diff(x);
dy = diff(y);
dd = sqrt(dx.^2 + dy.^2);

% Identify unique points
ui = find(dd > 0);
if numel(ui) < Nx
    x = x(ui);
    y = y(ui);
    dx = diff(x);
    dy = diff(y);
    dd = sqrt(dx.^2 + dy.^2);
    Nx = numel(x);
end

% Identify continuous segments based on max separation
segbreak = find(dd >= maxsep);
Nseg = numel(segbreak)+1;
seg0 = [1; segbreak+1];
seg1 = [segbreak; Nx];

% (Over) Estimate number of points it will take to join the dots
Nj_guess = 2*round(Nx*ceil(mean(dd(dd < maxsep))));

% Pre-assign output arrays
jx = zeros(Nj_guess,1);
jy = zeros(Nj_guess,1);


% -------------------------------------------------------------------------
% Start connecting!
j = 0;  % Keep count of total number of joined up points
for s=1:Nseg

    % Skip this segment if it's too short
    if (seg1(s) - seg0(s) <= minseg)
        continue
    end
    
    % Set coordinate of first point in segment
    i = seg0(s);
    j = j + 1;
    jx(j) = x(i);
    jy(j) = y(i);

    % Go through and connect remaining points in segment
    while i < seg1(s)
        % Get x/y coordinates of pixels between vertices
        txy = slie_get_tsectxy(x(i), y(i), x(i+1), y(i+1));
        Nt = numel(txy.x);
        jj = j + (1:Nt);
        jx(jj) = txy.x;
        jy(jj) = txy.y;
    
        % Increment i and j counters
        i = i + 1;
        j = j + Nt;
    end

    notalot = 0;
end

% Trim output arrays to length
jx = jx(1:j);
jy = jy(1:j);


end