function t =  slie_get_tsectxy(x0,y0,x1,y1,sz)
% Function to calculate x,y coordinates for each pixel along transect 
% for a given pair of start/end coorindates

% USAGE:
%    t =  slie_get_tsectxy(x0,y0,x1,y1)
%    t =  slie_get_tsectxy(__, sz)
%
% INPUTS
%    x0 : x-coordinate of origin of transect
%    y0 : y-coordinate of origin of transect
%    x1 : x-coordinate of end of transect
%    y1 : y-coordinate of end of transect
%    sz : [optional] size of array in which transect is being taken
%
% OUTPUTS
%     t : a structure with the following fields:
%       - x: x-coordinates of unique points along transect
%       - y: y-coordinates of unique points along transect
%       - d: cartesian distance along transect of each point along transect
%       - i: (if sz is specified in function call) the 1-D index for each
%            point along transect



% Calculate azimuth and length of transect
dx = x1-x0;
dy = y1-y0;
mag = round(sqrt(dx^2 + dy^2));
az = atan2(dy, dx);
r = 0:mag;

% Define x,y coords along transect
t.x = round(x0 + r*cos(az));
t.y = round(y0 + r*sin(az));

% Calculate distance along transect in pixels
t.d = sqrt((t.x-t.x(1)).^2 + (t.y-t.y(1)).^2);


% Check for duplicates caused by rounding
dupi = (diff(t.x) == 0) & (diff(t.y) == 0);
if any(dupi)
    keepi = ones(size(r));
    keepi(dupi) = 0; % This will keep the first of each series of duplicates
    t.x = t.x(keepi == 1);
    t.y = t.y(keepi == 1);
    dupf = find(dupi);
    t.d(dupi) = mean([t.d(dupf); t.d(dupf+1)], 1);
    t.d = t.d(keepi == 1);
    notalot = 0;
end

% If array size is specified in function call:
if exist('sz', 'var')

    % Exclude any points outside specified limits
    inbounds = (t.x > 0) & (t.x <= sz(2)) & ...
               (t.y > 0) & (t.y <= sz(1));
    t.x = t.x(inbounds);
    t.y = t.y(inbounds);

    % Calculate 1-D index of transect points
    t.i = sub2ind(sz, t.y, t.x);
end   

end