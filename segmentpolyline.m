function S = segmentpolyline(L, maxsep)
% A fairly simple function to identify discontinuities
% in a polyline based on a maximum distance between vertices
% Returns a structure containing a list of x and y cooordinates
% with NaNs inserted between continuous segments
%
% USAGE
%   S = segmentpolyline(L)
%   S = segmentpolyline(L, maxsep)
%
% INPUTS
%        L : a structure containing fields x and y:
%            - x : m-by-1 vector of x coordinates
%            - y : m-by-1 vector of y coordinates
%   maxsep : maximum allowable separation between vertices of a
%            continuous line segment
%
% OUTPUT
%        S : a structure containing the following fields:
%            -    x : (m+s-1)-by-1 vector of x coordinates
%            -    y : (m+s-1)-by-1 vector of y coordinates
%            -    i : (m+s-1)-by-1 vector identifying the original, 
%                     unsegmented indices of each vertex in the segmented line
%            - seg0 : s-by-1 list of indices specifying 
%                     the start of each continous segment
%            - seg1 : s-by-1 list of indices specifying 
%                     the end of each continous segment

% ------------------------------------------------------------------------

% Set default value of maxsep if not spec'd in functuion call
if ~exist('maxsep', 'var')
    maxsep = 100;
end


% Calculate distance between each vertex in line
d = sqrt(diff(L.x).^2 + diff(L.y).^2);  

% Identify continuous segments between breaks
segbreak = find(d >= maxsep);
seg0 = [1; segbreak+1];
seg1 = [segbreak; numel(L.x)];

% Idenfity and exclude any single-vertex segments
singlev = seg0 == seg1;
seg0 = seg0(~singlev);
seg1 = seg1(~singlev);

% Initialize arrays to hold NaN-segmented coordinates
Nseg = numel(seg0);
Ns = sum(1+seg1-seg0) + Nseg-1;
S.x = zeros(Ns,1) + NaN;
S.y = zeros(Ns,1) + NaN;
S.i = zeros(Ns,1) + NaN;
S.seg0 = zeros(Nseg,1);
S.seg1 = zeros(Nseg,1);
% Construct line leaving gaps (i.e. NaNs) at segment breaks
j = 0;
for s=1:Nseg
    % Get indices for this segment from unsegmented line
    si = seg0(s):seg1(s);

    % Detemine corresponding indices in NaN-segmented line
    Nsi = numel(si);
    sj = j + (1:Nsi);

    % Assign points from unsegment line to NaN-segmented line
    S.x(sj) = L.x(si);
    S.y(sj) = L.y(si);
    S.i(sj) = si;
    S.seg0(s) = sj(1);
    S.seg1(s) = sj(Nsi);

    % Updated j index to leave a one-NaN-gap between segments
    % (except for last point in polyline)
    j = j + Nsi + (s<Nseg);
    notalot = 0;
end


end