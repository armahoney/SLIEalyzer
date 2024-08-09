function [ix, iy, iw] = getSLIEintersection(slie, cvec)
% Function to determine the intersection point of a coast vector
% with a SLIE. 
%
% USAGE:
%     [ix, iy] = getSLIEintersection(slie, cvec)
%
% INPUTS
%     slie : A : a structure containing x and y fields defining the 
%            NaN-segmented SLIE
%     cvec : a 4-element array providing the start and end coordiantes
%            of the coast vector: [x0,y0,x1,y1[
%
% OUTPUTS
%       ix : a scalar value indicating the x-coordinate of the intersection
%       iy : a scalar value indicating the x-coordinate of the intersection
%       iw : a scalar value indicating the SLIE width at intersection
%
%

% -------------------------------------------------------------------------

% FIRST CHECK TO SEE IF SLIE LIES AT THE ORIGIN OF THE COAST VECTOR
if any((cvec(1) == slie.x) & (cvec(2) == slie.y))
    ix = cvec(1);
    iy = cvec(2);
    iw = 0;
    return;
end


% Calculate azimuth of coastvec
vaz = atan2(cvec(4)-cvec(2),cvec(3)-cvec(1));


% Determine location of any intesection between cvec and SLIE
[ix, iy, ii] = polyxpoly(slie.x, slie.y, cvec([1,3]), cvec([2,4]));
if any(ix)

    % There may multiple intersections, in which case
    % - pick only intersections coming from the left side of the SLIE
    % - pick the nearest of these
    % Calculcate azimuth of SLIE at intsersction point
    si = ii(:,1);
    dx = slie.x(si+1) - slie.x(si);
    dy = slie.y(si+1) - slie.y(si);
    saz = atan2(dy,dx);

    % Calculate left-hand component of the coast vector
    % when rotated to SLIE azimuth at intersection point
    leftcomponent = cos(saz)*sin(vaz) - sin(saz)*cos(vaz);

    % Keep only those with a positive left component
    keepi = leftcomponent < 0;
    ix = ix(keepi);
    iy = iy(keepi);
    ii = ii(keepi);
    
    % Calculate distance along coast vector for each intersection
    % and keep the furthest
    dx = ix - cvec(1);
    dy = iy - cvec(2);
    d = sqrt(dx.^2 + dy.^2);
    [iw, furthest] = max(d);
    ix = ix(furthest);
    iy = iy(furthest);
else
    iw = [];
end

end


