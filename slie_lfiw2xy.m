function [x,y] = slie_lfiw2xy(lfiw, cvecs, pixsz, csegs)
% Function to calculate the pixel x/y coordinates of a SLIE corresponding
% to landfast ice widths along specified coast vectors
%
% USAGE
%    [x,y] = slie_lfiw2xy(lfiw, cvecs, pixsz, csegs)
%
% INPUTS
%    lfiw : Nv-by-1 array of landfast ice width measurements
%           for a single SLIE
%   cvecs : Nv-by-4 array of coat vector coordinates
%   pixsz : pixel size
%   csegs : coastal segmentation information 
%           (see slie_read_coastsegmentsfile.m)
%
% OUTPUTS
%     x,y : Nv-by-1 arrays of x- and y-coordinates for each point where
%           a coast vector intersects the SLIE
%
% Andy Mahoney - October 2023
%

% -----------------------------------------------------------------------
% Calculate azimuth of each coast vector ...
x0 = cvecs(:,1);
y0 = cvecs(:,2);
dx = cvecs(:,3) - x0;
dy = cvecs(:,4) - y0;
az = atan2(dy,dx);

% ... and determine x,y position at every width measurement
x = x0 + lfiw.*cos(az)/pixsz;
y = y0 + lfiw.*sin(az)/pixsz;


% Now handle any overlapping segments
Nseg = numel(csegs);
segdone = zeros(Nseg,1);
for s=1:Nseg


    if strcmp(csegs(s).overlap, "none") || segdone(s) == 1
        % Skip to next segment if this is not an overlapping segment
        % or this segment has aleady one has already been dealt with
        % (i.e. this segment overlaps with another)
        continue;
    end

    % Extract lfiw x/y points for this coastal segment
    if strcmp(csegs(s).direction, 'right')
        csi = csegs(s).v0:csegs(s).v1;
    else
        csi = csegs(s).v1:-1:csegs(s).v0;
    end
    seg.x = x(csi);
    seg.y = y(csi);
    
    % And do same for overlapping segment
    o = find(strcmp({csegs(:).segname}, csegs(s).overlap));
    if strcmp(csegs(o).direction, 'right')
        osi = csegs(o).v0:csegs(o).v1;
    else
        osi = csegs(o).v1:-1:csegs(o).v0;
    end
    oseg.x = x(osi);
    oseg.y = y(osi);

    % Remove any NaNs from segment coordinates
    % - not entirely sure how they would have go there
    %   but we can effectively interpolate over them 
    nani = isnan(seg.x) | isnan(seg.y);
    seg.x = seg.x(~nani);
    seg.y = seg.y(~nani);
    nani = isnan(oseg.x) | isnan(oseg.y);
    oseg.x = oseg.x(~nani);
    oseg.y = oseg.y(~nani);
    
    % Combine the overlapping SLIE segments
    combseg = slie_combineoverlappingSLIEs(seg, oseg, cvecs(csi,:), ...
                                                      cvecs(osi,:));


    % Reconstruct x- and y-coordinate arrays to replace both
    % segments with the combined segment
    % with those of combined SLIE
    presegi = 1:(min(csi)-1);
    intersegi = (max(csi)+1):(min(osi)-1);
    x = [x(presegi); combseg.x; x(intersegi)];
    y = [y(presegi); combseg.y; y(intersegi)];

    % Mark these segments as done
    segdone([s,o]) = 1;


end

end