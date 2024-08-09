function C = slie_combineoverlappingSLIEs(slieA, slieB, cvecA, cvecB)
% A function to find the best SLIE to use when there are different
% sets of coast vectors measuring landfast ice width along the same
% section of coast - such as those to measure lfiw in Inner and Outer
% Kotzebue Sound
%
% Breifly the function uses both sets of coast vectors to find
% intersections with both SLIEs. Where there are intersections with
% both SLIEs, the mid-point is taken. 
% This results in the SLIEs converging for both sets of coast vectors.
% We then select the SLIE that covers the greatest length and is 
% defined by the greatest number of unique vertices / intersections
%
%
% USAGE
%   C = slie_combineoverlappingSLIEs(slieA, slieB, cvecA, cvecB)
%
% INPUTS
%   slieA, slieB : structures containing x- and y-fields defining
%                  the two SLIEs
%   cvecA, cvecB : Nvec-by-4 arrays specifiying end-points of the
%                  two sets of coast vectors
%
% OUTPUT
%              C : a structure containing x- and y-fields defining
%                  the "best" SLIE after combining
%
% Andy Mahoney - October 2023

% -------------------------------------------------------------------

% Identify gaps / discontinuities in each SLIE based on max separation
% and insert NaNs and exclude single vertices 
% for purposes of finding intersections
A = segmentpolyline(slieA);
B = segmentpolyline(slieB);
Aii = find(~isnan(A.x));
Bii = find(~isnan(B.x));

% Diagnostics
% plot(cvecA(:,1), cvecA(:,2), 'color', [0.75,0.5,0.5]);
% set(gca, 'ydir', 'reverse');
% line(cvecB(:,1), cvecB(:,2), 'color', [0.5,0.5,0.75]);
% line(A.x, A.y, 'color', 'red', 'marker', '.');
% line(B.x, B.y, 'color', 'blue', 'marker', '.');

% ------------------------------------------------------------------------
% Go through both sets of coast vectors and:
% - find intersections with both SLIEs
% - take mid-point whereever there are 2 intersections
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Go through coast vectors for SLIE A 
Nva = size(cvecA, 1);
A2.x = zeros(Nva,1) + NaN;
A2.y = zeros(Nva,1) + NaN;
for v=1:Nva
    % Determine intersection with SLIE A
    [Aix, Aiy, Aw] = getSLIEintersection(A, cvecA(v,:));

    % Determine intersection with SLIE B
    [Bix, Biy, Bw] = getSLIEintersection(B, cvecA(v,:));

   
    % If there were intersections with both SLIEs:
    % - of one of them has zero width, keep that one
    % - otherwise take the average position
    if any(Aix) && any(Bix)
        if (Aw == 0)
            A2.x(v) = Aix;
            A2.y(v) = Aiy;
        elseif (Bw == 0)
            A2.x(v) = Bix;
            A2.y(v) = Biy;
        else
            %Av.x(v) = Aix*(Aw>=Bw) + Bix*(Bw>Aw);
            %Av.y(v) = Aiy*(Aw>=Bw) + Biy*(Bw>Aw);
            A2.x(v) = (Aix + Bix)/2;
            A2.y(v) = (Aiy + Biy)/2;
        end
    elseif any(Aix)
        A2.x(v) = Aix;
        A2.y(v) = Aiy;
    elseif any(Bix)
        A2.x(v) = Bix;
        A2.y(v) = Biy;
    end

    if A2.x(v) == 0
        wtf = 99;
    end

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Go through coast vectors for SLIE B 
Nvb = size(cvecB, 1);
B2.x = zeros(Nvb,1) + NaN;
B2.y = zeros(Nvb,1) + NaN;
for v=1:Nvb
    % Determine intersection with SLIE A
    [Aix, Aiy, Aw] = getSLIEintersection(A, cvecB(v,:));

    % Determine intersection with SLIE B
    [Bix, Biy, Bw] = getSLIEintersection(B, cvecB(v,:));

    % If there were intersections with both SLIEs:
    % - of one of them has zero width, keep that one
    % - otherwise take the average position
    if any(Aix) && any(Bix)
        if (Aw == 0)
            B2.x(v) = Aix;
            B2.y(v) = Aiy;
        elseif (Bw == 0)
            B2.x(v) = Aix;
            B2.y(v) = Aiy;
        else                
            %Bv.x(v) = Aix*(Aw>=Bw) + Bix*(Bw>Aw);
            %Bv.y(v) = Aiy*(Aw>=Bw) + Biy*(Bw>Aw);
            B2.x(v) = (Aix + Bix)/2;
            B2.y(v) = (Aiy + Biy)/2;
        end
    elseif any(Aix)
        B2.x(v) = Aix;
        B2.y(v) = Aiy;
    elseif any(Bix)
        B2.x(v) = Bix;
        B2.y(v) = Biy;
    end
end

% Diagnostics
%line(A2.x, A2.y, 'color', 'magenta', 'marker', '.');
%line(B2.x, B2.y, 'color', 'cyan', 'marker', '.');

    
% Pick the best one
finiA = ~isnan(A2.x);
[~,uiA] = unique([A2.x, A2.y], 'rows');
Ascore = sum(sqrt(diff(A2.x(finiA)).^2 + diff(A2.y(finiA)).^2))*numel(uiA); 
finiB = ~isnan(B2.x);
[~,uiB] = unique([B2.x, B2.y], 'rows');
Bscore = sum(sqrt(diff(B2.x(finiB)).^2 + diff(B2.y(finiB)).^2))*numel(uiB); 

if Ascore > Bscore
    C.x = round(medianfilter(A2.x(finiA), 10));
    C.y = round(medianfilter(A2.y(finiA), 10));
else
    C.x = round(medianfilter(B2.x(finiB), 10));
    C.y = round(medianfilter(B2.y(finiB), 10));
end

notalot = 0;

end






