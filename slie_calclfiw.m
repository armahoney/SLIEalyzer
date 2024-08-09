function [lfiw, wdep] = slie_calclfiw(slieimg, cvecs, cdist, varargin)
% Takes a SLIE image and a list of coast vectors
% and measures the width of the landfast ice along each vector
% also determines water depth at this point if given a bathymetry image

% USAGE:
%   lfiw = slie_calclfiw(slieimg, cvecs, cdist)
%   lfiw = slie_calclfiw(__, 'bathy', bathy)
%   lfiw = slie_calclfiw(__, 'pixsz', pixsz)
%   lfiw = slie_calclfiw(__, 'cvecmap', cvecmap)
%   lfiw = slie_calclfiw(__, 'diagnost', disgnost)
%  [lfiw, wdep] = slie_calclfiw(__)
%
% INPUTS
%     slieimg : 2-D matrix representing a raster binarization of 
%               landfast ice inside of the seaward landfast ice edge (SLIE)
%       cvecs : Nv-by-4 array of coordinates defining beginning and end
%               of the coast vectors along which landfast ice width is
%               measured
%       cdist : 2-D matrix of same dimensions as slieimg specifying the 
%               distance from the coast every non-land pixel          
%       bathy : 2-D matrix of same dimensions as slieimg specifying the 
%               water depth at every non-land pixel          
%       pixsz : size of each pixel in slieimg
%               Default: 100 m (pixels are assumed to be square)
%     cvecmap : ...
%               Default: ...
%    diagnost : A scalar flag specifying whether to plot diagnostic
%               information useful for debugging
%               Default: 0 (don't plot diagnostics)
%
%
% OUTPUTS:
%   lfiw : A Nv-by-1 array of landfast ice widths, 
%          where Nv = number of coast vectors
%   wdep : A Nv-by-1 array of water depths at the SLIE
%          if a bathymetry array is provided with the function call
%
%
% Andy Mahoney October 2023


% -----------------------------------------------------------------------
% Handle variable arguments specified in function call
Nvararg = numel(varargin);
a = 1;
while (a <= Nvararg)
    switch (varargin{a})
        case 'pixsz'
            pixsz = varargin{a+1};
            a = a + 2;
        case 'bathy'
            bathy = varargin{a+1};
            a = a + 2;
        case 'cvecmap'
            cvecmap = varargin{a+1};
            a = a + 2;  
        case 'diagnost'
            diagnost = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a + 1;
    end
end

% Set default values of not spec'd in function call
if ~exist('pixsz', 'var')
    pixsz = 0.1;
end
if ~exist('diagnost', 'var')
    diagnost = 0;
end

% ------------------------------------------------------------------------
% Determin various dimensions from input arguments
sz = size(slieimg);
numvecs = numel(cvecs)/4;

% Pre-assign array to hold lfiw data
lfiw = zeros(numvecs,1);

% Pre-assign array to hold water depths if bathymetry is provided
if exist('bathy', 'var')
    wdep = zeros(numvecs,1);
end

% ========================================================================
% NOW GET TO BUSINES

% Go through each coast vector and determine intersection point with SLIE 
for v =1:numvecs
    width = 0.0;
    wdepth = 0.0;
    
    % Store endpoints of current coast vector in more convenient variables
    x0 = cvecs(v,1);
    y0 = cvecs(v,2);
    x1 = cvecs(v,3);
    y1 = cvecs(v,4);
    
    % Determine x-y coords of all points along coast vector
    % and extract SLIE image data at those points
    txy = slie_get_tsectxy(x0,y0,x1,y1, sz);
    tsect = slieimg(txy.i);
    tlen = numel(tsect);
    
    % Missing SAR mosaic regions in the M2014 dataset are coded with 111
    nodatacheck = sum(tsect == 111);
    if (nodatacheck == 0)
        
        % If there are islands in the landmask, some transects will intersect 
        % land pixels. If this happens, we fill in the landvalues based on
        % the value of the first pixel seaward of each island
        landpix = tsect == 128;
        if any(landpix)
            isl0 = find([landpix(1), diff(landpix)] == 1);
            isl1 = find([diff(landpix), -landpix(end)] == -1);
            Nisl = numel(isl0);

            % Diagnostic check for wierd island stuff
            if numel(isl1) ~= Nisl
                fprintf('Something wierd island stuff going on in cvec %0.0f\n', v);
                lfiw = NaN;
                wdep = NaN;
                notalot = 0;
                return
            end
            
            %islands_i = find(landpix);
            for i=1:Nisl
                     tsect(isl0(i):isl1(i)) = tsect(isl1(i)+1);
            end
            
        end
        
        % % Find landfast ice and non-landfast ice pixels along transect
        lfi = 1*(tsect == 255);
        Nlfi = sum(lfi);
        
        % Ideally, we would just find first non-landfast ice pixel
        % And set the edge one pixel shoreward of this
        
        % However, sometimes there's some non landfast ice near the coast, 
        % which we want to ignore (e.g. due to imperfectly drawn polygons
        % in ice chart data)
        % Also, the coast vectors occasionally intersect landfast ice again
        % at their shorward end, which creates a similar phenomena that
        % that we don't want to include in the lfiw measurement.
        % This tends to happen around Kotzebue Sound if it happens
        
        % So, if there's more than one SLIE, choose the one that is 
        % furthest from land, as determined from the cdistmap array
        % NOTE: this is not necessarily the furthest from the origin of 
        %       the coast vector in question        
        switch true
            case (Nlfi == 0)
                % No landfast ice: set slie_i to 0 and go on to next
                slie_i = 1;
            case (Nlfi == tlen)
                %  No end to landfast ice: set slie_i to transect length  
                slie_i = tlen;
            otherwise
                % Now find SLIEs and anti-SLIEs
                lfidif = [diff(lfi), -lfi(end)];
                slies = find(lfidif == -1);
                antislies = find(lfidif == 1);

                Nslies = numel(slies);
                
                % If there's only onne SLIE, it's pretty easy ...
                if Nslies == 1
                    slie_i = slies(1);
                else
                    % But if there's more than one SLIE, this might mean
                    % the coastline is a bit more complex, so we need a
                    % more sophisticated way of picking the right edge
                    % e.g., in Kotzebue Sound, the coastvector might
                    %       interesect fast ice attached to a different
                    %       section of the coastline
                    
                    % We should try to pick the edge that's farthest from
                    % land, which isn't necessarily the one that's farthest
                    % along the coast vector:
                    % ** in places like Kotzebue Sound, coastal distance
                    %    doesn't necessarily increase monotonically
                    %    along a coast vector
                    
                    % Find distance from coast along transect 
                    cd_tsect = cdist(txy.i);

                    % Zero-out regions of tsect where (smoothed) coastal 
                    % distance isn't increasing
                    % (i.e. tsect is heading toward land again)
                    smwid = min([100,round(tlen/2)]);
                    smkern = ones(1,smwid)/smwid;
                    cd_tsect_sm = conv(cd_tsect, smkern, 'same');
                    cd_trend = diff([0, cd_tsect_sm]);
                    cd_incr = cd_tsect_sm.*(cd_trend >= 0);
                    slie_cd = cd_incr(slies);
                                       
                    % Find the SLIE that is farthest from the coast
                    % *** and where the distance from coast is increasing ***
                    [~, farthest_i] = max(slie_cd);
                    slie_i = slies(farthest_i);         
                end
        end
        
        
        % Convert SLIE position to landfast ice width
        % Note: this is NOT the same as location of the edge along the
        %       coast vector transect, since this has duplicate coordinates
        %       removed (See slie_get_tsectxy.m)
        lfiw(v) = round(txy.d(slie_i))*pixsz;
        if exist('bathy', 'var')
            if slie_i > 1
                wdep(v) = -bathy(txy.y(slie_i), txy.x(slie_i));
            else
                wdep(v) = 0;
            end
        end

    else 
        lfiw(v) = NaN;
        if exist('bathy', 'var')
            wdep(v) = NaN;
        end
    end
                        
end


end

