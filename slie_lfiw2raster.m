function slieimg = slie_lfiw2raster(lfiw, cvecs, lmask, varargin)
% Function to draw and rasterize SLIE edge from lfiw data points


a = 1;
while a < numel(varargin)
    switch varargin{a}
        case 'pixsz'
            pixsz = varargin{a+1};
            a = a + 2;
        case 'csegs'
            csegs = varargin{a+1};
            a = a + 2;
        case 'shadregs'
            shadregs = varargin{a+1};
            a = a + 2;
        case 'noshadows'
            noshadows = varargin{a+1};
            a = a + 2;
        case 'filtwid'
            filtwid = varargin{a+1};
            a = a + 2;
        otherwise
            disp('Argument not recognized:');
            disp(varargin{a});
            a = a + 1;
    end
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set defaults where appropriate
% Default pixel size
if ~exist('pixsz', 'var')
    pixsz = 100;
end

% Close polygons using bottom left/right corners
% And fill from a random point on land
if ~exist('csegs', 'var')
    sz = size(lmask);
    csegs.segname = 'Default';
    csegs.v0 = 1;
    csegs.v1 = size(cvecs, 1);
    csegs.xclose = [sz(2), 1];
    csegs.yclose = [sz(1), sz(1)];
    foundland = 0;
    while (foundland == 0)
        randx = 2 + round(rand(1)*(sz(2)-2));
        randy = 2 + round(rand(1)*(sz(1)-2));        
        foundland = ~any(lmask(randy+(-1:1), randx+(-1:1)) < 255);
    end
    csegs.fillxy = [randx,randy];
    csegs.overlap = 'none';
end

% Add shadows by default
if ~exist('noshadows', 'var')
    noshadows = 0;
end

% Apply no median filtering by default
if ~exist('filtwid', 'var')
    filtwid = 0;
end

% Get size of image from coast mask
sz = size(lmask);

% Generate shadows region properties from coast mask
if ~exist('shadregs', 'var') && (noshadows == 0)
    cvshadows = lmask == 64;
    shadregs = regionprops(cvshadows, 'boundingbox', 'pixelidxlist', 'image');
    Nshad = numel(shadregs);
    for sh = 1:Nshad
        % Place feature image into large image 
        % with extra pixel all the way around
        shadimgX = zeros(size(shadregs(sh).Image)+2);
        shadimgX(2:end-1,2:end-1) = shadregs(sh).Image;
    
        % Grow object image with a 3x3 convolution
        shadconv = conv2(shadimgX, ones(3), 'same');
    
        % Identify the outer perimeter of feature
        outerperim = (shadconv > 0) & (shadimgX == 0);
        opi = find(outerperim);
        [opy,opx] = ind2sub(size(shadconv), opi);
    
        % Translate sub-image pixel index values into whole-image values
        % checking for oub-of-bounds (oob) pixels
        x0 = floor(shadregs(sh).BoundingBox(1));
        y0 = floor(shadregs(sh).BoundingBox(2));
        opxx = opx+x0-1;
        opyy = opy+y0-1;
        oob = (opxx < 1) | (opxx > sz(2)) | (opyy < 1) | opyy > sz(1);
        opxx = opxx(~oob);
        opyy = opyy(~oob);
        shadregs(sh).adjpix = sub2ind(sz, opyy, opxx);
        notalot = 0;
    end
end

% Force shadregs to empty array if noshadows is specified
if (noshadows == 1)
    shadregs = [];
end

% Initialize array to hold SLIE image
slieoutline = zeros(sz);
slieimg = zeros(sz);



% Intialize arrays to hold SLIE coordinates
% at intersections with coast vectors
Nv = numel(cvecs)/4;
sx = zeros(Nv,1);
sy = zeros(Nv,1);

% Calculate x,y pixel coordinates for landfast ice width measurement
% i.e., determine the point in image space where 
%       each coast vectors intersects the SLIE
for v=1:Nv
    % Figure out length and direction of coast vector
    dx = cvecs(v,3) - cvecs(v,1);
    dy = cvecs(v,4) - cvecs(v,2);
    az = atan2(dy,dx);
    d = sqrt(dx^2 + dy^2);
    
    % Calculate coordinates of each pixel along vector
    r = 0:round(d);
    vx = cvecs(v,1) + round(r*cos(az));
    vy = cvecs(v,2) + round(r*sin(az));

    % Find pixel corresponding SLIE intersection along vector
    lfi = find(r <= lfiw(v)/pixsz, 1, 'last');
    
    if any(lfi)        
        sx(v) = vx(lfi);
        sy(v) = vy(lfi);
    else
        sx(v) = NaN;
        sy(v) = NaN;
    end
    
    if (sy(v) < 1) || (sy(v) > sz(1)) || (sx(v) < 1) || (sx(v) > sz(2))
       wtf = 1;
    end
    
    notalot = 0;
end


% Now go through each coastal segment and:
% - add pre-defined points to close each segment
% - interpolate between LFIW points to create unbroken SLIE line
% - use flood fill operation to fill-in landfast ice area
Nseg = numel(csegs);
segdone = zeros(Nseg,1);
for s=1:Nseg

    if segdone(s)
        % Skip to next segment if this one has already been done
        % (i.e. this segment overlaps with another)
        continue;
    end

    % Extract lfiw x/y points for this coastal segment
    if strcmp(csegs(s).direction, 'right')
        csi = csegs(s).v0:csegs(s).v1;
    else
        csi = csegs(s).v1:-1:csegs(s).v0;
    end
    seg.x = sx(csi);
    seg.y = sy(csi);
    seg.w = lfiw(csi);

    % Remove any NaNs from segment coordinates
    % - not entirely sure how they would have go there
    %   but we can effectively interpolate over them 
    nani = isnan(seg.x) | isnan(seg.y);
    seg.x = seg.x(~nani);
    seg.y = seg.y(~nani);
    seg.w = seg.w(~nani);


    % Check to see if this segment overlaps with another
    if ~strcmp(csegs(s).overlap, "none") && segdone(s) == 0
        % If so, combine sx,sy coords 
        % and find the shortest path connecting them all
        
        o = find(strcmp({csegs(:).segname}, csegs(s).overlap));
        if strcmp(csegs(o).direction, 'right')
            osi = csegs(o).v0:csegs(o).v1;
        else
            osi = csegs(o).v1:-1:csegs(o).v0;
        end
        oseg.x = sx(osi);
        oseg.y = sy(osi);
        oseg.w = lfiw(osi);
        
        % Combine the overlapping SLIE segments
        seg = slie_combineoverlappingSLIEs(seg, oseg, cvecs(csi,:), ...
                                                      cvecs(osi,:));

        % Mark these segments as done
        segdone([s,o]) = 1;
    end

    % Check to see if this segment abuts a previous segment
    % - if so, this can create an unrealistic jump in the SLIE
    if ~strcmp(csegs(s).abut0, 'none')
        % Get abutting x,.y coords
        abutx = abut.(csegs(s).abut0).x;
        abuty = abut.(csegs(s).abut0).y;

        % Calculate distance from SLIE start point for this segment
        dx = abutx(end) - seg.x(1);
        dy = abuty(end) - seg.y(1);
        dd = sqrt(dx^2 + dy^2);
        if dd < 10
                % If the SLIEs are close, force them to abut
                seg.x = [abutx; seg.x];
                seg.y = [abuty; seg.y];
        end
    end

    % Check to see if this segment abuts next segment
    % - if so, assign abutting points
    if ~strcmp(csegs(s).abut1, 'none')
        % Get abutting x,.y coords
        abut.(csegs(s).abut1).x = [cvecs(csegs(s).v1,1); seg.x(end)];
        abut.(csegs(s).abut1).y = [cvecs(csegs(s).v1,2); seg.y(end)];
    end

    % Apply median filter to remove spikes if specified in function call
    if (filtwid > 0)
        seg.x = round(medianfilter(seg.x, filtwid));
        seg.y = round(medianfilter(seg.y, filtwid));
    end

    % Append closure points defined in coastal segments file
    seg.x = [seg.x; csegs(s).xclose; seg.x(1)];
    seg.y = [seg.y; csegs(s).yclose; seg.y(1)];

    % Pre-assign arrays to hold interpolated points
    % (overestimate required length)
    seg.xi = repmat(seg.x, 10, 1)*0;
    seg.yi = repmat(seg.y, 10, 1)*0;

    % Interpolate between all points to create an unbroken 
    % line of pixels outlining the landfast ice
    seg.xi(1) = seg.x(1);
    seg.yi(1) = seg.y(1);
    i = 1;
    for v=2:numel(seg.x)
    
        % Figure out length and direction of SLIE segment
        dx = seg.x(v) - seg.x(v-1);
        dy = seg.y(v) - seg.y(v-1);
        az = atan2(dy,dx);
        d = round(sqrt(dx^2 + dy^2));

        % Calculate coordinates of each pixel along SLIE segmemt
        r = 1:d;
        seg.xi(i+r) = seg.x(v-1) + round(r*cos(az));
        seg.yi(i+r) = seg.y(v-1) + round(r*sin(az));

        i = i + d;
    end


    % Trim arrays of interpolated coordinates
    seg.xi = seg.xi(1:i);
    seg.yi = seg.yi(1:i);

    % "Draw" outline onto SLIE image
    sii = sub2ind(sz, seg.yi,seg.xi);
    slieoutline(sii) = s;

            
    % Use flood fill to fill-in landfast ice area
    fillx0 = csegs(s).fillxy(1);
    filly0 = csegs(s).fillxy(2);
    sliefill = imfill(slieoutline==s, [filly0, fillx0]);

    % Assign this segment with a unique power-of-two value
    % This allows us to separate which pixel were filled by which segments.
    % Which is helpful for diagnostics
    slieimg = slieimg + (2^(s-1))*sliefill;
    
    
    % Diagnostics
    %imagesc(slieimg.*double(lmask==0) + double(lmask/2))
    %line(seg.xi, seg.yi, 'linestyle', 'none', ...
    %                 'marker', '.', ...
    %                 'color', 'blue');
    %line(fillx0, filly0, 'linestyle', 'none', ...
    %                    'marker', '.', ...
    %                    'color', 'red');
    
    notalot = 0;
end


% Binarize SLIE coverage and apply landmask, etc
slieimg = 255*(slieimg>0).*(lmask==0) + ...
          128*(lmask == 255) + ...
          64*(lmask == 128) + ...
          32*(lmask == 64);
      

% % DO NOT FILL-IN SHADOW ZONES. IT GETS AWFULLY COMPLICATED
% % ESPECIALLY IN KOTZEBUE SOUND WHEN TRYING TO COMBINE SLIES
% % FROM THE TWO DIFFERENT COAST VECTOR REGIONS

% % % Now fill-in coast vector shadow regions
% % for sh = 1:numel(shadregs)
% %     lfi_on_perim = slieimg(shadregs(sh).adjpix) == 255;
% %     if any(lfi_on_perim)
% %         slieimg(shadregs(sh).PixelIdxList) = 255;
% %     end
% % %     if numel(shad(sh).adjpix) > 50
% % %         x = slieimg;
% % %         x(shad(sh).adjpix) = 500;
% % %         imagesc(x);
% % %         notalot = 0;
% % %     end
% % 
% % end
      
      
      
slieimg = uint8(slieimg);

notalot = 0;


end