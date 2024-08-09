function mp = slie_findmissingpolygondata(lfiw, cvecs, pixsz)
% Identify cases of "missing polygons" that have been identified
% in ice chart data. These most likely originate from mis-coding of 
% polygons that ought to be identified as landfast ice.
% They are characterized by brief drops to zero lfiw that:
% - last for 3 or fewer consectutive datapoints
% - have similiar lfiw values before and after the drop
% - affect at least 25 consecutive coast vectors (~ 5km)
%
% USAGE
%     mp = slie_findmissingpolygondata(lfiw, cvecs, pixsz)
%
% INPUTS
%    lfiw : Nvec-by-Nslie array of landfast ice width values
%   cvecs : Nvec-by-4 array of coast vector coordinates
%   pixsz : pixel size
%
% OUTPUT
%      mp : Nvec-by-Nslie array of 0s and 1s indicating whether the
%           corresponding lfiw data point is the result of a missing
%           polygon (1 = missing)
%
%
% Andy Mahoney - October 2023


% -----------------------------------------------------------------------
% Get dimensions from input data
Nv = size(lfiw,1);
Ns = size(lfiw,2);

% Pre-assign array to hold missing polygon
mp = zeros(Nv,Ns);


minwid = 1;      % Minimum width to qualify for a missing polygon
minpolylen = 25; % Minimum missing polygon length in km

% Analyze timeseries of landfast ice width for each coast vector
% - identify any brief dropouts that meet the criteria above
for v=1:Nv
    w = lfiw(v,:);

    % For these purposes, any datapoints already identified as NaN values 
    % will be set to zero
    % - these represent regions missing SAR mosaic coverage
    w(isnan(w)) = 0;

    % Find datapoints greater than zero
    gtz  = w > 0;

    % Skip to next of there's no landfast at all
    if ~any(gtz)
        continue;
    end

    % Find instances where lfiw drops to zero (from > 0)
    % and where lfiw rises back up again (to > 0)
    alldrop = find(diff(gtz) == -1) + 1;
    allrise = find(diff(gtz) == 1);

    % Skip to next of there's no appearance / re-apperance of landfast
    % - i.e. landfast is present in first SLIE and does not subsequently
    %   reappear after a disappearance
    if numel(allrise) == 0
        continue;
    end

    % Skip to next of there's no disappearance of landfast
    % - i.e. landfast never goes away completely
    if numel(alldrop) == 0
        continue;
    end

    
    % We're only interested in rises that follow drop
    % (not the first time landfast ice appears at the coast)
    if (allrise(1) < alldrop(1))
        allrise = allrise(2:end);
    end

    % And we're only interested in drops that are followed by a rise
    if numel(allrise) > 0
        if alldrop(end) > allrise(end)
            alldrop = alldrop(1:end-1);
        end
    end

    % Calculate lengths of dropouts
    zlen = 1 + allrise - alldrop;

    % Exclude dropouts that last longer than 3 SLIEs
    shortdrop = alldrop(zlen <= 3);
    shortrise = allrise(zlen <= 3);

    % Calculate change either side of dropout
    beforedrop = w(shortdrop-1);
    afterrise = w(shortrise+1);
    change = abs(beforedrop-afterrise);
    
    % Determine a cut-off amount out change based on 
    % how much landfast ice is present before and after
    maxchange = 0.5*max([beforedrop,afterrise]);

    % Find those dropouts that:
    % 1. show little change
    % 2. have > minwid amount of landfast ice before and after
    % - each one of these represents one of more ice charts
    %   with a missing polyon that is intersected by this coast vector
    littlechange = find((change < maxchange) & ...
                         (beforedrop > minwid) & ...
                         (afterrise > minwid));

    % Assign a value of 1 to each data point in each drop out
    m = 0;
    for c=1:numel(littlechange)
        cc = littlechange(c);

        % Missing polygon data points in this drop out
        mpi = shortdrop(cc):shortrise(cc);

        mp(v,mpi) = 1;
        m = m + numel(mpi);
        notalot = 0;
    end
    notalot = 0;

    if any(isnan(mp(v,:)))
        notalot = 0;
    end

end


% Now exclude any dropouts that affect only a short length of coast
% - fewer than minpolylen consecutive coast vectors
for s=1:Ns
    if any(mp(:,s))
        mp0 = find([mp(1,s); diff(mp(:,s))] == 1);
        mp1 = find([diff(mp(:,s)); -mp(end,s)] == -1);
        Nmpseg = numel(mp0);
        if numel(mp1) ~= Nmpseg
            somethingsawry = 1;
        end
    
        % Calculate coastal distance 
        % between start and end of each segment
        dx = cvecs(mp1,1) - cvecs(mp0,1);
        dy = cvecs(mp1,2) - cvecs(mp0,2);
        mpseglen = sqrt(dx.^2 + dy.^2)*pixsz;
        for ss = 1:Nmpseg
            if mpseglen(ss) < minpolylen
                mp(mp0(ss):mp1(ss),s) = 0;
            end
        end
    end
end



end