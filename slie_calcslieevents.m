function ev = slie_calcslieevents(lfiw, wdep, sliedates, varargin)
% A function to determine key events within a landfast ice season 
% along each coast vector
%
% USAGE:
%      ev = calc_slieevents(lfiw, wdep, sliedates)
%      ev = calc_slieevents(__, 'min_ice', min_ice)
%      ev = calc_slieevents(__, 'stab_dep', stab_dep) 
%      ev = calc_slieevents(__, 'specwid', specwid)
%
% INPUTS:
%      lfiw : Nvec-by-Nslie array of landfast ice width measurements
%             for a single season
%      wdep : Nvec-by-Nslie array of water depth at SLIE measurements
%             for a single season
% sliedates : Nvec-by-Nslie array of dates for each lfiw and wdep
%             data point
%   min_ice : scalar specifiying minimum width to qualify for first/last
%             ice
%             default: 0.5 (km)
%  stab_dep : scalar specifiying water depth at SLIE to qualify for 
%             onset of landfsat ice stability
%             default: 15 (m)
%   specwid : scalar specifiying a width for determine date range when
%             landfast ice continuously exceeded this given width
%             default: 10 (km)
%
% OUTPUTS:
%    ev : a structure containing the following fields each with a 
%         1-by-Nvecs vector of dates
%       - firstice    : first appearance of landfast ice that grows
%                       to exceed min_wid in width
%       - stableice   : the date on which the ice becomes stable
%                       (water depth >= stab_dep)
%       - maxice      : the date on which the maximum extent occurs
%       - unstableice : the date on which the ice stability ceases
%                       (water depth < stab_dep)
%       - breakup     : onset of breakup (most abrupt loss of width during
%                       period of continuous retreat)
%       - icefree     : the start of the ice-free season
%                       (width drops below min_wid and stays there)
%       - firstwid    : first date landfast ice reliably wider than 
%                       specified width
%       - lastwid     : last date landfast ice reiably wider than 
%                       specified width
%
% Andy Mahoney - November 2023
%              - January 2023 : amended definitions of first and last ice
%                               to return NaNs when the beginning or end
%                               of season is not captured
% ------------------------------------------------------------------------

% Handle any variable arguments specified with function call
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'min_ice'
            min_ice = varargin{a+1};
            a = a + 2;
        case 'stab_dep'
            stab_dep = varargin{a+1};
            a = a + 2;
        case 'specwid'
            specwid = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end  


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set default values for parameters not specified in function call
% Default value of specwid for events 6 and 7
if ~exist('min_ice', 'var')
    min_ice = 0.5;
end

if ~exist('stab_dep', 'var')
    stab_dep = 15;
end

if ~exist('specwid', 'var')
    specwid = 10;
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Figure out size of data arrays
datasz = size(lfiw);
Nslies = datasz(2);
Nvecs = datasz(1);

% Preassign output arrays and data structure
ev.firstice = zeros(Nvecs,1);
ev.stableice = zeros(Nvecs,1);
ev.maxice = zeros(Nvecs,1);
ev.unstableice = zeros(Nvecs,1);
ev.breakup = zeros(Nvecs,1);
ev.icefree = zeros(Nvecs,1);
ev.firstwid = zeros(Nvecs,1);
ev.lastwid = zeros(Nvecs,1);

% Assign numbers to each event
evnum = {'1','2','3','4','5','6','7','8'};

% Convert sliedates to days since Jan 1 of year 0 of the season
[yr, mo, da] = datevec(sliedates);
year0 = min(yr(:));
sliedoy = sliedates - datenum(year0-1,12,31);

if year0 == 2007 
 notmuch = 0;
end

for v=1:Nvecs
    
    w = lfiw(v,:);
    d = wdep(v,:);
    
    % Exclude NaNs from record
    % ** Don't forget to index using fini later!!
    fini = find(~isnan(w) & ~isnan(d));
    Nfin = numel(fini);
    w = w(fini);
    d = d(fini);
    
    % If there's no ice on the coast, everything is NaN
    if ~any(w > 0)
        ev.firstice(v) = NaN;
        ev.icefree(v) = NaN;
        ev.stableice(v) = NaN;
        ev.unstableice(v) = NaN;
        ev.maxice(v) = NaN;
        ev.breakup(v) = NaN;
        ev.firstwid(v) = NaN;
        ev.lastwid(v) = NaN;
        continue;
    end

        
    % FIRST ICE ----------------------------------------------
    % - Find the first time landfast ice width > min_ice  
    % - Find the last time width was zero before this time
    % - Find the first non-zero width between these two times
    firstGEmin = find(w >= min_ice, 1, 'first');
    lastzero = find(w(1:firstGEmin) == 0, 1, 'last');
    firsti = (lastzero-1) + find(w(lastzero:firstGEmin) > 0, 1, 'first');
    if any(firsti)
        ev.firstice(v) = sliedoy(v,fini(firsti));
    else
        % Set date to NaN if there are no non-nan width observations 
        % before first observation > min_wid
        ev.firstice(v) = NaN;
    end

    % LAST ICE ----------------------------------------------
    % - Find the last time landfast ice width > min_ice  
    % - Find the first time width was zero after this time
    % - Find the last non-zero width between these two times
    lastGEmin = find(w >= min_ice, 1, 'last');
    firstzero = (lastGEmin-1) + find(w(lastGEmin:end) == 0, 1, 'first');
    lasti = (lastGEmin -1) + find(w(lastGEmin:firstzero) > 0, 1, 'last');
    if any(lasti)
        ev.icefree(v) = sliedoy(v,fini(lasti));
    else
        % Set date to NaN if there are no non-nan width observations 
        % after last observation > min_wid
        ev.icefree(v) = NaN;
    end

    % STABLE / UNSTABLE ICE ----------------------------------------------
    % Determine when SLIE was in water >= 15m deep
    % stable ice starts and ends with the longest conitinuous period
    deep = (d >= stab_dep);
    deepdiff = diff(deep);
    stab0 = find([deep(1),deepdiff] == 1);       % Find where water deep increases to >= 15m
    stab1 = find([deepdiff, -deep(end)] == -1);  % Find where water deep decreases to <= 15m
    Nstab0 = numel(stab0);
    Nstab1 = numel(stab1);
    if (Nstab0 > 0) && (Nstab0 == Nstab1)
        [stab_dur, maxi] = max(stab1-stab0);
        stabi = stab0(maxi);
        unstabi = stab1(maxi);
        ev.stableice(v) = sliedoy(v,fini(stabi));
        ev.unstableice(v) = sliedoy(v,fini(unstabi));
    else
        if Nstab0 == 0
            stabi = NaN;
            unstabi = NaN;
            ev.stableice(v) = NaN;      % For help in distinguishing why
            ev.unstableice(v) = NaN;    % a stablization date wasn't found
        else                    % Use + and - NaN values
            stabi = NaN;
            unstabi = NaN;
            ev.stableice(v) = -NaN;
            ev.unstableice(v) = -NaN;
        end
    end
        
    % MAX WIDTH -----------------------------------------------------------
    % Determine when maximum lfiw occurred
    [~, maxi] = max(w);
    ev.maxice(v) = sliedoy(v, fini(maxi(1)));
        

    % BREAKUP -------------------------------------------------------------
    % Determine when breakup occurred
    % Define this by the time of steepest negative gradient
    % that is followed by exclusively negative gradients or less than minimum ice
    if ~isnan(ev.icefree(v))
        dw = diff(w(1:lasti+1));
        Ndw = numel(dw);
        tailneg = zeros(Ndw,1);
        for t=1:Ndw
            dwtail = dw(t:Ndw);                              % Tailneg is true if whole tail from t onwards is negative
            tailneg(t) = sum(dwtail < 0.5) == numel(dwtail); % (or less than a 0.5km increase)     
        end
        allnegtail = find(tailneg);
        if numel(allnegtail > 0)
            %[mostneg, mostnegi] = min(dw(allnegtail));
            %bkupi = allnegtail(mostnegi(1));
            taildrop = -dw(allnegtail) ./ w(allnegtail);
            firstbigdrop = find(taildrop > 0.25, 1, 'first');
            if ~any(firstbigdrop)
                firstbigdrop = find(taildrop > 0, 1, 'first');
            end
            bkupi = allnegtail(firstbigdrop);
            ev.breakup(v) = sliedoy(v, fini(bkupi));
        else
            bkupi = NaN;
            ev.breakup(v) = NaN;
        end
    else
        bkupi = NaN;
        ev.breakup(v) = NaN;
    end
        
        
    % SPECIFIED WIDTH DATES -----------------------------------------------
    % Determine when SLIE was at least specwid km from shore
    % and identify beginning and end of 
    % longest conitinuous period when this was so
    gespecwid = w >= specwid;
    gechange = diff(gespecwid);
    upi = find([gespecwid(1),gechange]==1);
    downi = find([gechange,-gespecwid(end)]==-1);
    Nup = numel(upi);
    Ndown = numel(downi);
    if (Nup > 0) && (Nup == Ndown)
        updoy = sliedoy(v,fini(upi));
        downdoy = sliedoy(v,fini(downi));
        [specwid_dur, maxi] = max(updoy-downdoy);
        ev.firstwid(v) = updoy(maxi);
        ev.lastwid(v) = downdoy(maxi);
    else
        if Nup == 0
            upi = NaN;
            downi = NaN;
            ev.firstwid(v) = NaN;      % For help in distinguishing why
            ev.lastwid(v) = NaN;       % a stablization date wasn't found
        else                           % Use + and - NaN values
            upi = NaN;
            downi = NaN;
            ev.firstwid(v) = -NaN;
            ev.lastwid(v) = -NaN;
        end
    end

        
        
    % % Some diagnostics        
    % plot(sliedates(v,:), w);
    % evi = [firsti,stabi,maxi,unstabi,bkupi,lasti];
    % fini = isfinite(evi);
    % evifin = evi(fini);
    % line(sliedates(v,evifin), w(evifin), 'marker', 'o', ...
    %                                      'color', 'red', ...
    %                                      'linestyle', 'none'); 
    % text(sliedates(v,evifin), w(evifin), evnum(fini), ...
    %                                    'verticalal', 'bottom', ...
    %                                    'horizontalal', 'center');     
    
    notmuch = 0;
        
  
end
 
end

  
  
