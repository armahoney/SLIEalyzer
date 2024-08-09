function slie_calc_sfctemp_at_cps(varargin)

Nvarargs = numel(varargin);
a = 1;
while a < Nvarargs
    switch varargin{a}
        case 'dataset'
            dataset = varargin{a+1};
            a = a + 2;
        case 'cfgfile'
            cfgfile = varargin{a+1};
            a = a + 2;
        case 'yr0'
            yr0 = varargin{a+1};
            a = a + 2;
        case 'yr1'
            yr1 = varargin{a+1};
            a = a + 2;
        case 'ncfile'
            ncfile = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end    
    
    
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Identify the location where data are stored for specfied region
% Set default values for parameters not specified in function call
if ~exist('dataset', 'var')
    dataset = 'beau';
end
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end


% Read SLIE config file
cfg = slie_readconfigfile(sfolder, cfgfile);


% Assign default NetCDF file
if ~exist('ncfile', 'var')
    %ncfile = '~/DATA/NCEP/SfcAirTemp_ChukBeau_1948-2011.nc';
    %ncfile = '~/DATA/NCEP/NCEP_airtemp_AK_1948-2023.nc';
    ncfile = [sfolder filesep '../ConfigData/NCEP_airtemp_AK_1948-2023.nc'];
end


% Figure out date stuff
if ~exist('yr0', 'var')
    % Use config-specified seasons to determine start/end dates
    seasons = cfg.seasons;
    Nseasons = numel(seasons);
    syr0 = zeros(Nseasons,1);
    for s = 1:Nseasons
        syr0(s) = str2double(seasons{s}(1:4));
    end
    yr0 = min(syr0);
    yr1 = max(syr0)+1;
    date0 = datenum(yr0,7,1);
    date1 = datenum(yr1,8,31);
else
    % Determine seasons from user-specified dates
    date0 = datenum([yr0,7,1]);
    date1 = datenum([yr1,8,31]);
    syr0 = yr0:(yr1-1);
    Nseasons = numel(syr0);
    seasons = cell(Nseasons,1);
    for s=1:Nseasons
        syr0str = num2str(syr0(s), '%4d');
        syr1str = num2str(syr0(s)+1, '%4d');
        seasons{s} = [syr0str '-' syr1str(3:4)];
    end
end
Nd = 1 + date1-date0;
t = date0:date1;

% Read land mask to get GeoTIFF information
[~, R] = geotiffread(cfg.lmaskfile);
gtinfo = geotiffinfo(cfg.lmaskfile);
geokeyDT = gtinfo.GeoTIFFTags.GeoKeyDirectoryTag;
% Remove any "unknown" fields from geokey
if isfield(geokeyDT, 'Unknown')
    geokeyDT = rmfield(geokeyDT, 'Unknown');
end

% Convert coast vector pixel coords to lat/lons
cvecs = slie_readcvecsfile(cfg.cvecfile);
[cpgx, cpgy] = R.intrinsicToWorld(cvecs(:,1),cvecs(:,2));
[cplat, cplon] = projinv(gtinfo, cpgx, cpgy);
cplon = cplon + 360*(cplon < 0);
Ncp = numel(cplat);


% Read in NCEP data file (and convert to degC)
stemp = nc_varget(ncfile,'air') - 273.15;

% stemp_info = nc_getvarinfo(ncfile, 'air');
% stemp_units = nc_attget(ncfile, 'air', 'units');
% stemp_offset = nc_attget(ncfile, 'air', 'add_offset');
% stemp_sfactor = nc_attget(ncfile, 'air', 'scale_factor');



        
% Get lat / lon data (and make lon 0-360)
nlat=nc_varget(ncfile,'lat');
nlon=nc_varget(ncfile,'lon');
nlon = nlon + 360*(nlon<0);

nlon_g = repmat(nlon', numel(nlat), 1); 
nlat_g = repmat(nlat, 1, numel(nlon)); 

% Get time and convert to Julian day
% (NCEP time is in "hours since 1-1-1")
ntime = nc_varget(ncfile, 'time');
jtime = datenum(1800,1,1,0,0,0)+(ntime/24);

% Determine time range of data in file
ntime_range = nc_attget(ncfile, 'time', 'actual_range');
ndate0 = datenum(1800,1,1,0,0,0)+(ntime_range(1)/24);
ndate1 = datenum(1800,1,1,0,0,0)+(ntime_range(2)/24);

% % Make sure this is daily data
% ntime_dt = datenum(nc_attget(ncfile, 'time', 'delta_t'));
% if ntime_dt ~= 1
%     disp('Warning: NetCDF data are not daily values');
% end

%Initialize temperature array
cptemp = zeros(Ncp, Nd);


% Go through day-by-day and determine temperature at each coast point
for d=1:Nd
    % Determine day index for NCEP data
    jd = date0 + d -1;
    jti = floor(jtime) == jd;
    
    %disp(datestr(jd, 'yyyy/mm/dd'));
    
    % Extract average temparature for day
    stemp_d = squeeze(mean(stemp(jti,:,:), 1));    
    
    % Interpolate temperature at each coastpoint
    stemp_p = interp2(nlon_g, nlat_g, stemp_d, cplon, cplat);
    
    cptemp(:,d) = stemp_p;
end

disp('Found daily temperatures at coast points');



% ----------------------------------------------------------------
% Now calculate onset of freeze / thaw and accumated FDDs and TDDs

% Define thaw and freeze windows
% (used to narrow down search for freeze/thaw onsets)
thawwin0 = datenum(1,03,01)-datenum(0,12,31);
thawwin1 = datenum(1,07,31)-datenum(0,12,31);
frzwin0 = datenum(1,08,01)-datenum(0,12,31);
frzwin1 = datenum(1,12,01)-datenum(0,12,31);

fwdlook = 5; % Determines how many days ahead to look for consecutive
             % freezing / thawing temperatures
kpad = floor(fwdlook/2);
kern = [zeros(kpad,1); ones(fwdlook,1); zeros(kpad,1)];
             
            
% Figure out first, last and number of calendar years in time period
[yr0, mo0, da0] = datevec(date0); 
[yr1, mo1, da1] = datevec(date1);
Nyr = 1+yr1-yr0;

% Define arrays to hold freeze/thaw onsets
thaw_onset = zeros(Ncp, Nseasons);
frz_onset = zeros(Ncp, Nseasons);


for s=1:Nseasons
    % Determine calendar year for freeze onset
    fyr = syr0(s);
        
    % Define freeze onset window
    % and extract data for this period
    fw0 = datenum(fyr,1,1) + frzwin0 - 1;
    fw1 = datenum(fyr,1,1) + frzwin1 - 1;
    fwi = find((t >= fw0) & (t <= fw1));
    Nfw = numel(fwi);
    if Nfw > 0
        for p = 1:Ncp
            
            % Identify all negative temperatures 
            negtemps = 1*(cptemp(p,fwi) < 0);
            
            % Calculate forward-looking running sum
            negtemps_conv = conv(negtemps, kern');
            negtemps_fwdsum = negtemps_conv(kpad:end); % This will be shotr array
            
            % Determine first occerence of fwdlook consecutive
            % negative temperatures
            frz_i = find(negtemps_fwdsum >= fwdlook, 1, 'first');
            
            if any(frz_i)
                frz_onset(p,s) = fw0 + frz_i;
            else
                frz_onset(p,s) = NaN;
            end
        end
    else
        %frz_onset(:,s) = -NaN;
        notmuch = 0;
    end
    
    % Determine calendar year for thaw onset
    tyr = syr0(s)+1;
    
    
    % Define thaw onset window
    % and extract data for this period
    tw0 = datenum(tyr,1,1) + thawwin0 - 1;
    tw1 = datenum(tyr,1,1) + thawwin1 - 1;
    twi = find((t >= tw0) & (t <= tw1));
    Ntw = numel(twi);

    if Ntw > 0 
        for p = 1:Ncp
            
            % Identify all positive temperatures 
            postemps = 1*(cptemp(p,twi) > 0);
            
            % Calculate forward-looking running sum
            postemps_conv = conv(postemps, kern', 'valid');
            postemps_fwdsum = postemps_conv(kpad:end);
            
            % Determine first occerence of fwdlook consecutive
            % postive temperatures
            thaw_i = find(postemps_fwdsum >= fwdlook, 1, 'first');
            
            if any(thaw_i)
                thaw_onset(p,s) = tw0 + thaw_i;
            else
                thaw_onset(p,s) = NaN;
            end
        end
    else
        %thaw_onset(:,s) = -NaN;
        notmuch = 0;
    end
    
    
    

end



% Now calculate regional means & stdevs of these dates
[cvecregs, regnames, regbounds, Nregs] = read_cvecregsfile(cfg.cvecregionsfile);
% Define output arrays
thaw_onset_regmean = zeros(Nregs, Nyr+1);
thaw_onset_regstdev = zeros(Nregs, Nyr+1);
frz_onset_regmean = zeros(Nregs, Nyr+1);
frz_onset_regstdev = zeros(Nregs, Nyr+1);
for r = 1:Nregs
    ri = cvecregs == r;
    ri = ri';
    
    % Calculate mean and stdev for each year
    for s = 1: Nseasons
        thaw_onset_rs = thaw_onset(ri,s);
        fini = isfinite(thaw_onset_rs);
        Nfin = sum(fini);
        
        if Nfin > 0                
            thaw_onset_regmean(r,s) = mean(thaw_onset_rs(fini));
        end
        
        if Nfin > 2
            thaw_onset_regstdev(r,s) = std(thaw_onset_rs(fini));
        end

        frz_onset_rs = frz_onset(ri,s);
        fini = isfinite(frz_onset_rs);
        Nfin = sum(fini);
        
        if Nfin > 0                
            frz_onset_regmean(r,s) = mean(frz_onset_rs(fini));
        end
        
        if Nfin > 2
            frz_onset_regstdev(r,s) = std(frz_onset_rs(fini));
        end
    end
    
    % Calculate mean and stdev for all seasons
    thaw_onset_r = dayofyear(thaw_onset(ri,:));
    fini = isfinite(thaw_onset_r);
    Nfin = sum(fini(:));
    
    if Nfin > 0
        thaw_onset_regmean(r,Nyr+1) = mean(thaw_onset_r(fini));
    end
    
    if Nfin > 2
        thaw_onset_regstdev(r,Nyr+1) = std(thaw_onset_r(fini));
    end
    
    frz_onset_r = dayofyear(frz_onset(ri,:));
    fini = isfinite(frz_onset_r);
    Nfin = sum(fini(:));
    
    if Nfin > 0
        frz_onset_regmean(r,Nseasons+1) = mean(frz_onset_r(fini));
    end
    
    if Nfin > 2
        frz_onset_regstdev(r,Nseasons+1) = std(frz_onset_r(fini));
    end
    
end


% Round regional means to nearest day
thaw_onset_regmean = round(thaw_onset_regmean);
frz_onset_regmean = round(frz_onset_regmean);


% -------------------------------------------------------------------------
% Write these data to ASCII files
% *** AND CONVERT JULIAN DAYS TO DAYS OF YEAR 
frzthaw_folder = [sfolder filesep 'AllSeasonsAnalysis2023' filesep 'FrzThaw_analysis'];
if ~exist(frzthaw_folder, 'file'), mkdir(frzthaw_folder); end

[labs, labs_i, Nlabs] = read_cveclabsfile(cfg.cveclabelsfile);
cveclabs = cell(Ncp,1);
cveclabs(labs_i) = labs;



% Thaw onset values
outfile = [frzthaw_folder filesep 'cp_thaw_onset.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, '%s', 'CP label, CP number');
fprintf(outfid, ',%s', seasons{:});
fprintf(outfid, '%s\n', '');
for p=1:Ncp
    fprintf(outfid, '%s\n', [cveclabs{p} ',' num2str(p, '%4d') ...
                             num2str(dayofyear(thaw_onset(p,:)), ',%3d')]);
end
fclose(outfid);
disp(['Written: ', outfile]);

% Regional mean thaw onset values
outfile = [frzthaw_folder filesep 'cp_thaw_onset_regmean.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, ',%s',  seasons{:});
fprintf(outfid, ',%s\n', 'All seasons');
for r=1:Nregs
    fprintf(outfid, '%s\n', [regnames{r} ...
                             num2str(dayofyear(thaw_onset_regmean(r,:)), ',%3d')]);
end
fclose(outfid);
disp(['Written: ', outfile]);

% Regional stdev thaw onset values
outfile = [frzthaw_folder filesep 'cp_thaw_onset_regstdev.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, ',%s',  seasons{:});
fprintf(outfid, ',%s\n', 'All seasons');
for r=1:Nregs
    fprintf(outfid, '%s\n', [regnames{r} ...
                             num2str(thaw_onset_regstdev(r,:), ',%5.1f')]);
end
fclose(outfid);
disp(['Written: ', outfile]);

% Freeze onset values
outfile = [frzthaw_folder filesep 'cp_frz_onset.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, '%s', 'CP label, CP number');
fprintf(outfid, ',%s', seasons{:});
fprintf(outfid, '%s\n', '');
for p=1:Ncp
    fprintf(outfid, '%s\n', [cveclabs{p} ',' num2str(p, '%4d') ...
                             num2str(dayofyear(frz_onset(p,:)), ',%3d')]);
end
fclose(outfid);
disp(['Written: ', outfile]);


% Regional mean freeze onset values
outfile = [frzthaw_folder filesep 'cp_frz_onset_regmean.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, ',%s',  seasons{:});
fprintf(outfid, ',%s\n', 'All seasons');
for r=1:Nregs
    fprintf(outfid, '%s\n', [regnames{r} ...
                             num2str(dayofyear(frz_onset_regmean(r,:)), ',%3d')]);
end
fclose(outfid);
disp(['Written: ', outfile]);

% Regional stdev freeze onset values
outfile = [frzthaw_folder filesep 'cp_frz_onset_regstdev.csv'];
outfid = fopen(outfile, 'wt');
fprintf(outfid, ',%s',  seasons{:});
fprintf(outfid, ',%s\n', 'All seasons');
for r=1:Nregs
    fprintf(outfid, '%s\n', [regnames{r} ...
                             num2str(frz_onset_regstdev(r,:), ',%5.1f')]);
end
fclose(outfid);
disp(['Written: ', outfile]);









% -----------------------------------------------------------------
% Now determine the number of accumulated TDDs and FDDs in each season
tdd = zeros(Ncp, Nd)+NaN;
fdd = zeros(Ncp, Nd)+NaN;

for s = 1:Nseasons
    
    yr = syr0(s);
    
    % Determine relative day indices of 
    % first and last days of this season
    sd0 = 1+ datenum(yr,8,1) - date0;
    sd1 = 1+ datenum(yr+1,7,31) - date0;
    
    disp([yr, sd0, sd1, Nd]);
    
    for p = 1:Ncp
        
        
        % Look up thaw and freeze dates for this season
        % and translate to indices of data time period 
        fd0 = 1+ frz_onset(p,s) - date0;
        td0 = 1+ thaw_onset(p,s) - date0;
        
        
        % If this seasons is consecutive with last
        % continue computing TDDs up unil freeze onset
        % otherwise set to NaN
        if s > 1
            if syr0(s) == syr0(s-1)+1
                consec = 1;
            else
                consec = 0;
            end
        else
            consec = 0;
        end
        
        if consec
            tdd_last = tdd(p,sd0-1);
        else
            tdd_last = NaN;
        end
        thaw_temps = cptemp(p,sd0:fd0);
        thaw_temps = thaw_temps.*(thaw_temps > 0);
        tdd(p,sd0:fd0) = cumsum(thaw_temps)+tdd_last;
        
      
        % Now calculate cumulative FDDs from freeze onset to thaw onset
        frz_temps = cptemp(p,fd0:td0);
        frz_temps = frz_temps.*(frz_temps < 0);
        fdd(p,fd0:td0) = cumsum(frz_temps);
        
            
        % and calculate cumulative TDDs from thaw onset to end of season        
        thaw_temps = cptemp(p,td0:sd1);
        thaw_temps = thaw_temps.*(thaw_temps > 0);
        tdd(p,td0:sd1) = cumsum(thaw_temps);
    end    
end

% Write these data to binary files
outfile = [frzthaw_folder filesep 'cp_TDDs.dat'];
outfid = fopen(outfile, 'w');
fwrite(outfid, [date0, date1, Ncp], 'double');
fwrite(outfid, tdd, 'double');
fclose(outfid);
disp(['Written: ', outfile]);


outfile = [frzthaw_folder filesep 'cp_FDDs.dat'];
outfid = fopen(outfile, 'w');
fwrite(outfid, [date0, date1, Ncp], 'double');
fwrite(outfid, fdd, 'double');
fclose(outfid);
disp(['Written: ', outfile]);
    

% Save all output data in a MatLab file
outname = [frzthaw_folder filesep 'cp_freezethaw.mat'];
save(outname, 'date0', 'date1', 'Ncp', 'yr0', 'yr1', ...
              'cptemp', 'thaw_onset', 'frz_onset', ...
              'thaw_onset_regmean', 'frz_onset_regmean', ...
              'thaw_onset_regstdev', 'frz_onset_regstdev', ...
              'tdd', 'fdd');
disp(['Written: ', outname]);
          

end
        
 


