function [sdates, pref, suff] = slie_name2date(slienames, pos)
% Function to convert SLIE filenames to a date number
% Handles both the "standard" file names of the form:
%     pyyyy_ddd-ddd_slie.tif
% As well as daily file names of the form:
%     pyyyymmdd_dailyslie.tif
% 
%
% USAGE 
%    sdates = slie_name2date(slienames)
%    sdates = slie_name2date(slienames, pos)
%    [sdates, pref] = slie_name2date(__)
%    [sdates, pref, suff] = slie_name2date(__)
%
% INPUTS
%    slienames : a string or cell array of strings containing one more
%                SLIE image filenames of the form p_yyyy_ddd-ddd_slie.tif
%          pos : a string specififying which date to output in the case
%                of non-daily SLIE filenames
%                - 'first': output first date
%                - 'mid': output mid-date between first and last
%                - 'last': output last date
%                Default: 'mid'
%
% OUTPUT
%       sdates : a Nslie-by-1 vector of date numbers
%         pref : a Nslie-by-1 cell array of filename prefices
%         suff : a Nslie-by-1 cell array of filename suffices
%
% Andy Mahoney - October 2023
%
% ------------------------------------------------------------------------


if ~exist('pos', 'var')
    pos = 'mid';
end

% Convert input sliename into a cell array if it isn't one already
if ~iscell(slienames)
    slienames = {slienames};
end

% Count number of SLIEs in input list
Nslie = numel(slienames);

% Check first file name to see if we're dealing with "standard" format
% SLIE names or daily filenames
dailycheck = ~isempty(regexp(slienames{1}, 'daily', 'once'));

if ~dailycheck
    % Initialize arrays to hold date values extracted from SLIE names
    yyyy = zeros(Nslie,1);
    ddd0 = zeros(Nslie,1);
    ddd1 = zeros(Nslie,1);
    pref = cell(Nslie,1);
    suff = cell(Nslie,1);
    
    % Define regular expression for parsing filenames
    sdaterexp = ['(?<prefix>[\s\w]+)' ...
                 '(?<yr>\d\d\d\d)[_-]' ...
                 '(?<doy0>\d\d\d)[_-]' ...
                 '(?<doy1>\d\d\d)[_-]' ...
                 '(?<suffix>.*)\.tif'];
    
    
    for s = 1:Nslie
        % Extract parts of filename using regular expression
        sdate = regexp(slienames{s}, sdaterexp, 'names');

        % Keep year and day info
        yyyy(s) = str2double(sdate.yr);
        ddd0(s) = str2double(sdate.doy0);
        ddd1(s) = str2double(sdate.doy1);

        % Keep prefix and suffix info
        pref{s} = sdate.prefix;
        suff{s} = sdate.suffix;

    end

    % Identify leap years
    leap = 1*(mod(yyyy, 4)==0);
    
    % Handle year wrap-arounds accounting for leap years
    ddd1x = ddd1 + (365+leap).*(ddd1 < ddd0);
    
    % Compute date number for Jan 1 of the year of the SLIE
    sdates = datenum(yyyy, ones(Nslie,1), ones(Nslie,1)) -1;
    
    
    switch pos
        case 'mid'
            sdates = round(sdates + (ddd0 + ddd1x)/2);
        case 'first'
            sdates = sdates + ddd0;
        case 'last'
            sdates = sdates + ddd1x;
    end
else
    % Initialize arrays to hold date values extracted from SLIE names
    yyyy = zeros(Nslie,1);
    mm = zeros(Nslie,1);
    dd = zeros(Nslie,1);
    pref = cell(Nslie,1);
    suff = cell(Nslie,1);

    % Define regular expression for parsing filenames
    sdaterexp = ['(?<prefix>[\s\w]+)' ...
                 '(?<yr>[1-2][0-9][0-9][0-9])' ...
                 '(?<mm>[0-1][0-9])' ...
                 '(?<dd>[0-3][0-9])[_-]' ...
                 '(?<suffix>.*)\.tif'];
        
    for s = 1:Nslie
        % Extract parts of filename using regular expression
        sdate = regexp(slienames{s}, sdaterexp, 'names');
        
        % Keep year and day info
        yyyy(s) = str2double(sdate.yr);
        mm(s) = str2double(sdate.mm);
        dd(s) = str2double(sdate.dd);

        % Keep prefix and suffix info
        pref{s} = sdate.prefix;
        suff{s} = sdate.suffix;
    end

    % Compute date numbers
    sdates = datenum([yyyy, mm, dd]);
end

end
