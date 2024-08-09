function bocat = slie_readbreakoutcatalog(infile)
% Simple function to read ASCII data file created by slie_find_breakouts.m
%
% USAGE
% bocat = read_breakoutcatalog(infile)
%
% INPUT
%     infile : string specifying full path and filename for ASCII file
%              containing output from slie_find_breakouts.m
%
% OUTPUT
%      bocat : a structure containing fields named according to each
%              column of data in the ASCII file and containing a Nb-by-1
%              vector of either numeric or text data
%
% NOTES: 
%    The ASCII file created by slie_find_breakouts.m contains 9 columns:
%    brkoutid : integer specifying unique index number of each breakout
%      season : string specifying season in which each breakout occurred
%        slie : the filename of the SLIE image in which the breakout was
%               detected
%    sliedate : string specifying data associated with SLIE image file
%        vec0 : index of first coastvector bounding breakout
%        vec1 : index of last coastvector bounding breakout
%  meanwid_km : mean width of landfast ice lost in breakout (in km)
% coastlen_km : length of coast between vec0 and vec 1 (in km)
%    area_km2 : area of landfast ice lost in breakout (in km^2) 
%    latitude : latitude of centroid of breakout region
%   longitude : longitude of centroid of breakout region
%
% Andy Mahoney - October 2023
%
% -----------------------------------------------------------------------

% Open input for reading as text
infid = fopen(infile, 'rt');

% Read first line and extract column headers
headline = fgetl(infid);
colheads = lower(regexp(headline, ',', 'split'));

% Specify format of text data and read remainder of file
infstr = '%f %s %s %s %f %f %f %f %f %f %f';
indata = textscan(infid, infstr, 'delimiter', ',');

% Assign data read from file to named fields in output structure
for h=1:numel(colheads)
    bocat.(colheads{h}) = indata{h};
end

% Close input file
fclose(infid);

end
