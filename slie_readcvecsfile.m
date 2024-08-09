function [cvecs, v] = slie_readcvecsfile(infile)
% Function to read comma-separated ASCII file containing
% coordinate pairs that define a series of coast vectors
%
% Ignores any lines beginning with '#'
%
% USAGE
%        cvecs = slie_readcvecsfile(infile)
%   [cvecs, v] = slie_readcvecsfile(infile)
%
%
% INPUTS
%    infile : string specifying full path and filename of a 
%             comma-separated ASCII file containing a headerless 
%             4-column table of integers specifying the coast vector 
%             coordinates in pixel space:
%             - column 1 : x0
%             - column 2 : x1
%             - column 3 : y0
%             - column 4 : y1

%
% OUTPUT
%     cvecs : Nvec-by-4 array of coast vector coordinates
%         v : Nvec-by-1 list of indices for each coast vector
%
% Andy Mahoney October 2023
% ---------------------------------------------------------------------


% Open file for reading as text file
infid = fopen(infile, 'rt');

% Count number of lines in file
Nl = 0;
while ~feof(infid)
    fgetl(infid);
    Nl = Nl + 1;
end
fseek(infid, 0, 'bof');

% Pre-assign array to hold coast vector coordinates
cvecs = zeros(Nl,4);

% Read coordinates from file
v = 0;
while ~feof(infid)
    vecline = fgetl(infid);

    % Some new f'd up problem with non-printing space-like
    % character appearing at beginning of text files
    if double(vecline(1)) == 65279
        vecline = vecline(2:end);
    end
    
    % Only read lines that don't begin with a #
    if ~strcmp(vecline(1), '#')
        v = v + 1;
        xyxy = textscan(vecline, '%d %d %d %d', 'delimiter', ',');
        cvecs(v,:) = [xyxy{:}]+1; % Add one so that pixel numbers begin at 1
    end
end

% Trim cvec array to actual number of vectors
cvecs = cvecs(1:v,:);


end
