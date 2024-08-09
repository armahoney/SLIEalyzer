function [evdates, Nvecs, Nseasons, veclabels] = slie_readeventdates(infile)
% Function to read in ASCII files generated by sliealyzer containing
% SLIE events for different seasons as functions of coastal location
% 
%
% USAGE
%            evdates = slie_readeventdates(infile)
%        [__, Nvecs] = slie_readeventdates(infile)
%     [__, Nseasons] = slie_readeventdates(infile)
%    [__, veclabels] = slie_readeventdates(infile)
%
% INPUTS
%    infile : string specifying full path and filename for ASCII file
%             containing event data 
%
% OUTPUT
%     evdates : Nvecs-by-Nseason array of datenumbers specifying the
%               occurence date for particular event
%       Nvecs : The number of coast vectors for which there are data
%               (i.e. number of rows in data file)
%    Nseasons : The number of seasons for which there are data
%               (i.e. number of columns in data file)
%   veclabels : Any coast vector labels read along the way
%
%
% Andy Mahoney - October 2023
%
%
% -----------------------------------------------------------------------


% Open file for reading as text
infid = fopen(infile, 'rt');


% Read in header to determine # of seasons
headline = fgetl(infid);
splitline = regexp(headline, ',', 'split');
Nseasons = numel(splitline)-2;


% Get file size in bytes in order to estimate # of lines
headlen = ftell(infid); % Note end of header position 
fseek(infid, 0, 'eof');
filesz = ftell(infid);
fseek(infid,headlen,'bof'); % Go back to end of header

% Assume 12 chars per season line
% (an underestimate, but this just means
%  we over-estimate the number of lines)
charperline = 12*Nseasons;
Nl_guess = ceil(filesz/(charperline));
%disp(['Reading ',infile,': ',num2str(filesz,'%10.0f')]);
%disp(['(guesstimated ',num2str(Nl_guess),' lines)']);


% Read in header to determine # of seasons
fseek(infid,0,'bof'); % Go back to beginning of file
headline = fgetl(infid);
splitline = regexp(headline, ',', 'split');
Nseasons = numel(splitline)-2;

% Initialize evdates array with guess
evdates = zeros(Nl_guess, Nseasons);
veclabels = cell(Nl_guess, 1);

l = 0;
while (feof(infid) == 0) 
    line = fgetl(infid);
    l = l + 1;
    
    splitline = regexp(line, ',', 'split');
        
    evdates(l,:) = str2double(splitline(3:end));
    veclabels{l} = splitline(1);
    
end

% Trim data arrays to actual number of coast vectors read
Nvecs = l;
evdates = evdates(1:Nvecs, :);
veclabels = veclabels(1:Nvecs);

% Close file
fclose(infid);


end

