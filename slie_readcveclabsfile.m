function [cveclabs, labs_i, Nlabs]  = slie_readcveclabsfile(infile)
% Function to read in the labels for each coast vector
%
% Ignores any lines beginning with '#', '%', or ';'
%
% USAGE
%                  cveclabs = slie_readcveclabsfile(infile)
%  [cveclabs,labs_i, Nlabs] = slie_readcveclabsfile(infile)
%
%
% INPUTS
%    infile : string specifying full path and filename of a 
%             comma-separated ASCII file containing a list of 
%             coast vector indices and their associated text labels

%
% OUTPUT
%  cveclabs : Nlabs-by-1 cell array of coast vector labels
%    labs_i : Nlabs-by-1 array providing the coast vector indices 
%             corresponding to each label
%     Nlabs : A scalar integer specifying the number of labels
%
%
% Andy Mahoney October 2023
% ---------------------------------------------------------------------


% Open file for reading as text file
infid = fopen(infile, 'rt');


% Get file size in bytes in order to estimate # of lines
fseek(infid, 0, 'eof');
filesz = ftell(infid);

% Assume 10 chars per line
% (an underestimate, but this just means
%  we over-estimate the number of linesa bit)
charperline = 10;
Nl_guess = ceil(filesz/(charperline));
%disp(['Reading ',infile,': ',num2str(filesz,'%10.0f')]);
%disp(['(guesstimated ',num2str(Nl_guess),' lines)']);

% Pre-allocate arrays base on this guess
cveclabs =  cell(Nl_guess,1);
labs_i = zeros(Nl_guess, 1);

% Read in file data
fseek(infid,0,'bof'); % Go back to beginning of file

% Skip comment lines
Nlabs = 0;
while (feof(infid) == 0)
    line = fgetl(infid);

    % Skip any line starting with ;, %, or #
    if ~((line(1) == ';') || (line(1) == '#') || ...
         (line(1) == '%') || (numel(line) == 0))      
        Nlabs = Nlabs + 1;
        splitline = regexp(line, ' *, *', 'split'); % Ignore white spaces before and after delimiter
        
        labs_i(Nlabs) = str2double(splitline{1})+1; % Have to add 1 to convert fropm IDL!
        cveclabs{Nlabs} = splitline{2};
    end 
end


%Close file
fclose(infid);

% Trim output arrays to actual number of labels found
cveclabs = cveclabs(1:Nlabs);
labs_i = labs_i(1:Nlabs);



end

