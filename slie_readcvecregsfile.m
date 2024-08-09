function [cvecregs, regnames, regbounds, Nregs] = slie_readcvecregsfile(infile)
% Function to read in the start and end vectors and names 
% for different SLIE regions

% Ignores any lines beginning with '#', '%', or ';'
%
% USAGE
%  cvecregs = slie_readcvecregsfile(infile)
%  [__, regnames] = slie_readcvecregsfile(infile)
%  [__, regbounds] = slie_readcvecregsfile(infile)
%  [__, Nregs] = slie_readcvecregsfile(infile)
%
%
% INPUTS
%    infile : string specifying full path and filename of a 
%             comma-separated ASCII file containing a table of values
%             specifying the start vector, end vector, and name
%             of each subregion
%
%
% OUTPUT
%  cvecregs : Nvecs-by-1 array of integer specifying the sub-region
%             number that each coast vector lies within
%  regnames : Nregs-by-1 cell array of sub-region names
% regbounds : Nregs-by-2 array of integers specifying the starting 
%             and ending coast vector indices for each sub-region
%     Nregs : A scalar integer specifying the number of sub-regions
%
%
% Andy Mahoney October 2023
% ---------------------------------------------------------------------


% Open the vector regions ASCII file
infid = fopen(infile, 'rt');

% Get file size in bytes in order to estimate # of lines
fseek(infid, 0, 'eof');
filesz = ftell(infid);

% Assume 10 chars per line
% (an underestimate, but this just means
%  we over-estimate the number of lines)
charperline = 10;
Nl_guess = ceil(filesz/(charperline));
%disp(['Reading ',infile,': ',num2str(filesz,'%10.0f')]);
%disp(['(guesstimated ',num2str(Nl_guess),' lines)']);

% Pre-allocate arrays base on this guess
regnames = cell(Nl_guess,1);
regbounds = zeros(Nl_guess, 2);

% Read in file data
fseek(infid,0,'bof'); % Go back to beginning of file

% Skip comment lines
Nregs = 0;
while (feof(infid) == 0) 
    line = fgetl(infid);

    % Ignore commented or blank lines
    if ~((line(1) == ';') || (line(1) == '#') || ...
         (line(1) == '%') || (numel(line) == 0))
        
        Nregs = Nregs + 1;
        
        splitline = regexp(line, ',', 'split');
        
        vec0 = str2double(splitline{1}); 
        vec1 = str2double(splitline{2});
        
        cvecregs(vec0:vec1) = Nregs;
        
        regnames{Nregs} = splitline{3};
        regbounds(Nregs,:) = [vec0,vec1];
    end 
end

regnames = regnames(1:Nregs);
regbounds = regbounds(1:Nregs,:);

fclose(infid);


end

