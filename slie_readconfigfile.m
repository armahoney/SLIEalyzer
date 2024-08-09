function cfg = slie_readconfigfile(spath, cfgfile)
% Function to read key information from SLIEalyzer configuration file
%
% USAGE
%   cfg = slie_readconfigfile(spath)
%   cfg = slie_readconfigfile(spath, cfgfile)
%
%
% INPUTS
%     spath : string specifying full path to folder in which the
%             configuration file is held. This is referred elsewhere
%             as the SLIE path or base path
%   cfgfile : string specifying name of the configuration file to use
%             Default: 'sliealyzer_config.txt'
%
% OUTPUT
%       cfg : a structure containing the following fields:
%           -       localpaths: ** OPTIONAL ** text flag indicating whether
%                               or not to treat filepaths as local to the
%                               specified sfolder
%                               default: 'FALSE'
%           -       slieprefix: character array of letters that may be used
%                               at the beginning of value SLIE image files
%           -         cvecfile: path specifying ASCII file containing
%                               coastvector coordinates
%           -        lmaskfile: path specifying GeoTIFF containing 
%                               landmask image
%           -    lmask_ni_file: path specifying GeoTIFF containing 
%                               island-free landmask image
%           -  headbaymaskfile: path specifying non-GeoTIFF containing 
%                               landmask image with headlands and bays
%                               encoded
%           -   lagoonmaskfile: path specifying non-GeoTIFF containing 
%                               landmask image with lagoons encoded
%           -        bathyfile: path specifying GeoTIFF containing 
%                               bathymetric grid
%           -        cdistfile: path specifying GeoTIFF containing 
%                               coastal distance grid
%           -               sz: dimensions of all images used for this
%                               region. 
%                               NOTE - this is specified in column-row
%                                      order unlike most Matlab stuff
%                               
%           -            pixsz: size of (square) pixels
%           -          min_ice: minimum ice width threshhold for
%                               determining occurence of first and last ice
%           -         stab_dep: water depth defining stablity
%           -   cveclabelsfile: a list of coast vector indices and labels
%           -  cvecregionsfile: a list of beginning and ending coast vector
%                               indices together with sub region names
%           - studyoutlinefile: a file defining the study area outline
%           -      coastxyfile: a file proving x-y coordinates for the
%                               coastline
%           -          seasons: a list of seasons to be processed. Where
%                               each seasnon is defined by a string of
%                               the form yyyy-yy
%           - coastsegmentfile: path specifying a text file
%                               containing information on how to segment
%                               the coastline for the purposes of
%                               rasterizing lfiw data
%
% Andy Mahoney - October 2023
%              - Jan 2024: updated to allow specification of config data
%                          folder, with all other paths then being local
%                          to that folder

% Assign default config filename if not specified in function call
if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

% Exit if no config file found
if ~exist([spath filesep cfgfile], 'file')
    disp(['No config file found in ',spath]);
    cfg = -1;
    return
end

% Open the config file for reading as text
infid = fopen([spath filesep cfgfile], 'rt');

% Read through line-by-line until the end
% - ignore blanks lines those beginning with ';', '#', or '%'
while (feof(infid) == 0) 
    line = fgetl(infid);
    
    % Ignore commented or blank lines
    if ~((line(1) == ';') || (line(1) == '#') || ...
         (line(1) == '%') || (numel(line) == 0))
                
        % Separate parts either side of equals sign
        splitline = regexp(line, ' = ', 'split');
        
        % Exit if badly formatted line is found
        if numel(splitline) ~= 2
            disp(['Bad line in config file: ', line]);
            cfg = -1;
            fclose(infid);
            return
        end

        % Assign the left and right sides to "var" and "val"
        var = splitline{1};
        val = splitline{2};
        
        % Now assign output fields accordingly
        switch var
            case 'localpaths'
                cfg.localpaths=val;
            case 'slieprefix'
                cfg.slieprefix=val;
            case 'lmaskfile'
                cfg.lmaskfile=val;
            case'lmaskfile_noisl'
                cfg.lmask_ni_file=val;
            case'headbaymaskfile'
                cfg.headbaymaskfile=val;
            case'lagoonmaskfile'
                cfg.lagoonmaskfile=val;
            case 'cvecfile'
                cfg.cvecfile = val;
            case 'bathyfile'
                cfg.bathyfile = val;
            case 'cdistfile'
                cfg.cdistfile = val;
            case 'sz'
                sz = str2double(regexp(val, ',', 'split'));
                cfg.sz = sz([2,1]);  
            case 'pixsz'
                cfg.pixsz = str2double(val);
            case 'min_ice'
                cfg.min_ice = str2double(val);
            case 'stab_dep'
                cfg.stab_dep = str2double(val);
            case 'cveclabelsfile'
                cfg.cveclabelsfile = val;
            case 'cvecregionsfile'
                cfg.cvecregionsfile = val;
            case 'studyoutlinefile'
                cfg.studyoutlinefile = val;
            case 'coastxyfile'
                cfg.coastxyfile = val;
            case 'seasons'
                cfg.seasons = strtrim(regexp(val, ',', 'split'));
            case 'coastsegmentfile'
                cfg.coastsegmentfile = val;
            otherwise
                disp(['Unrecognized variable in config file: ', var]);
        end
    end 
end

% Assign a default slieprefix value if not specified in config file
if ~isfield(cfg, 'slieprefix')
    cfg.slieprefix = 'r';
end

% If localpaths is specified as TRUE, use sfolder as prefix for all 
% file path designations
if ~isfield(cfg, 'localpaths')
    cfg.localpaths = 'FALSE';
end

if strcmp(cfg.localpaths, 'TRUE')
    cfgflds = fieldnames(cfg);
    for f=1:numel(cfgflds)
        if contains(cfgflds{f}, 'file')
            cfg.(cfgflds{f}) = [spath filesep cfg.(cfgflds{f})];
        end
    end
end


% Close the input file
fclose(infid);

end

