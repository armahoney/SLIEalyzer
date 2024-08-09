function slie_stackslieimages(varargin)

% Function to stack SLIE images over a specified time period to compute
% landfast ice occurrence frequency / probability
%
% USAGE
%      slie_stackslieimages
%      slie_stackslieimages('dataset', dataset)
%      slie_stackslieimages('asadir', asadir)
%      slie_stackslieimages('cfgfile', cfgfile)
%      slie_stackslieimages('mo2stack', mo2stack)
%      slie_stackslieimages('yr0', yr0)
%      slie_stackslieimages('yr1', yr1)
%
% INPUTS
%   dataset : string specifying which dataset to use.
%             See slie_get_basefolder.m for accepted values.
%             default: 'beau'
%   cfgfile : string specifying the name of the configuraion file stored
%             in the region folder containing essential information for
%             performing all analysis of SLIE data. 
%             See slie_readconfigfile.m for more details
%             default: 'slieconfig.txt'
%    asadir : string specifying folder name in region folder where
%             the data files from the all-seasons analysis are held
%             default: 'AllSeasonsAnalysis'
%  mo2stack : array of integers specifying which calendar months to include in the
%             in the SLIE image stack
%             default: [10:12, 1:7];
%       yr0 : integer specifying start year
%             default: first year in record;
%       yr1 : integer specifying end year
%             default: last year in record;
%
%
% Andy Mahoney - May 2024

% ------------------------------------------------------------------------
% Handle any variable arguments specified with function call
Nvarargs = numel(varargin);
a = 1;
while a <= Nvarargs
    switch varargin{a}
        case 'dataset'
            dataset = varargin{a+1};
            a = a + 2;
        case 'cfgfile'
            cfgfile = varargin{a+1};
            a = a + 2;
        case 'asadir'
            asadir = varargin{a+1};
            a = a + 2;
        case 'mo2stack'
            mo2stack = varargin{a+1};
            a = a + 2;
        case 'yr0'
            yr0 = varargin{a+1};
            a = a + 2;
        case 'yr1'
            yr1 = varargin{a+1};
            a = a + 2;
        otherwise
            disp([varargin{a} ': argument not recognized']);
            a = a  +1;
    end
end    
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set default values for parameters not specified in function call
if ~exist('dataset', 'var')
    dataset = 'beau';
end
% Identify the location where data are stored for specfied region
sfolder = slie_get_basefolder(dataset);

if ~exist('cfgfile', 'var')
    cfgfile = 'sliealyzer_config.txt';
end

% Set default events data folder if not spec'd in function call
if ~exist('asadir', 'var')
    asadir = 'AllSeasonsAnalysis2023';
end

% Specify which events we're plotting
if ~exist('mo2stack')
    mo2stack = [10:12, 1:7];
end
Nmo = numel(mo2stack);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Read the SLIE configuration file to get subregion and vector label info
cfg = slie_readconfigfile(sfolder, cfgfile);
Nseasons = numel(cfg.seasons);
if ~exist('yr0', 'var')
    yr0 = str2double(cfg.seasons{1}(1:4));
end
if ~exist('yr1', 'var')
    yr1 = str2double(cfg.seasons{end}(1:4))+1;
end



sliestack_all = zeros(cfg.sz);
sliestack_yrly = zeros(cfg.sz);
scount = 0;
for ss=1:Nseasons
    sliestackS = zeros(cfg.sz);
    sfiles = filedigger([sfolder filesep cfg.seasons{ss}], '*slie.tif');
    
    sdates = slie_name2date(sfiles);
    [syr,smo,~] = datevec(sdates);
    
    for s=1:numel(sfiles)

        momatch = any(smo(s) == mo2stack);
        yr0match = syr(s) >= yr0 + (smo(s) < 10);
        yr1match = syr(s) <= yr1 - (smo(s) >= 10);

        if ~momatch || ~yr0match || ~yr1match
            continue;
        end
            
        fprintf('Reading %s ...', sfiles{s});
        sliefile = [sfolder filesep cfg.seasons{ss} filesep sfiles{s}];
        [slieimg, R] = geotiffread(sliefile);

        sliestackS = sliestackS + double(slieimg == 255);
        scount = scount + 1;
        fprintf(' done.\n');
    end
    sliestack_all = sliestack_all + sliestackS;
    sliestack_yrly = sliestack_yrly + (sliestackS > 0);

    notalot = 0;
end

notalot = 0;

gtinfo = geotiffinfo(sliefile);
geokeyDT = gtinfo.GeoTIFFTags.GeoKeyDirectoryTag;
% Remove any "unknown" fields from geokey
if isfield(geokeyDT, 'Unknown')
    geokeyDT = rmfield(geokeyDT, 'Unknown');
end

% Write stack images to GeoTIFFs
outname = sprintf('SLIEimagestack_%04d-%04d_mo%02d-%02d_%s.tif', ...
                   yr0, yr1, mo2stack(1), mo2stack(end), dataset);
outtif = [sfolder filesep outname];
geotiffwrite(outtif, sliestack_all, R, 'GeoKeyDirectoryTag',geokeyDT);

% Write probability images to GeoTIFFs
outname = sprintf('SLIEimageprob_%04d-%04d_mo%02d-%02d_%s.tif', ...
                   yr0, yr1, mo2stack(1), mo2stack(end), dataset);
outtif = [sfolder filesep outname];
geotiffwrite(outtif, sliestack_all/scount, R, 'GeoKeyDirectoryTag',geokeyDT);


% Write yearly presence images to GeoTIFFs
outname = sprintf('SLIEyearlypres_%04d-%04d_mo%02d-%02d_%s.tif', ...
                   yr0, yr1, mo2stack(1), mo2stack(end), dataset);
outtif = [sfolder filesep outname];
geotiffwrite(outtif, sliestack_yrly, R, 'GeoKeyDirectoryTag',geokeyDT);

notalot = 0;

end

        
