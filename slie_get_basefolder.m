function sfolder = slie_get_basefolder(dataset)
% A simple function used at the beginning of a number of tools
% to identify the location where data are stored for different regions
%
% This simplifies editting the code to run on different machines 
% with data stored in different locations. 
%
% The intent is that the base paths below can be updated for whatever
% system the sliealyzer code is being run on
%
% USAGE
%   sfolder = slie_get_regionfolder(region)
%
% INPUT
%   region : a character array specifying the name of the desired
%            region / dataset to analyze
%
% OUTPUT
%  sfolder : the full path to the base folder or SLIE folder,
%            which serves as the base path for all the data and 
%            analsys that follows
%
% Andy Mahoney - October 2023

% ---------------------------------------------------------------------
switch dataset
    case 'fib'
        sfolder = '~/DATA//CMI-MMS_DATA/SLIEs/FoggyIslandBay';
    case 'beau'
        sfolder = '~/Google Drive/My Drive/Share/DATA/M2014_LFIdata/SLIEs/Beau';
    case 'chuk'
        sfolder = '~/Google Drive/My Drive/Share/DATA/M2014_LFIdata/SLIEs/Chuk';
    case 'beauM2014'
        sfolder = '~/Google Drive/My Drive/Share/DATA/M2014_LFIdata/SLIEs/Beau';
    case 'chukM2014'
        sfolder = '~/Google Drive/My Drive/Share/DATA/M2014_LFIdata/SLIEs/Chuk';
    case 'fibNIC'
        %sfolder = '~/DATA/NICcharts/FoggyIslandBay/GeoTIFFs_BOEM';
        sfolder = '~/DATA/NICcharts/GeoTIFFs_BOEM/FoggyIslandBay';
    case 'beauNIC'
        %sfolder = '~/DATA/NICcharts/Beau/GeoTIFFs_BOEM';
        sfolder = '~/DATA/NICcharts/GeoTIFFs_BOEM/Beau';
    case 'chukNIC'
        %sfolder = '~/DATA/NICcharts/Chuk/GeoTIFFs_BOEM';
        sfolder = '~/DATA/NICcharts/GeoTIFFs_BOEM/Chuk';
    case 'fibASIP'
        sfolder = '~/DATA/ASIP_IceCharts/BOEM_GeoTIFFs/FoggyIslandBay';
    case 'beauASIP'
        %sfolder = '~/DATA/ASIP_IceCharts/BOEM_GeoTIFFs/Beau';
        sfolder = '~/Google Drive/My Drive/Share/DATA/ASIP_IceCharts/BOEM_GeoTIFFs/Beau';
    case 'chukASIP'
        %sfolder = '~/DATA/ASIP_IceCharts/BOEM_GeoTIFFs/Chuk';    
        sfolder = '~/Google Drive/My Drive/Share/DATA/ASIP_IceCharts/BOEM_GeoTIFFs/Chuk';    
    case 'fibHilcorp'
        sfolder = '~/Documents/HillCorp/Sentinel1_SLIEs';
    case 'beauDaily'
        sfolder = '~/DATA//CMI-MMS_DATA/SLIEs/Beaufort_Daily';
    case 'chukDaily'
        sfolder = '~/DATA//CMI-MMS_DATA/SLIEs/Chukchi_Daily';
    case 'beauEM24'
        sfolder = '~/Google Drive/My Drive/Share/DATA/EM24/Beau';
    case 'chukEM24'
        sfolder = '~/Google Drive/My Drive/Share/DATA/EM24/Chuk';
    otherwise
        % Checks to see if dataset name is of the form
        % rrrrMOM6[rmxx]_mm_ttt_loc
        % where
        %    rrrr specifies region (either "chuk" or "beau")
        %    rmxx optionally specifies ax xx-day runnming minimum
        %      mm specifies model run (51 or 56)
        %     ttt specifies v threshold (01, 005, or 001)
        %     loc specifies location ("gdrive" or "local")
        if ~isempty(regexp(dataset, 'MOM6', 'once'))
            rexp = ['(?<reg>(beau)|(chuk))MOM6' ...
                    '(?<dly>(daily)?)_' ...
                    '(?<run>[0-9]+)_(?<vthresh>[0-9]+)' ...
                    '_*(?<loc>[a-z]*)'];                   
            d = regexp(dataset, rexp, 'names');
            Reg = [upper(d.reg(1)) d.reg(2:end)];
            
            if numel(d.dly) > 0
                dailysuff = '_Daily';
            else
                dailysuff = '';
            end

            % Construct model data folder name
            modeldir = ['ModelRun_' d.run '_Threshold_' d.vthresh];

            % Figure out where to look for this folder
            if numel(d.loc) == 0
                d.loc = 'gdrive';
            end
            switch d.loc
                case 'gdrive'
                    MOM6dir = '~/Google Drive/My Drive/Share/DATA/MOM6_FICE';
                case 'local'
                    MOM6dir = '~/DATA/MOM6_FICE_local';
                otherwise
                    disp([d.loc ': location not recognized. Trying Google Drive.']);
            end
            
            % Assemble this into a path to find SLIE data
            sfolder = [MOM6dir filesep modeldir filesep ...
                       'BOEM_GeoTIFFs' filesep Reg dailysuff];
        else
            disp(['Region not recognized: ', dataset]);
            disp('Quitting.');
            return
        end
end

end