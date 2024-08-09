function sliedates = slie_calcsliedates(lfiw, mosdates)
% Determines the appropriate mosaic date to assign to each SLIE point
% according to whether the SLIE was advancing or retreating at the time

% Some along-coast smoothing is done so that advance/retreat is considered on
% the basis of neighbouring rather than individual vectors

% NOTE: This was the approach taken in the analysis described by
%       Mahoney et al (2007; 2014), but inconsistencies and innaccuracies
%       were found and we no longer recommend this "time shifting"
%       Specifically, we found a small number of cases where the lfiw
%       drops to a minimum that lasts fewer than 3 consecutive SLIEs. 
%       Such short-lived minima should not occur if each SLIE truly
%       represents the minimu landfast ice extent in 3 consecutive
%       mosaics. Hence, these occurences demonstrate some level of
%       inconsitency in our original analysis. 
%       Also, the approach requires some modification when applied to 
%       data from ice charts, which are not explicitly derived from 
%       three consecutive epochs of imagery
%
% USAGE
%    sliedates = slie_calcsliedates(lfiw, mosdates)
%
% INPUTS
%       lfiw : a Nvec-by-Nslie array of landfast ice width measuremets
%   mosdates : a Nslie+2 list of dates corresponding to the SAR mosaics
%              used to delineate the SLIEs
%
% OUTPUT
%  sliedates : a Nvec-by-Nslie array of dates for each lfiw data point
%
%
% Andy Mahoney - October 2023

% ----------------------------------------------------------------------

% Get dimensions from input data
datasz = size(lfiw);
Nslies = datasz(2);
Nvecs = datasz(1);

% Pre-assign array to hold SLIE dates
sliedates = zeros(Nvecs,Nslies);

% Assign criteria for shifting dates
dif_thresh = 1;
smwid = 10;
sm_kern = zeros(smwid,1)+1.0/smwid;

% Determine the difference in width over time at each coast vector
dw = [lfiw(:,2)-lfiw(:,1), ...                    % First SLIE diff (forward looking)
      lfiw(:,3:Nslies)-lfiw(:,1:Nslies-2), ...    % Subsequent SLIE diffs (centred)
      lfiw(:,Nslies)-lfiw(:,Nslies-1)];           % Last SLIE diff (back looking)

% Smooth along the coast (i.e. cvec dimension)
dw_sm = conv2(dw, sm_kern, 'same');


% Determine whether SLIE was retreating, advancing or staying still
% and assign index to approipriate parent mosaic
% - assign to first mosaic (pmi = 1) if ice is advancing  (dw_sm +ve and > thresh)
% - assign to last mosaic  (pmi = 3) if ice is retreating (dw_sm -ve and > thresh)
% - assign to middle mosaic (pmi = 2) otherwise
pmi = 2 + 1*(dw_sm < -dif_thresh) - 1*(dw_sm >= dif_thresh);
 

% Calculate the mid-way date between 
%   1) the first and middle 
%   2) the first and last 
%   3) the middle and last
middates = zeros(Nvecs,Nslies,3);
middates(:,:,1) = repmat((mosdates(:,1)+mosdates(:,2))'/2, Nvecs,1);
middates(:,:,2) = repmat((mosdates(:,1)+mosdates(:,3))'/2, Nvecs,1);
middates(:,:,3) = repmat((mosdates(:,2)+mosdates(:,3))'/2, Nvecs,1);
 

% 2-D array of slie index numbers corresponding to each element in pmi
slie_i = repmat(1:Nslies, Nvecs, 1);
cvec_i = repmat((1:Nvecs)', 1, Nslies);
ii = sub2ind([Nvecs,Nslies,3], cvec_i, slie_i, pmi); 
sliedates = middates(ii);



end
