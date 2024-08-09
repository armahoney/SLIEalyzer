function slie_writemisssingpolygondata(mp, seasondir, sliefiles)
% A short function to output the locations of missing polyons identified
% from landfast ice width (lfiw) data
% See slie_findmissingpolygondata.m for more info
%
% USAGE
%       slie_writemisssingpolygondata(mp, seasondir, sliefiles)
%
% INPUTS
%        mp : Nvec-by-Nslie array of 0s and 1s indicating whether the
%             corresponding lfiw data point is the result of a missing
%             polygon (1 = missing)
%             This is the output from slie_findmissingpolygondata.m
% seasondir : String specifying the full path for the SLIE season
%             in which any missing polygons were found
% sliefiles : A cell array containg the filenames of the SLIE images
%             from which the landfast ice widths were calculated
%             
%
% OUTPUTS
%         An comma-separated ASCII file named 'missingpolygons.csv', 
%         which is created within the folder specified by seasondir.
%         Each row of this file provides a SLIE filenames followed by
%         two integers specifying the indices of the coast vectors that
%         bound the geographic region of a unique breakout event
         

% Andy Mahoney - October 2023
%
% -----------------------------------------------------------------------
outfid = fopen([seasondir filesep 'missingpolygons.csv'], 'wt');

fprintf(outfid, 'SLIEfile, cvec0, cvec1\n');
outfstr = '%s, %4.0f,%4.0f\n';

slieout = find(sum(mp,1) > 0);
for s=1:numel(slieout)
    ss = slieout(s);
    mpi0 = find([mp(1,ss); diff(mp(:,ss))] == 1); 
    mpi1 = find([diff(mp(:,ss)); -mp(end,ss)] == -1); 

    Nmpi = numel(mpi0);
    for i = 1:Nmpi
        fprintf(outfid, outfstr, sliefiles(ss).name, mpi0(i), mpi1(i));
    end
end

fclose(outfid);

end