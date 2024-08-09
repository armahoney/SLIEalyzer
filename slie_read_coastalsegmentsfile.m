function csegs = slie_read_coastalsegmentsfile(csegfile)
% Function to read textg file defining coastal segments
% and how to close them for rasterizing landfast ice width data
% See lfiw2raster


infid = fopen(csegfile, 'rt');

s = 0;
while ~feof(infid)
    line = fgetl(infid);
    splitline = strtrim(regexp(line, ':', 'split'));
    switch true
        case strcmp(splitline{1}, 'Segment')
            s = s + 1;
            csegs(s).segname = splitline{2};
        case regexp(splitline{1}, 'v[01]', 'once')
            csegs(s).(splitline{1}) = str2double(splitline{2});            
        case strcmp(splitline{1}, 'direction')
            csegs(s).(splitline{1}) = splitline{2};            
        case regexp(splitline{1}, '[xy]close', 'once')
            listsplit = regexp(splitline{2}, ',', 'split');
            csegs(s).(splitline{1}) = str2double(listsplit)'+1;
        case strcmp(splitline{1}, 'fillxy')
            listsplit = regexp(splitline{2}, ',', 'split');
            csegs(s).(splitline{1}) = str2double(listsplit);
        case strcmp(splitline{1}, 'overlap')
            csegs(s).(splitline{1}) = splitline{2};            
        case strcmp(splitline{1}, 'abut0')
            csegs(s).(splitline{1}) = splitline{2};            
        case strcmp(splitline{1}, 'abut1')
            csegs(s).(splitline{1}) = splitline{2};            
        otherwise
            notalot = 0;
    end
end


fclose(infid);
    
end
