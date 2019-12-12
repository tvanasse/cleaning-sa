function mff_copy_info_wo_calibration(oldfile,newfile)

oldid      = fopen(oldfile, 'r');
newid      = fopen(newfile, 'w');

section = '';

while ~feof(oldid)
    
    tline = fgets(oldid);
    
    [variables] = regexp(tline, '<(?<section>calibrations)( .+)?>', 'names');
    
    if ~isempty(variables)
        section = variables.section;
        
        continue;
    end
    
    [variables] = regexp(tline, '</(?<section>calibrations)>', 'names');
    
    if ~isempty(variables)
        section = '';
        
        continue;
    end
    
    switch section
        case 'calibrations'
            continue
        otherwise
        %    fseek(newid,0,'eof');
            fwrite(newid,tline);
    end  
end

fclose(oldid);
fclose(newid);

end