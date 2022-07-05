folder = 'E:\combined_FICs_for_each_yeast_cell';
files = dir(folder);
filesNames = {files.name};

for fileNameCell = filesNames
    fileName = fileNameCell{1};
    if strcmp(fileName, '.') || strcmp(fileName, '..')
        continue
    end
    
    load(fileName);
    file_subset;
    mat_struct;
    fileName = [fileName(1:end-3), 'xlsx'];
    for k = 1:length(mat_struct)
        force_a = mat_struct(k).force_a.';
        absci_a = mat_struct(k).absci_a.';
        force_r = mat_struct(k).force_r.';
        absci_r = mat_struct(k).absci_r.';
        T = table(force_a, absci_a, force_r, absci_r);
        
        vitesse = mat_struct(k).vitesse;
        cell_number = mat_struct(k).cell_number;
        
        writecell({file_subset}, fullfile(folder, fileName), 'Sheet', k, 'Range', 'A1');
        writetable(table(cell_number, vitesse), fullfile(folder, fileName), 'Sheet', k, 'Range', 'A3');
        writetable(T, fullfile(folder, fileName), 'Sheet', k, 'Range', 'A6');
    end
end