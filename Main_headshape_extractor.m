% This function is designed to read the fiducials from every subject file.
% Since these are different for each block, we propose to average them in
% order to get a common representation for our timelock data.
% In the end of this script a file is written for each of the subjects
% containing such average fiducial information into the following folder:
% ../Results/Subject_fiducials/
%
% Visitor: 
% Antonio Rodriguez Hidalgo 
% Dept. of Signal Theory and Communications
% Universidad Carlos III de Madrid
% arodh91@gmail.com
%
% Principal Investigator:
% Maria Chait 
% Ear Institute
% University College London
% m.chait@ucl.ac.uk
%
% Last update: 10/January/2020

clear
addpath(genpath('D:\fieldtrip-20190819'))

subject_list = [2:13, 15];
% subject_list2 = [18:24];

for subject = subject_list

    output_path = fullfile('..','Results','Subject_fiducials');
%     mkdir(output_path)

    folder = fullfile('..','MEGGAP','data_longshort',sprintf('Subj%d', subject));
    
    % First, we find all the available files. In order to do so, we focus
    % exclusively on getting the names of the folders where MEG data is
    % stored.
    files_global = dir(folder);
    files_global = files_global(3:end);    
    
    files = {};
    c = 1;
    for ind_files = 1:length(files_global)
        if files_global(ind_files).isdir == true
            files{c} = files_global(ind_files).name;
            c = c+1;
        end
    end
    
    n_files = length(files);

    % Once we have the folder names, we read the headshape from each of
    % them and get their average fiducials. Notice that there might be
    % unreadable files, in which case we will not consider it in the
    % averaging process. We consider that units and coordinates are common
    % for every file.

    average_headshape = [];
    for ind = 1:n_files
        file_path = files{ind};
        try
            headshape = ft_read_headshape(sprintf('../MEGGAP/data_longshort/Subj%d/%s', subject, file_path));
            if ind == 1
                average_headshape.fid.pos = headshape.fid.pos;
            else
                average_headshape.fid.pos = average_headshape.fid.pos + headshape.fid.pos;
            end    
            
        catch
            n_files = n_files-1;
            fprintf('Problem reading fidutials.SUBJ:%d. FILE:%d. FILENAME:%s\n', subject, ind, file_path);
        end
    end
    
    average_headshape.fid.pos = average_headshape.fid.pos/n_files;
    average_headshape.pos = headshape.pos;
    average_headshape.fid.label = headshape.fid.label;
    average_headshape.coordsys = headshape.coordsys;
    average_headshape.unit = headshape.unit;

    
    % We also get the standard deviation, but just for fun (it doesn't get
    % stored).
    for ind = 1:n_files
        file_path = files{ind};
        try
            headshape = ft_read_headshape(sprintf('../MEGGAP/data_longshort/Subj%d/%s', subject, file_path));
            if ind == 1
                stdev = pow2(headshape.fid.pos-average_headshape.fid.pos);
            else
                stdev = stdev + pow2(headshape.fid.pos-average_headshape.fid.pos);
            end    
            
        catch
            n_files = n_files-1;
            fprintf('Problem reading fidutials.SUBJ:%d. FILE:%d. FILENAME:%s\n', subject, ind, file_path);
        end
    end
    
    stdev = sqrt(stdev/(n_files-1));

    save(fullfile(output_path, sprintf('Fid_Subj%d.mat', subject)), 'average_headshape')

    clear average_headshape subject_structure file_path output_path files
    
end
clear
