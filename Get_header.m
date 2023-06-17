function header_info = Get_header(subject_number)
    filename = fullfile('Detrended',sprintf('SUBJ_%d',subject_number),'BLOCK_2.mat');
    load(filename)
    
    header_info.hdr = data.hdr;
    header_info.grad = data.grad;
    
end

