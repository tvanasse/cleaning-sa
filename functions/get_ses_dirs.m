function [inputlist] = get_ses_dirs()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

more_files = 'Yes';
first_iter = 0;
while strcmp(more_files, 'No') ~= 1
if first_iter == 0;
    inputlist = uigetfile_n_dir([pwd '/../']);
else
    inputlist = [inputlist uigetfile_n_dir([pwd '/../'])];
end

celldisp(inputlist)

more_files = questdlg('Add more session folders?', ...
	'??', ...
	'Yes','No','No');
first_iter = first_iter + 1;
end

end

