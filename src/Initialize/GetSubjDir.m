% GetSubjDir.m ------------------------------------------------------------
% This function prompts the user to select the directory containing the
% data to be processed. The user is then prompted to select the
% participants to be processed from a list of all participant folders in
% the data directory. The function returns a structure containing the
% handles to all selected subject directories.
%
% The function requires that participants' raw motion capture data are
% stored in one data folder with each participant's individual data
% contained in one subfolder per participant. It is recommended that the
% participant identifiers are used to name these subfolders for clarity.
%                                                                           
%--------------------------------------------------------------------------
% Author: Torstein Daehlin, PhD. August, 2024.                              
%--------------------------------------------------------------------------

function [subj_dir, data_path, start_idx] = GetSubjDir()
% Get data path
start_path = regexp(pwd, filesep, 'split');
start_path = fullfile(start_path{1:end-2});

data_path = uigetdir(start_path, 'Select data directory');

% Get handle to path contents
subj_dir = dir(data_path);
subj_dir = subj_dir(~contains({subj_dir.name},'.'));

% Select trials of interest
subj_idx = listdlg('ListString',{subj_dir.name});

% Shorten directory list and return
subj_dir = subj_dir(subj_idx);

% Output index of first selected participant
start_idx = subj_idx(1);
end