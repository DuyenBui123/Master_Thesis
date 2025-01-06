clc
clear
close all
cd 'C:/Master_Thesis/'
root_dir = 'C:/Master_Thesis/';
%% 

data_dir = sprintf('%COLD/Data/rstudio-export/', root_dir);
addpath(data_dir)
code_dir = sprintf('%sGitHub_code/COLD',root_dir);
addpath(code_dir)
code_dir_extra = sprintf('%sCOLD/Code',root_dir);
addpath(code_dir_extra)
%% 
main_COLD_ori(1,1)