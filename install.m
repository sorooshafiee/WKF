%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Distributionally Robust Kalman Filter
% Soroosh Shafieezadeh-Abadeh, Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Install function
% Add the src and utils folder into the working directory of matlab


path1 = pwd;
path_src = strcat(path1, filesep, 'src');
addpath(genpath(path_src)); 
path_utils = strcat(path1, filesep, 'utils');
addpath(genpath(path_src)); 
clear path1 path_src path_utils
