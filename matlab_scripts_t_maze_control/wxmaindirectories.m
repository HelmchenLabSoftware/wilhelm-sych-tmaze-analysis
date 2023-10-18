%% wxmaindirectories
 % Assigns main directory paths
 % JL Alatorre Warren

function [pathRoot, ...
          pathData, ...
          pathRois, ...
          pathTrials, ...
          pathTrialList, ...
          pathSmoothStop] = wxmaindirectories
 
pathRoot           = 'I:';
pathData           = [pathRoot   '/' '2017_m001-m006_backup'];
pathRois           = [pathData   '/' 'rois'];
pathTrials         = [pathData   '/' 'trials'];
pathTrialList      = [pathTrials '/' 'list'];
pathSmoothStop     = [pathData   '/' 'smoothstop'];