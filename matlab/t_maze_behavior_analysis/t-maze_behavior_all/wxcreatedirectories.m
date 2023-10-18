%% wxcreatedirectories
 % Creates trial-specific directories
 % JL Alatorre Warren

function [pathRoot, ...
          pathData, ...
          pathRois, ...
          pathTrials, ...
          pathTrialList, ...
          pathCurrentTrial, ...
          pathMouseName, ...
          pathDateString, ...
          pathExperimentName, ...
          pathTimeString, ...
          pathFlags, ...
          pathSignals, ...
          pathFrames, ...
          pathLogs, ...
          pathSmoothStop] = wxcreatedirectories(timestampVector, ... 
                                                mouseName, ...
                                                dateString, ...
                                                experimentName, ...
                                                timeString)

% Trial name
trialName = ['trial_' datestr(timestampVector,'yyyymmddHHMMSSFFF')];        

% Get main directories
[pathRoot, ...
 pathData, ...
 pathRois, ...
 pathTrials, ...
 pathTrialList, ...
 pathSmoothStop] = wxmaindirectories;

% Assign additional directory paths
pathCurrentTrial   = [pathTrialList '/' trialName];
pathMouseName      = [pathTrials    '/' mouseName];
pathDateString     = [pathTrials    '/' mouseName '/' dateString];
pathExperimentName = [pathTrials    '/' mouseName '/' dateString '/' experimentName];
pathTimeString     = [pathTrials    '/' mouseName '/' dateString '/' experimentName '/' timeString];
pathFlags          = [pathTrials    '/' mouseName '/' dateString '/' experimentName '/' timeString '/' 'flags'];
pathSignals        = [pathTrials    '/' mouseName '/' dateString '/' experimentName '/' timeString '/' 'signals'];
pathFrames         = [pathTrials    '/' mouseName '/' dateString '/' experimentName '/' timeString '/' 'frames'];
pathLogs           = [pathTrials    '/' mouseName '/' dateString '/' experimentName '/' timeString '/' 'logs'];

% Create directory: current trial                                
if ~exist(pathCurrentTrial, 'dir')
  mkdir(pathCurrentTrial);
end

% Create directory: mouse
if ~exist(pathMouseName, 'dir')
  mkdir(pathMouseName);
end

% Create directory: date (of the trial)
if ~exist(pathDateString, 'dir')
  mkdir(pathDateString);
end

% Create directory: experiment type
if ~exist(pathExperimentName, 'dir')
  mkdir(pathExperimentName);
end

% Create directory: time of the day
if ~exist(pathTimeString, 'dir')
  mkdir(pathTimeString);
end

% Create directory: flags
if ~exist(pathFlags, 'dir')
  mkdir(pathFlags);
end

% Create directory: signals
if ~exist(pathSignals, 'dir')
  mkdir(pathSignals);
end

% Create directory: frames
if ~exist(pathFrames, 'dir')
  mkdir(pathFrames);
end

% Create directory: logs
if ~exist(pathLogs, 'dir')
  mkdir(pathLogs);
end

% Create directory: smooth stop
if ~exist(pathSmoothStop, 'dir')
  mkdir(pathSmoothStop);
end