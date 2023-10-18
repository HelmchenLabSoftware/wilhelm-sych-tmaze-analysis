%% wxtmaze
 % Main program for the automation of the T-maze
 % JL Alatorre-Warren

%% Initial setup

% Request user input
[mouseName, ...
 experimentName, ...
 maxNumberOfRunsShapingPeriod, ...
 maxNumberOfRunsLearningPeriod, ...
 secondsInsideStartBoxAfterForcedRun, ...
 secondsInsideStartBoxAfterFreeRun, ... 
 optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
 optionWhiteNoiseInsideStartBoxAfterFreeRun, ...  
 optionUseStartGateAfterFreeingTheMouse, ...
 optionLickToOpenValve] = wxuserinput;

% Create main logbook
mainLogbook = cell(50*(maxNumberOfRunsShapingPeriod + ...
                       maxNumberOfRunsLearningPeriod),5);

% Alpha: start time tracking
timestampVector = clock;
dateString = datestr(timestampVector,'yyyymmdd');
timeString = datestr(timestampVector,'HHMMSSFFF');
mainLogbook = wxupdatelogbook(mainLogbook, 'alpha', timestampVector);

% Directory paths, create new trial, and save relevant variables
[pathRoot, ...
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
 pathLogs] = wxcreatedirectories(timestampVector, ... 
                                 mouseName, ...
                                 dateString, ...
                                 experimentName, ...
                                 timeString); %#ok<*ASGLU>

% Save relevant variables to pass to wxrecordsignals and wxrecordvideo
save([pathCurrentTrial '/' 'variables.mat'], ...
     'mouseName', ...
     'experimentName', ...
     'timestampVector', ...
     'dateString', ...
     'timeString');
   
%% DIPimage and DIPlib
%{
  The function |dipstart| initializes DIPimage (MATLAB)
  DIPimage is a MATLAB toolbox for scientific image processing and
  analysis.
  Most image processing operations are deferred to DIPlib, a dedicated
  library written in C. In this way we extend MATLAB with new core
  functions, and thereby overcoming speed limitations of interpreted
  scripts.
  DIPlib is a platform independent scientific image processing library
  written in C. It consists of a large number of functions for processing
  and analysing multi-dimensional image data. The library provides
  functions for performing transforms, filter operations, object
  generation, local structure analysis, object measurements and
  statistical analysis of images.
%}
dipstart;
   
%% National Instruments DAQ (Multifunction Data Acquisition)
%  National Instruments USB-6210

% Start DAQ session
[daqDevice, daqSession] = wxnibox6210;
 
%% Check water valves

% Check valves and load water in them
[daqDevice,daqSession] = wxcheckvalves(daqDevice,daqSession);

%% Test Phidgets's Advanced Servo Motor
%  1061_1 - PhidgetAdvancedServo 8-Motor

% Open and close all servos |numberOfCycles| times
numberOfCycles = 2;
mainLogbook = wxtestservos(numberOfCycles,mainLogbook);
clear numberOfCycles

%% Test Lumenera camera

% Assign camera number
% Display video preview
cameraNumber = 1; 
LucamShowPreview(cameraNumber);

%% Regions of interest

% Open and display regions of interest
[roi01, infoRoi01, ...
 roi02, infoRoi02, ...
 roi03, infoRoi03, ...
 roi04, infoRoi04] = wxgetrois(cameraNumber, pathRois);

%% Parallel MATLAB sessions

% Start additional MATLAB sessions
 !matlab -nosplash -minimize -r "wxrecordvideo" &
 !matlab -nosplash -minimize -r "wxrecordsignals" &
 !matlab -nosplash -minimize -r "wxdisplaysignals" &

% Wait until all flags are found
flagsFound = false;
while flagsFound == false
  flagDisplaySignals = exist([pathFlags '/' 'flagdisplaysignalsready.mat'],'file');
  flagRecordSignals  = exist([pathFlags '/'  'flagrecordsignalsready.mat'],'file');
	flagRecordVideo    = exist([pathFlags '/'    'flagrecordvideoready.mat'],'file');
  if (flagRecordSignals == 2) && (flagDisplaySignals == 2) && (flagRecordVideo == 2)
    flagsFound = true;
  else
    pause(1)
    disp(['Waiting for the flags: ' datestr(now,'HH:MM:SS')])
  end
end

%% Launch status check

% Check valves and launch the T-maze sequence
wxlaunchstatuscheck;

%% Sequence 01: shaping period

if maxNumberOfRunsShapingPeriod == 0
  counterOfRuns = 1;
else
  counterOfRuns = 0;
end

while counterOfRuns < maxNumberOfRunsShapingPeriod
  
  % Display trial number
  disp(['Welcome to trial number ' num2str(counterOfRuns+1)])
  
  for jj = 1:2
    
    if jj == 1
      currentSide = 'right';
    elseif jj == 2
      currentSide = 'left';
    end
    
    % Shaping sequence
    [daqDevice, ...
     daqSession, ...
     mainLogbook, ...
     startTimeForcedRun, ...
     endTimeForcedRun] = wxshapingsequence(currentSide, ...
                                           cameraNumber, ...
                                           mainLogbook, ...
                                           daqDevice, ...
                                           daqSession, ...
                                           roi01, ...
                                           roi02, ...
                                           secondsInsideStartBoxAfterForcedRun, ...
                                           optionWhiteNoiseInsideStartBoxAfterForcedRun, ...  
                                           optionUseStartGateAfterFreeingTheMouse, ...
                                           optionLickToOpenValve);
      
  end
  
  counterOfRuns = counterOfRuns+1;
  
  % Omega: end of full run
  % Save mainLogbook
  if counterOfRuns < 10
    counterStamp = ['0' num2str(counterOfRuns)]; 
  else
    counterStamp = num2str(counterOfRuns);
  end
  timestampVector = clock;
  mainLogbook = wxupdatelogbook(mainLogbook, ['omega' '_' counterStamp], timestampVector);
  save([pathLogs '/' 'logs.mat'], 'mainLogbook');         
  
end
clear jj

%% Sequence 02: learning period

% Generate pseudorandom sequence of sides
sequenceOfSides = wxgeneraterandomsequence(maxNumberOfRunsLearningPeriod);

counterOfRuns = 0;
while counterOfRuns < maxNumberOfRunsLearningPeriod
  
  % Calculate and display trial number
  counterOfRuns = counterOfRuns+1;
  disp(['Welcome to trial number ' num2str(counterOfRuns)])
  
  if sequenceOfSides{counterOfRuns,1} == 'R'
    currentSide = 'right';
  elseif sequenceOfSides{counterOfRuns,1} == 'L'
    currentSide = 'left';
  end
  
  % Learning sequence
  [daqDevice, ...
   daqSession, ...
   mainLogbook, ...
   sequenceOfSides, ...
   startTimeForcedRun, ...
   endTimeForcedRun, ...
   startTimeFreeRun, ...
   endTimeFreeRun] = wxlearningsequence(sequenceOfSides, ...
                                        mainLogbook, ...
                                        currentSide, ...
                                        counterOfRuns, ...
                                        cameraNumber, ...
                                        daqDevice, ...
                                        daqSession, ...
                                        roi01, ...
                                        roi02, ...
                                        roi03, ...
                                        roi04, ...
                                        secondsInsideStartBoxAfterForcedRun, ...
                                        secondsInsideStartBoxAfterFreeRun, ...
                                        optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
                                        optionWhiteNoiseInsideStartBoxAfterFreeRun, ...  
                                        optionUseStartGateAfterFreeingTheMouse, ...
                                        optionLickToOpenValve);
                                                              
  % Omega: end of full run
  if counterOfRuns < 10
    counterStamp = ['0' num2str(counterOfRuns)];
  else
    counterStamp = num2str(counterOfRuns);
  end
  timestampVector = clock;
  mainLogbook = wxupdatelogbook(mainLogbook, ['omega' '_' counterStamp], timestampVector);
  
  % Update trial times
  sequenceOfSides = wxupdatetrialtimes(sequenceOfSides, ...
                                       counterOfRuns, ...
                                       startTimeForcedRun, ...
                                       endTimeForcedRun, ...
                                       startTimeFreeRun, ...
                                       endTimeFreeRun);
                                     
  % Clear current trial times
  clear startTimeForcedRun
  clear endTimeForcedRun
  clear startTimeFreeRun
  clear endTimeFreeRun
  
  % Save mainLogbook and sequenceOfSides
  save([pathLogs '/' 'logs.mat'], 'mainLogbook', 'sequenceOfSides');
  
  % Update figure with mouse choices and trial times
  wxvisualizechoicesandtrialtimes(sequenceOfSides, counterOfRuns, maxNumberOfRunsLearningPeriod)
                                       
end
clear counterStamp
clear timestampVector

%% End MATLAB sessions

% Wait 10 seconds
pause(10)

% Send flag to wxrecordvideo, wxrecordsignals, and wxdisplaysignals
flagSequenceCompleted = true; %#ok<NASGU>
save([pathFlags '/' 'flagsequencecompleted.mat'], 'flagSequenceCompleted');

% Close MATLAB session
exit