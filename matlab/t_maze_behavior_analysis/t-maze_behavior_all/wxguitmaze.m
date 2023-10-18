%% wxguitmaze
 % Main program for the automation of the T-maze
 % JL Alatorre-Warren

function [sequenceOfSides, mainLogbook] = wxguitmaze(mouseName, ...
                                                     experimentName, ...
                                                     maxNumberOfRunsShapingPeriod, ...
                                                     maxNumberOfRunsLearningPeriod, ...
                                                     secondsInsideStartBoxAfterForcedRun, ...
                                                     secondsInsideStartBoxAfterFreeRun, ... 
                                                     optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
                                                     optionWhiteNoiseInsideStartBoxAfterFreeRun, ...  
                                                     optionUseStartGateAfterFreeingTheMouse, ...
                                                     optionLickToOpenValve, ...
                                                     optionSaveFrames, ...
                                                     optionSaveVideo, ...
                                                     laserSelection, ...
                                                     imageFormat) %#ok<INUSD,INUSL>

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
 pathLogs, ...
 pathSmoothStop] = wxcreatedirectories(timestampVector, ... 
                                       mouseName, ...
                                       dateString, ...
                                       experimentName, ...
                                       timeString);

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
 roi04, infoRoi04, ...
 roi06, infoRoi06] = wxgetrois(cameraNumber, pathRois);

%% Parallel MATLAB sessions

% Save flag variables to be used in wxrecordvideo
save([pathFlags '/' 'optionsaveframes.mat'], 'optionSaveFrames');
save([pathFlags '/' 'optionsavevideo.mat' ], 'optionSaveVideo' );
save([pathFlags '/' 'imageformat.mat'     ], 'imageFormat'     );

% Start additional MATLAB sessions
 !matlab -nosplash -minimize -r "wxrecordvideo" &
 !matlab -nosplash -minimize -r "wxrecordsignals" &
 !matlab -nosplash -minimize -r "wxdisplaysignals" &
 !matlab -nosplash -minimize -r "wxlasermodulation" & 

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

% Create cell array sequenceOfSides
sequenceOfAlternation = cell(maxNumberOfRunsShapingPeriod,13);

if maxNumberOfRunsShapingPeriod == 0
  counterOfRuns = 1;
else
  counterOfRuns = 0;
end

% Set inhibition status as 'off' (laser off)
currentInhibitionStatus = 'off';

while counterOfRuns < maxNumberOfRunsShapingPeriod
  
  % Display trial number
  disp(['Welcome to trial number ' num2str(counterOfRuns+1)])
  
  for jj = 1:2
    
    % This creates the alternation sequence: R, L, R, L, etc.
    if jj == 1
      currentSide = 'right';
      sequenceOfAlternation{counterOfRuns+1,1} = 'R';
    elseif jj == 2
      currentSide = 'left';
      sequenceOfAlternation{counterOfRuns+1,1} = 'L';      
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
                                           roi06, ...
                                           secondsInsideStartBoxAfterForcedRun, ...
                                           optionWhiteNoiseInsideStartBoxAfterForcedRun, ...  
                                           optionUseStartGateAfterFreeingTheMouse, ...
                                           optionLickToOpenValve, ...
                                           laserSelection, ...
                                           pathFlags, ...
                                           currentInhibitionStatus); %#ok<*ASGLU>
      
  end
  
  % Just increment the counter
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
  
  % Update trial times
  sequenceOfAlternation = wxupdatetrialtimes(sequenceOfAlternation, ...
                                             counterOfRuns, ...
                                             startTimeForcedRun, ...
                                             endTimeForcedRun, ...
                                             startTimeForcedRun, ...
                                             endTimeForcedRun);
                                           
  % Clear current trial times
  clear startTimeForcedRun
  clear endTimeForcedRun                                           

  % Update figure with alternation sides and trial times
  wxvisualizealternationandtrialtimes(sequenceOfAlternation, counterOfRuns, maxNumberOfRunsShapingPeriod)                                           
                                           
  % Search for "smooth stop" flag
  % If the flag is found, break the loop and delete the flag
  flagSmoothStop = exist([pathSmoothStop '/' 'flagsmoothstop.mat'],'file');  
  if (flagSmoothStop == 2)
    counterOfRuns = maxNumberOfRunsShapingPeriod;
    delete([pathSmoothStop '/' 'flagsmoothstop.mat'])
  end
  
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
  
  % Get current side
  if sequenceOfSides{counterOfRuns,1} == 'R'
    currentSide = 'right';
  elseif sequenceOfSides{counterOfRuns,1} == 'L'
    currentSide = 'left';
  end
  
  % Get current inhibition status (laser on vs laser off)
  if sequenceOfSides{counterOfRuns,14} == 1
    currentInhibitionStatus = 'on';
  elseif sequenceOfSides{counterOfRuns,14} == 0
    currentInhibitionStatus = 'off';
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
                                        roi06, ...
                                        secondsInsideStartBoxAfterForcedRun, ...
                                        secondsInsideStartBoxAfterFreeRun, ...
                                        optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
                                        optionWhiteNoiseInsideStartBoxAfterFreeRun, ...  
                                        optionUseStartGateAfterFreeingTheMouse, ...
                                        optionLickToOpenValve, ...
                                        laserSelection, ...
                                        pathFlags, ...
                                        currentInhibitionStatus);
                                                              
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
  
  % Search for "smooth stop" flag
  % If the flag is found, break the loop and delete the flag
  flagSmoothStop = exist([pathSmoothStop '/' 'flagsmoothstop.mat'],'file');  
  if (flagSmoothStop == 2)
    counterOfRuns = maxNumberOfRunsLearningPeriod;
    delete([pathSmoothStop '/' 'flagsmoothstop.mat'])
  end
                                       
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