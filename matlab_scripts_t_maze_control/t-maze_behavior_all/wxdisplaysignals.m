%% wxdisplaysignals
 % Display Ca 2+ and licking signals using the National Instruments USB-6218 device
 % JL Alatorre-Warren
 
%% Initial setup

% Get relevant information from wxtmaze
[pathRoot,...
 pathData,...
 pathRois,...
 pathTrials,...
 pathTrialList,...
 structureOfTrials,...
 listOfTrials,...
 trialName,...
 pathCurrentTrial,...
 mouseName,...
 experimentName,...
 timestampVector,...
 dateString,...
 timeString] = wxlistentowxtmaze; %#ok<*ASGLU>

% Get the rest of the information
[~,...
 ~,...
 ~,...
 ~,...
 ~,...
 ~,...
 pathMouseName,...
 pathDateString,...
 pathExperimentName,...
 pathTimeString,...
 pathFlags,...
 pathSignals,...
 pathFrames] = wxcreatedirectories(timestampVector,... 
                                   mouseName,...
                                   dateString,...
                                   experimentName,...
                                   timeString);
                                  
% Send flag to wxtmaze
flagDisplaySignalsReady = true; %#ok<*NASGU>
save([pathFlags '/' 'flagdisplaysignalsready.mat'], 'flagDisplaySignalsReady');

%% Display signals
%  National Instruments DAQ (Multifunction Data Acquisition)
%  National Instruments USB-6218

% Set basic parameters
daqSessionRate               = 2000;
samplesPerDaqSession         = 200;
numberOfSecondsToBeDisplayed = 30;
frequencyRangeLaser490nm     = [489.5 490.5];
frequencyRangeLaser573nm     = [572.5 573.5];

% Compute secondary parameters
numberOfModulatedSamplesToBeDisplayed   = numberOfSecondsToBeDisplayed * daqSessionRate;
numberOfDemodulatedSamplesToBeDisplayed = numberOfSecondsToBeDisplayed*(daqSessionRate/samplesPerDaqSession);
recentDaqSignals                        = zeros(numberOfModulatedSamplesToBeDisplayed,3);
recentDemodulatedCaSignalsLaser490nm    = zeros(numberOfDemodulatedSamplesToBeDisplayed,1);
recentDemodulatedCaSignalsLaser573nm    = zeros(numberOfDemodulatedSamplesToBeDisplayed,1);

% Get a list with all currently saved DAQ signal files
structureOfDaqSignals          = dir(pathSignals);
listOfDaqSignals               = ({structureOfDaqSignals(3:end).name})';
currentNumberOfSavedDaqSignals = length(listOfDaqSignals);

% Prepare figure parameters
figurePlots = figure;
figurePlots.OuterPosition = [2150 41 411 1400];

ii = 1;
maxPeakToPeak = 0;
maxsLaser490nm = 0;
minsLaser490nm = 0;
upperLimitLaser490nm = 0;
lowerLimitLaser490nm = 0;
maxsLaser573nm = 0;
minsLaser573nm = 0;
upperLimitLaser573nm = 0;
lowerLimitLaser573nm = 0;


% Display signals until flag sent from wxtmaze is found
flagFound = false;
while flagFound == false
  
  % Track time
  tic
    
  % Check whether the MATLAB file daqSignalsMostRecent.mat exists
  if (exist([pathSignals '/' 'daqSignalsMostRecent.mat'], 'file')) == 2
    
    % If the file daqSignalsMostRecent.mat is being updated by
    % wxrecordsingals, an error exception is permitted and the file is load
    % again
    try
      load([pathSignals '/' 'daqSignalsMostRecent.mat'], 'daqSignals');
    catch
      warning('The file daqSignalsMostRecent.mat is not currently available. Trying again.');
    end
    
  end
  
  % Check whether the MATLAB variable daqSignals exists
  % In other words, check whether the file daqSignalsMostRecent.mat was successfully loaded
  if (exist('daqSignals', 'var')) == 1
    
    % Update ii if and only if the variable daqSignals exists
    ii = ii + 1;
    
    % DAQ signals to be displayed
    recentDaqSignals((end-samplesPerDaqSession+1):end,:) = daqSignals;
    
    % Demodulate Ca 2+ signals
    tic, demodulatedCurrentCaSignalLaser490nm = bandpower(daqSignals(:,1),daqSessionRate,frequencyRangeLaser490nm); toc;
    tic, demodulatedCurrentCaSignalLaser573nm = bandpower(daqSignals(:,1),daqSessionRate,frequencyRangeLaser573nm); toc;

    % History of demodulated Ca 2+ signals
    recentDemodulatedCaSignalsLaser490nm(end) = demodulatedCurrentCaSignalLaser490nm;
    recentDemodulatedCaSignalsLaser573nm(end) = demodulatedCurrentCaSignalLaser573nm;
    
    % Dynamically resize the vertical axis of the plot of the demodulated
    % signals recentDemodulatedCaSignalsLaser490nm and recentDemodulatedCaSignalsLaser573nm
    currentmaxLaser490nm = max(recentDemodulatedCaSignalsLaser490nm);
    currentminLaser490nm = min(recentDemodulatedCaSignalsLaser490nm);
    currentmaxLaser573nm = max(recentDemodulatedCaSignalsLaser573nm);
    currentminLaser573nm = min(recentDemodulatedCaSignalsLaser573nm);
    if (currentmaxLaser490nm >= maxsLaser490nm) || (currentminLaser490nm <= minsLaser490nm)
      upperLimitLaser490nm = max(recentDemodulatedCaSignalsLaser490nm) + 0.10*max(recentDemodulatedCaSignalsLaser490nm);
      lowerLimitLaser490nm = min(recentDemodulatedCaSignalsLaser490nm) - 0.10*(abs(min(recentDemodulatedCaSignalsLaser490nm)));
    end
    if (currentmaxLaser573nm >= maxsLaser573nm) || (currentminLaser573nm <= minsLaser573nm)
      upperLimitLaser573nm = max(recentDemodulatedCaSignalsLaser573nm) + 0.10*max(recentDemodulatedCaSignalsLaser573nm);
      lowerLimitLaser573nm = min(recentDemodulatedCaSignalsLaser573nm) - 0.10*(abs(min(recentDemodulatedCaSignalsLaser573nm)));    
    end
    
    subplot(5,1,1)
    plot(recentDemodulatedCaSignalsLaser490nm)
    axis([1 numberOfDemodulatedSamplesToBeDisplayed lowerLimitLaser490nm upperLimitLaser490nm])
    title({'Ca 2+ Demodulated (490 Hz)',''})
    
    subplot(5,1,2)
    plot(recentDemodulatedCaSignalsLaser573nm)
    title({'Ca 2+ Demodulated (573 Hz)',''})
    axis([1 numberOfDemodulatedSamplesToBeDisplayed lowerLimitLaser573nm upperLimitLaser573nm])
    
    subplot(5,1,3)
    plot(recentDaqSignals(:,1))
    axis([1 numberOfModulatedSamplesToBeDisplayed -2.000 1.000])
    title('Ca 2+ Modulated')
    
    subplot(5,1,4)
    plot(recentDaqSignals(:,2))
    axis([1 numberOfModulatedSamplesToBeDisplayed -1.000 10.000])
    title('Licking Detector - Left Arm')
    
    subplot(5,1,5)
    plot(recentDaqSignals(:,3))
    axis([1 numberOfModulatedSamplesToBeDisplayed -1.000 10.000])
    title('Licking Detector - Right Arm')
    
    % Pause allows the plots to be properly displayed
    pause(0.030)
    
    % Update recentDaqSignals and recentDemodulatedCaSignalsLaser490nm by shifting
    % their values (remove oldest values to create space for new values)
    recentDaqSignals(1:(end-samplesPerDaqSession),:) = recentDaqSignals((1+samplesPerDaqSession):(end),:);
    recentDemodulatedCaSignalsLaser490nm(1:(end-1)) = recentDemodulatedCaSignalsLaser490nm(2:end);
    recentDemodulatedCaSignalsLaser573nm(1:(end-1)) = recentDemodulatedCaSignalsLaser573nm(2:end);
    
    % Stop stopwatch
    toc
    
    % Search flag and break the loop if found
    flagSequenceCompleted = exist([pathFlags '/' 'flagsequencecompleted.mat'],'file');  
    if (flagSequenceCompleted == 2)
      flagFound = true;
    end    
    
  end
  
end

% Close MATLAB session
exit