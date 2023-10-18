%% wxlasermodulation
%  Switch on/off laser by using digital modulation using the National Instruments USB-6218 device
%  Maria Chernysheva

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

%% Search for "flagLaserOn" and "flagLaserOff" flags

% Prepare DAQ device and exit loop flag
daqDevice = daq.getDevices;
daqSession = daq.createSession('ni');
flagEndOfExperiment = false;

% Make sure that the inhibition laser is off before starting
inhibitionLaserStatus = 'Inhibition Laser Off';
[channelDigitalOutputForLaser, indexDigitalOutputForLaser] = addDigitalChannel(daqSession, daqDevice(1,1).ID, 'port1/line0', 'OutputOnly');
outputSingleScan(daqSession,0);
removeChannel(daqSession,indexDigitalOutputForLaser)

% %first let's switch off the laser:
% 
% pathLaserOFF = 'I:\2017_m001-m006_backup/LaserOFF';
%   flagLaserOFF = true; %#ok<NASGU>
%   save([pathLaserOn '/' 'flagLaserOFF.mat'], 'flagLaserOn');
% 
%  % Turn inhibition laser off
%   flagLaserOff = exist([pathLaserOn '/' 'flagLaserOFF.mat'],'file');  
%   if (flagLaserOff == 2)
%     [channelDigitalOutputForLaser, indexDigitalOutputForLaser] = addDigitalChannel(daqSession, daqDevice(1,1).ID, 'port1/line0', 'OutputOnly');
%     outputSingleScan(daqSession,0);
%     removeChannel(daqSession,indexDigitalOutputForLaser)
%     delete([pathLaserOn '/' 'flagLaserOFF.mat']) 
%     inhibitionLaserStatus = 'End trigger for the Miniscope';
%   end
  
%start the recording of the miniscope data (send the trigger, by creating the flag)
pathLaserOn = 'I:\2017_m001-m006_backup/LaserOn';
  flagLaserOn = true; %#ok<NASGU>
  save([pathLaserOn '/' 'flagLaserOn.mat'], 'flagLaserOn');
% Track time
tic

% Search for flags
% If the flag is found, break the loop and delete the flag
while flagEndOfExperiment == false
  
  % Display current laser status
  %disp([inhibitionLaserStatus ' :: ' num2str(toc)])
    
  % Turn inhibition laser on
  flagLaserOn = exist([pathLaserOn '/' 'flagLaserOn.mat'],'file'); 

  if (flagLaserOn == 2)
    [channelDigitalOutputForLaser, indexDigitalOutputForLaser] = addDigitalChannel(daqSession, daqDevice(1,1).ID, 'port1/line0', 'OutputOnly');
    outputSingleScan(daqSession,1);
    removeChannel(daqSession,indexDigitalOutputForLaser)
    delete([pathFlags '/' 'flagLaserOn.mat'])
    inhibitionLaserStatus = 'Miniscope trigger';
  end
  
 
  
  % Search for end of experiment flag
  flagSequenceCompleted = exist([pathFlags '/' 'flagsequencecompleted.mat'],'file');
  if (flagSequenceCompleted == 2)
    flagEndOfExperiment = true;
  end
  
  % Just wait a bit
  pause(0.025)
  
end

pathLaserOFF = 'I:\2017_m001-m006_backup/LaserOFF';
  flagLaserOFF = true; %#ok<NASGU>
  save([pathLaserOn '/' 'flagLaserOFF.mat'], 'flagLaserOn');

 % Turn inhibition laser off
  flagLaserOff = exist([pathLaserOn '/' 'flagLaserOFF.mat'],'file');  
  if (flagLaserOff == 2)
    [channelDigitalOutputForLaser, indexDigitalOutputForLaser] = addDigitalChannel(daqSession, daqDevice(1,1).ID, 'port1/line0', 'OutputOnly');
    outputSingleScan(daqSession,0);
    removeChannel(daqSession,indexDigitalOutputForLaser)
    delete([pathLaserOn '/' 'flagLaserOFF.mat']) 
    inhibitionLaserStatus = 'End trigger for the Miniscope';
  end
% Close MATLAB session
exit