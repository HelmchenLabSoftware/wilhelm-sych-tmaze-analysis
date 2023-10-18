%% wxopenvalve
 % Confines the mouse in the start box and opens the corresponding valve
 % and the corresponding gates
 % JL Alatorre-Warren
 
function [daqDevice, ...
          daqSession, ...
          mainLogbook] = wxopenvalve(currentSide, ...
                                     daqDevice, ...
                                     daqSession, ...
                                     mainLogbook)
 
switch currentSide

  case 'right'

    % Right valve
    % Add a channel that generate pulses on a counter/timer subsystem
    [channelDigitalOutputForRightValve,...
     indexDigitalOutputForRightValve] = addCounterOutputChannel(daqSession,daqDevice(1,2).ID,1,'PulseGeneration'); %#ok<ASGLU>

 daqSession.DurationInSeconds = 0.040;
    % Important:
    % The default frequency is 100.
    % Change it to 1 so that there is a single continuous pulse in the
    % daqSession, so that the valves can open.
    channelDigitalOutputForRightValve.Frequency = 1;


  case 'left'

    % Left valve
    % Add a channel that generate pulses on a counter/timer subsystem
    [channelDigitalOutputForLeftValve,...
     indexDigitalOutputForLeftValve] = addCounterOutputChannel(daqSession,daqDevice(1,2).ID,0,'PulseGeneration');
    
    % The left valve is opened for a little bit shorter time than the right
    % one, as more water seems to be delivered on the left arm of the T-maze
    daqSession.DurationInSeconds = 0.040;

    % Important:
    % The default frequency is 100.
    % Change it to 1 so that there is a single continuous pulse in the
    % daqSession, so that the valves can open.
    channelDigitalOutputForLeftValve.Frequency = 1;

end

% Open valve & timestamp
daqSession.startForeground;
timestampVector = clock;
mainLogbook = wxupdatelogbook(mainLogbook, ['valve_opens' '_' currentSide], timestampVector);

switch currentSide

  case 'right'

    % Remove channel (right valve)
    removeChannel(daqSession,indexDigitalOutputForRightValve)

  case 'left'

    % Remove channel (left valve)
    removeChannel(daqSession,indexDigitalOutputForLeftValve)
    
    % Change the duration of the session to the default value
    daqSession.DurationInSeconds = 0.040;

end