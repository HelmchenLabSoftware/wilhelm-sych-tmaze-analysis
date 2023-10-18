%% wxloadwater
 % Checks that the valves are delivering water and allows the user to fill
 % the water containers, if needed.
 % JL Alatorre-Warren
 
function [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession)

% User input: number of pulses
correctInput = false;
while correctInput == false
  prompt = 'Provide the number pulses: ';
  userInput = input(prompt, 's');
  userInput = str2double(userInput);
	if (numel(userInput)==1) && (isnan(userInput)==0)
    numberOfPulses = uint64(userInput);
    correctInput = true;
    clear userInput
  else
    errorMessage = [newline 'Sorry.' ...
                    newline 'That''s not a valid number.' ...
                    newline 'Please try again.' ...
                    newline];
    disp(errorMessage);
  end
end

% Right valve
% Add a channel that generate pulses on a counter/timer subsystem
[channelDigitalOutputForRightValve,...
 indexDigitalOutputForRightValve] = addCounterOutputChannel(daqSession,daqDevice(1,2).ID,1,'PulseGeneration');
display(channelDigitalOutputForRightValve)
display(indexDigitalOutputForRightValve)
daqSession.DurationInSeconds = 0.050;
% Right valve
% The default frequency is 100.
% Change it to 1 so that there is a single continuous pulse in the
% daqSession, so that the valves can open.
channelDigitalOutputForRightValve.Frequency = 1;
display(channelDigitalOutputForRightValve.Frequency)

% Open the right valve n times
for ii = 1:numberOfPulses
  daqSession.startForeground;
end

% Remove channel associated with the right valve
removeChannel(daqSession,indexDigitalOutputForRightValve)

% Left valve
% Add a channel that generate pulses on a counter/timer subsystem
[channelDigitalOutputForLeftValve,...
 indexDigitalOutputForLeftValve] = addCounterOutputChannel(daqSession,daqDevice(1,2).ID,0,'PulseGeneration');
display(channelDigitalOutputForLeftValve)
display(indexDigitalOutputForLeftValve)

% The left valve seems to deliver less water than the right one
% That's why we set this valve to be open for a longer period of time
daqSession.DurationInSeconds = 0.050;

% Left valve
% The default frequency is 100.
% Change it to 1 so that there is a single continuous pulse in the
% daqSession, so that the valves can open.
channelDigitalOutputForLeftValve.Frequency = 1;
display(channelDigitalOutputForLeftValve.Frequency)

% Open left valve n times
for ii = 1:numberOfPulses
  daqSession.startForeground;
end

% Remove channel associated with the left valve
removeChannel(daqSession,indexDigitalOutputForLeftValve)

% Set the duration of the DAQ session to its default value
daqSession.DurationInSeconds = 0.050;