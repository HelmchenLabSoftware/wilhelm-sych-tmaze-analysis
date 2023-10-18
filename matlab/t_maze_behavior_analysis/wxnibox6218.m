%% wxnibox6218
 % Creates and starts a DAQ session using the NI USB-6218 device
 % JL Alatorre-Warren
 
function [daqDevice, daqSession] = wxnibox6218

% Discover Available Devices 
% Use the |daq.getDevices| command to display a list of devices available
% to your machine and MATLAB(R).
daqDevice = daq.getDevices;

% Display relevant information
% daqDevice(1,1) corresponds to NI USB-6218
% daqDevice(1,2) corresponds to NI USB-6008
% Here, we are only interested in NI USB-6218
display(daqDevice(1,1))
display(daqDevice(1,1).Terminals)
display(daqDevice(1,1).Vendor)
display(daqDevice(1,1).ID)
display(daqDevice(1,1).Model)
display(daqDevice(1,1).Subsystems)
display(daqDevice(1,1).Description)

% Create a DAQ session
% The |daq.createSession| command creates a session. The session contains
% information describing the hardware, scan rate, duration, and other
% properties associated with the acquisition.
% The duration in seconds is set to 0.050 
% The duration in seconds in LabVIEW was 0.035
daqSession = daq.createSession('ni');
daqSession.DurationInSeconds = 0.050;
daqSession.Rate = 2000;
display(daqSession.Channels)

% Add an analog input channel from this device to the session
% The channel AI0 is used for the Ca 2+ signal
% The channel AI1 is used for the left licking detector
% The channel AI2 is used for the right licking detector
[channelAnalogInputForCalciumSignal,...
 indexAnalogInputForCalciumSignal]       = addAnalogInputChannel(daqSession,daqDevice(1,1).ID,13,'Voltage');
[channelAnalogInputForLeftLickingDetector,...
 indexAnalogInputForLeftLickingDetector] = addAnalogInputChannel(daqSession,daqDevice(1,1).ID, 1,'Voltage');
[channelAnalogInputForRightLickingDetector,...
 indexAnalogInputForRighLickingDetector] = addAnalogInputChannel(daqSession,daqDevice(1,1).ID, 2,'Voltage');

% Display channels and indices
display(channelAnalogInputForCalciumSignal)
display(channelAnalogInputForLeftLickingDetector)
display(channelAnalogInputForRightLickingDetector)
display(indexAnalogInputForCalciumSignal)
display(indexAnalogInputForLeftLickingDetector)
display(indexAnalogInputForRighLickingDetector)