%% wxnibox6210
 % Creates and starts a DAQ session using the NI USB-6210 device
 % JL Alatorre-Warren
 
function [daqDevice, daqSession] = wxnibox6210

% Discover Available Devices 
% Use the |daq.getDevices| command to display a list of devices available
% to your machine and MATLAB(R).
daqDevice = daq.getDevices;

% Display relevant information
% daqDevice(1,1) corresponds to NI USB-6218
% daqDevice(1,2) corresponds to NI USB-6210
% Here, we are only interested in NI USB-6210
display(daqDevice(1,2))
display(daqDevice(1,2).Terminals)
display(daqDevice(1,2).Vendor)
display(daqDevice(1,2).ID)
display(daqDevice(1,2).Model)
display(daqDevice(1,2).Subsystems)
display(daqDevice(1,2).Description)

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

% Add analog input channels
% The channel AI1 is used for the left licking detector
% The channel AI2 is used for the right licking detector
[channelAnalogInputForLeftLickingDetector ,...
 indexAnalogInputForLeftLickingDetector]  = addAnalogInputChannel(daqSession,daqDevice(1,2).ID, 1,'Voltage');
[channelAnalogInputForRightLickingDetector,...
 indexAnalogInputForRightLickingDetector] = addAnalogInputChannel(daqSession,daqDevice(1,2).ID, 2,'Voltage');

% Display channels and indices
display(channelAnalogInputForLeftLickingDetector)
display(channelAnalogInputForRightLickingDetector)
display(indexAnalogInputForLeftLickingDetector)
display(indexAnalogInputForRightLickingDetector)