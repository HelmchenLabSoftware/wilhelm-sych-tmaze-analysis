%% wxdelayperiod
% Keeps the mouse inside the start box and adds white noise
% JL Alatorre-Warren
 
function wxdelayperiod(secondsInsideStartBox,optionWhiteNoiseInsideStartBox)

% Create (uniformily distributed) white noise
% Play it in the speakers
if optionWhiteNoiseInsideStartBox == 1
  whiteNoise = rand(1000000,1);
  sound(whiteNoise)
end

% Waiting time in seconds inside the start box
pause(secondsInsideStartBox);

% Stop sound
clear sound