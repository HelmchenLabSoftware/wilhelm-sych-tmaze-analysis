%% wxlickingdetection
 % Detects licking in valve and opens start gate
 % JL Alatorre-Warren
 
function mainLogbook = wxlickingdetection(daqSession,currentSide,mainLogbook)
 
% Assign the appropriate channel from the NI-6218 device
switch currentSide
  case 'right'
    currentSideChannel = 2;
  case 'left'
    currentSideChannel = 1;
end

% Detect mouse licking
ii = 0;
endLoopCondition = 10000;
while ii < endLoopCondition

  % Acquire data
  % The maximum datum from the data acquired will be used to break the loop
  [lickingVoltage,~] = daqSession.startForeground;
  lickingVoltageMax = max(lickingVoltage(:,currentSideChannel));

  % If a licking is detected, create a timestamp and exit the loop
  if lickingVoltageMax > 0.100
    ii = endLoopCondition;
    timestampVector = clock;
    mainLogbook = wxupdatelogbook(mainLogbook, ['licking' '_' currentSide], timestampVector);
  end
  
	% Display detection score and ROI name
  disp(['Licking voltage (' currentSide '): ' num2str(lickingVoltageMax)])

  ii = ii+1;
end