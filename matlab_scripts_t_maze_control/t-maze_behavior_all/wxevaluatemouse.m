%% wxevaluatemouse
 % Evaluates whether the mouse deserves a reward
 % JL Alatorre-Warren

function [daqDevice, ...
          daqSession, ...
          mainLogbook] = wxevaluatemouse(currentSide, ...
                                         detectedSide, ...
                                         mainLogbook, ...
                                         daqDevice, ...
                                         daqSession, ...
                                         optionLickToOpenValve)

% To get the reward, the mouse has to choose the opposite side
if strcmp(currentSide, 'right') == 1
  rewardSide = 'left';
elseif strcmp(currentSide, 'left') == 1
  rewardSide = 'right';
end  

% Evaluate and reward mouse
if strcmp(detectedSide,currentSide) == 1
  
  mouseStatus = 'stupidMouse';
  disp(mouseStatus)
  mainLogbook = wxservo(rewardSide,'close',mainLogbook);

elseif strcmp(detectedSide,currentSide) == 0

  mouseStatus = 'cleverMouse';
  disp(mouseStatus)
  
  % Close the opposite gate (currentSide)
  mainLogbook = wxservo(currentSide,'close',mainLogbook);
    
  % If optionLickToOpenValve is 1, the mouse opens the valve by licking it
  if optionLickToOpenValve == 1
    mainLogbook = wxlickingdetection(daqSession,rewardSide,mainLogbook);
  end	
  
  % Open valve
  [daqDevice, daqSession, mainLogbook] = wxopenvalve(rewardSide, ...
                                                   daqDevice, ...
                                                   daqSession, ...
                                                   mainLogbook);
  
  % Close 
  mainLogbook = wxlickingdetection(daqSession,rewardSide,mainLogbook);

end