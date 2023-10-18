%% wxuserinput
 % Request input from the user
 % JL Alatorre-Warren
 
function [mouseName, ...
          experimentName, ...
          maxNumberOfRunsShapingPeriod, ...
          maxNumberOfRunsLearningPeriod, ...
          secondsInsideStartBoxAfterForcedRun, ...
          secondsInsideStartBoxAfterFreeRun, ...
          optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
          optionWhiteNoiseInsideStartBoxAfterFreeRun, ...
          optionUseStartGateAfterFreeingTheMouse, ...
          optionLickToOpenValve] = wxuserinput

% User input: mouse name
prompt = 'Provide the name of the mouse: ';
mouseName = input(prompt, 's');

% User input: experiment name
prompt = 'Provide the name of the experiment: ';
experimentName = input(prompt, 's');

% Maximum number of full runs during the shaping period
correctName = false;
while correctName == false
  prompt = 'Provide the number of runs for the shaping period: ';
  userInput = input(prompt, 's');
  userInput = str2double(userInput);
	if (numel(userInput)==1) && (isnan(userInput)==0)
    maxNumberOfRunsShapingPeriod = uint64(userInput);
    correctName = true;
    clear userInput
  else
    errorMessage = [newline 'Sorry.' ...
                    newline 'That''s not a valid number.' ...
                    newline 'Please try again.' ...
                    newline];
    disp(errorMessage);
  end
end

% Maximum number of full runs during the learning period
correctName = false;
while correctName == false
  prompt = 'Provide the number of runs for the learning period: ';
  userInput = input(prompt, 's');
  userInput = str2double(userInput);
	if (numel(userInput)==1) && (isnan(userInput)==0)
    maxNumberOfRunsLearningPeriod = uint64(userInput);
    correctName = true;
    clear userInput
  else
    errorMessage = [newline 'Sorry.' ...
                    newline 'That''s not a valid number.' ...
                    newline 'Please try again.' ...
                    newline];
    disp(errorMessage);
  end
end

% Waiting time in seconds inside start box after forced run
correctName = false;
while correctName == false
  prompt = 'Provide the waiting time in seconds inside the start box after the forced run: ';
  userInput = input(prompt, 's');
  userInput = str2double(userInput);
	if (numel(userInput)==1) && (isnan(userInput)==0)
    secondsInsideStartBoxAfterForcedRun = userInput;
    correctName = true;
    clear userInput
  else
    errorMessage = [newline 'Sorry.' ...
                    newline 'That''s not a valid number.' ...
                    newline 'Please try again.' ...
                    newline];
    disp(errorMessage);
  end
end

% Waiting time in seconds inside start box after free (open) run
correctName = false;
while correctName == false
  prompt = 'Provide the waiting time in seconds inside the start box after the free run: ';
  userInput = input(prompt, 's');
  userInput = str2double(userInput);
	if (numel(userInput)==1) && (isnan(userInput)==0)
    secondsInsideStartBoxAfterFreeRun = userInput;
    correctName = true;
    clear userInput
  else
    errorMessage = [newline 'Sorry.' ...
                    newline 'That''s not a valid number.' ...
                    newline 'Please try again.' ...
                    newline];
    disp(errorMessage);
  end
end

% User input: option white noise after forced run
correctInput = false;
while correctInput == false
  prompt = 'Do you want to reproduce white noise inside the start box after the forced run : ';
  userInput = input(prompt, 's');
  [correctInput, optionWhiteNoiseInsideStartBoxAfterForcedRun] = wxcheckyesorno(userInput);
end

% User input: option white noise after free run
correctInput = false;
while correctInput == false
  prompt = 'Do you want to reproduce white noise inside the start box after the free run : ';
  userInput = input(prompt, 's');
  [correctInput, optionWhiteNoiseInsideStartBoxAfterFreeRun] = wxcheckyesorno(userInput);
end

% User input: option ignore the start gate after freeing the mouse
correctInput = false;
while correctInput == false
  prompt = 'Do you want to use the start gate after freeing the mouse from the start box: ';
  userInput = input(prompt, 's');
  [correctInput, optionUseStartGateAfterFreeingTheMouse] = wxcheckyesorno(userInput);
end

% User input: option lick to open valve (optionLickToOpenValve)
correctInput = false;
while correctInput == false
  prompt = 'Do you want to use the LickToOpenValve option: ';
  userInput = input(prompt, 's');
  [correctInput, optionLickToOpenValve] = wxcheckyesorno(userInput);
end