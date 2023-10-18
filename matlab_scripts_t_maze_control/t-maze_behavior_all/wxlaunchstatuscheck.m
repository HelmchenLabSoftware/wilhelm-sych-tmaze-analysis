%% wxlaunchstatuscheck
 % Launches the T-maze when the user is ready to do so
 % JL Alatorre-Warren
 
% User input: launch status check
correctInput = false;
while correctInput == false
  prompt = 'Do you want to start the experiment: ';
  userInput = input(prompt, 's');
  switch(userInput)
  case 'yes'
    correctInput = true;
  case 'y'
    correctInput = true;
  case 'Yes'
    correctInput = true;
  case 'Y'
    correctInput = true;
  case 'YES'
    correctInput = true;
  case '1'
    correctInput = true;
  otherwise
    errorMessage = [newline 'Sorry, not a valid answer. ' newline];
    disp(errorMessage);
  end
end