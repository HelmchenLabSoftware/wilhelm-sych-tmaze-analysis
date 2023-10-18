%% wxcheckvalves
 % Reminds the user to fill the T-maze with water and launches the T-maze
 % when the user is ready to do so
 % JL Alatorre-Warren
 
function [daqDevice,daqSession] = wxcheckvalves(daqDevice,daqSession)

% User input: check the water containers
correctInput = false;
while correctInput == false
  prompt = 'Do you want to check the water containers: ';
  userInput = input(prompt, 's');
  switch(userInput)
  case 'yes'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case 'y'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case 'Yes'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case 'Y'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case 'YES'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case '1'
    [daqDevice,daqSession] = wxloadwater(daqDevice,daqSession);
  case 'N'
    correctInput = true;
  case 'n'
    correctInput = true;
  case 'no'
    correctInput = true;
  case 'No'
    correctInput = true;
  case 'NO'
    correctInput = true;
  case '0'
    correctInput = true;
  otherwise
    errorMessage = [newline 'Sorry, not a valid answer. ' newline];
    disp(errorMessage);
  end
end