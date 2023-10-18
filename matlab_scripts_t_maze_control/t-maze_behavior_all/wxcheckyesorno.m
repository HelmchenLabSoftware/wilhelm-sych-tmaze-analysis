%% wxcheckyesorno
 % Checks whether the user correctly assigned a yes or no answer
 % JL Alatorre-Warren
 
function [correctInput, yesOrNo] = wxcheckyesorno(userInput)
 
switch(userInput)
case 'yes'
  yesOrNo = 1;
  correctInput = true;
case 'y'
  yesOrNo = 1;
  correctInput = true;
case 'Yes'
  yesOrNo = 1;
  correctInput = true;
case 'Y'
  yesOrNo = 1;
  correctInput = true;
case 'YES'
  yesOrNo = 1;
  correctInput = true;
case '1'
  yesOrNo = 1;
  correctInput = true;
case 'N'
  yesOrNo = 0;
  correctInput = true;
case 'n'
  yesOrNo = 0;
  correctInput = true;
case 'no'
  yesOrNo = 0;
  correctInput = true;
case 'No'
  yesOrNo = 0;
  correctInput = true;
case 'NO'
  yesOrNo = 0;
  correctInput = true;
case '0'
  yesOrNo = 0;
  correctInput = true;
case 'nope'
  yesOrNo = 0;
  correctInput = true;
case 'Nope'
  yesOrNo = 0;
  correctInput = true;
case 'NOPE'
  yesOrNo = 0;
  correctInput = true;
otherwise
  errorMessage = [newline 'Sorry, not a valid answer. ' newline];
  disp(errorMessage);
end
