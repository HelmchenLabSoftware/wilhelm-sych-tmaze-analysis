%% wxsetthemousefree
 % Opens the start box and the gate that leads to the opened valve (reward)
 % JL Alatorre-Warren

function mainLogbook = wxsetthemousefree(currentSide, mainLogbook)
 
switch currentSide

  case 'right'

    % Move right servo
    mainLogbook = wxservo('right','open', mainLogbook);
    mainLogbook = wxservo('start','open', mainLogbook);

  case 'left'

    % Move left and start servos
    mainLogbook = wxservo('left','open', mainLogbook);
    mainLogbook = wxservo('start','open', mainLogbook);

end