%% wxservo
 % Function used to move the T-maze's servomotor gates and to check if the
 % servomotor actually moved
 % JL Alatorre Warren

function mainLogbook = wxservo(nameServo, gateState, mainLogbook)

ii = 1;
while ii < 10
  
  % Display the attempt number
  disp(['Servo iteration: ' num2str(ii)]);

  % Workhorse function to move the servomotor gates
  [mainLogbook, successServo] = wxservoworkhorse(nameServo, ...
                                                 gateState, ...
                                                 mainLogbook);

  if successServo == true
    ii = 11;
  else
    ii = ii+1;
  end
  
end