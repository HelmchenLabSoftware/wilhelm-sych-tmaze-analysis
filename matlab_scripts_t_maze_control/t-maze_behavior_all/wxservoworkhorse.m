%% wxservoworkhorse
 % Workhorse function for the manipulation of the T-maze's servomotor gates
 % JL Alatorre Warren

function [mainLogbook, successServo] = wxservoworkhorse(nameServo, gateState, mainLogbook)

% This loads the phidget21 library according to the OS
loadphidget21;  

% The function libpointer creates a pointer object
% for use with shared libraries
handle = libpointer('int32Ptr');

% The function calllib calls the Phidgets functions (shared library)
% phidget21 is the library name
% CPhidgetAdvancedServo_create and CPhidget_open are functions in phidget21
calllib('phidget21', 'CPhidgetAdvancedServo_create', handle);
calllib('phidget21', 'CPhidget_open', handle, -1);

% The function libpointer creates a pointer object
% for use with shared libraries
valPtr = libpointer('doublePtr', 0);

% Assigns servomotor index to 0 if the user chooses an incorrect one
% Indexes: start is 0, left is 2, and right is 4
switch(nameServo)
  case 'start'
    indexServo = 0;
  case 'right'
    indexServo = 4;
  case 'left'
    indexServo = 2;
  otherwise
    indexServo = 0;
end

% Assigns the correct servomotor position according to the specific
% servomotor and the gate status (open or close)
% Indexes: start is 0, left is 2, and right is 4
switch(gateState)
  case 'close'
    if indexServo == 0
      posn = 110;
    elseif indexServo == 4
      posn = 140;
    elseif indexServo == 2
      posn = 105;
    end
  case 'open'
    if indexServo == 0
      posn = 57;
    elseif indexServo == 4
      posn = 195;
    elseif indexServo == 2
      posn =40;
    end
end

% attemptNumber = 1;
% successServo = false;
% while successServo == false
% disp(['Number of attempts: ' num2str(attemptNumber)])

% CPhidget_waitForAttachment waits for attachment to happen.
% This can be called wirght after calling CPhidget_open,
% as an alternative to using the attach handler.
if calllib('phidget21', 'CPhidget_waitForAttachment', handle, 2500) == 0

  disp('Opened Advanced Servo')

  % The variable t is used to control the timing of the opening and 
  % closing of the gates
  t = timer('TimerFcn','disp(''moving servo...'')', 'StartDelay', 0.300);

  % Get and maximum velocity
  calllib('phidget21', 'CPhidgetAdvancedServo_getVelocityMax', handle, indexServo, valPtr);
  maxVelocity = get(valPtr, 'Value');
  
  % Get and maximum accelaration
  calllib('phidget21', 'CPhidgetAdvancedServo_getAccelerationMax', handle, indexServo, valPtr);
  maxAccel = get(valPtr, 'Value');
  
  % Set maximum velocity and maximum accelaration
  calllib('phidget21', 'CPhidgetAdvancedServo_setAcceleration',  handle, 0, maxAccel);
  calllib('phidget21', 'CPhidgetAdvancedServo_setVelocityLimit', handle, 0, maxVelocity);

  % Power the servomotor on
  calllib('phidget21', 'CPhidgetAdvancedServo_setEngaged', handle, indexServo, 1);

  start(t);
  wait(t);
  calllib('phidget21', 'CPhidgetAdvancedServo_setPosition', handle, indexServo, posn);

% % Get current position
% % If the current position is correct, break the loop and create timestamp
% calllib('phidget21', 'CPhidgetAdvancedServo_getPosition', handle, indexServo, valPtr);
% checkPosn = valPtr.Value;
% if checkPosn == posn
%   %successServo = true;
%   disp(['Desired position: ' num2str(posn)])
%   disp(['Current position: ' num2str(checkPosn)])
% end

    successServo = true;

else
  
    % Subtle servomotor error notification
    disp('Could not open advanced servo')
    
    % Let the servo take a breath (necessary to prevent MATLAB to crash)
    pause(1.000);
    
    % Dramatic servomotor error notification
    wxnotifyservoerror;
    
    % Request user to press any key to continue
    pause;
    
    % Pass this flag as a function output to notify MATLAB to try to move the servomotor again
    successServo = false;
    
end

timestampVector = clock;
mainLogbook = wxupdatelogbook(mainLogbook,['gate' '_' gateState 's' '_' nameServo],timestampVector);

calllib('phidget21', 'CPhidget_close', handle);
calllib('phidget21', 'CPhidget_delete', handle);