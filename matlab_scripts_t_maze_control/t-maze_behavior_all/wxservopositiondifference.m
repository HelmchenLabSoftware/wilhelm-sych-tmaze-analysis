%% wxservopositiondifference
 % Workhorse function used to move the T-maze's servomotor gates
 % JL Alatorre Warren

function expectationVsRealityDifference = wxservopositiondifference(nameServo, gateState)

% Time window in seconds used to open the servo motor
timeWindow = 0.500;

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
      expectedServoPosition = 105;
    elseif indexServo == 4
      expectedServoPosition = 135;
    elseif indexServo == 2
      expectedServoPosition = 87;
    end
  case 'open'
    if indexServo == 0
      expectedServoPosition = 42;
    elseif indexServo == 4
      expectedServoPosition = 195;
    elseif indexServo == 2
      expectedServoPosition = 35;
    end
end

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

% CPhidget_waitForAttachment waits for attachment to happen.
% This can be called wirght after calling CPhidget_open,
% as an alternative to using the attach handler.
if calllib('phidget21', 'CPhidget_waitForAttachment', handle, 2500) == 0

  % The variable t is used to control the timing of the opening and 
  % closing of the gates
  t = timer('TimerFcn','disp(''Moving servo'')', 'StartDelay', timeWindow);
  
  % Power the servomotor on
  calllib('phidget21', 'CPhidgetAdvancedServo_setEngaged', handle, indexServo, 1);
  
  % Get the current servomotor position and compare this to its expected value
  calllib('phidget21', 'CPhidgetAdvancedServo_getPosition', handle, indexServo, valPtr);
  currentServoPosition = valPtr.Value;
  expectationVsRealityDifference = expectedServoPosition - currentServoPosition;
  disp(['Difference: ' num2str(expectationVsRealityDifference)])
  
  % Turn the servomotor off
  calllib('phidget21', 'CPhidgetAdvancedServo_setEngaged', handle, indexServo, 0);

else
  disp('Nel')
end

% Close and delete the phidget21 session
calllib('phidget21', 'CPhidget_close', handle);
calllib('phidget21', 'CPhidget_delete', handle);