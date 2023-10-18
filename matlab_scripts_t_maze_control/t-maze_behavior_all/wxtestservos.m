%% wxtestservos
 % This function tests the Phidget Advanced Servomotor, which controls the 
 % opening and closing of the gates
 % JL Alatorre Warren

%{
  Servomotor details: 
  Phidgets's Advanced Servo Motor
  1061_1 - PhidgetAdvancedServo 8-Motor
%}

function mainLogbook = wxtestservos(numberOfCycles,mainLogbook)

if isscalar(numberOfCycles) == 0
  numberOfCycles = 1;
end

n=1;
while n <= numberOfCycles 
  mainLogbook = wxservo('start','open' ,mainLogbook);
  mainLogbook = wxservo('start','close',mainLogbook);
  mainLogbook = wxservo('right','open' ,mainLogbook);
  mainLogbook = wxservo('right','close',mainLogbook);
  mainLogbook = wxservo('left' ,'open' ,mainLogbook);
  mainLogbook = wxservo('left' ,'close',mainLogbook);
  n=n+1;
end