%% ccc
 % Close all T-maze's servomotor gates in case of emergency
 % JL Alatorre Warren
 
emergencyLogbook = cell(10,5);
emergencyLogbook = wxservo('start', 'close', emergencyLogbook);
emergencyLogbook = wxservo('right', 'close', emergencyLogbook);
emergencyLogbook = wxservo('left',  'close', emergencyLogbook);