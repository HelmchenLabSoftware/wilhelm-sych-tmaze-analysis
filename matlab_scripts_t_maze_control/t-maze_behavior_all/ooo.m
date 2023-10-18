%% ooo
 % Open all T-maze's servomotor gates in case of emergency
 % JL Alatorre Warren
 
emergencyLogbook = cell(10,5);
emergencyLogbook = wxservo('start', 'open', emergencyLogbook);
emergencyLogbook = wxservo('right', 'open', emergencyLogbook);
emergencyLogbook = wxservo('left',  'open', emergencyLogbook);