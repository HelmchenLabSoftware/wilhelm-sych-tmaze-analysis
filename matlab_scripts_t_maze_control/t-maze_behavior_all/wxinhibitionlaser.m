%% wxinhibitionlaser
 % Turns on and off the laser for inhibition experiments when needed
 % JL Alatorre-Warren
 
function mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,laserStatus)

switch laserStatus
  case 'on'
    flagInhibitionLaserOn = true; %#ok<NASGU>
    save([pathFlags '/' 'flaginhibitionlaseron.mat' ],'flagInhibitionLaserOn');
    mainLogbook = wxupdatelogbook(mainLogbook, 'inhibitionlaser_on', clock);
  case 'off'
    flagInhibitionLaserOff = false; %#ok<NASGU>
    save([pathFlags '/' 'flaginhibitionlaseroff.mat'],'flagInhibitionLaserOff');
    mainLogbook = wxupdatelogbook(mainLogbook, 'inhibitionlaser_off', clock);
end