function wxwritesignals(~,event,pathSignals)

  % Get the data that is currently available in the DAQ device
  daqSignals = event.Data; %#ok<NASGU>

  % Get timestamp
  timestampVector = clock;
  
  % Save the data that is currently available
  save([pathSignals '/' 'daqSignals' ...
        datestr(timestampVector,'yyyymmddHHMMSSFFF') '.mat'], ...
        'daqSignals');

  % Save the data that is currently available
  save([pathSignals '/' 'daqSignalsMostRecent.mat'], 'daqSignals');
  
  % Display
  disp(datestr(timestampVector,'yyyymmddHHMMSSFFF'))
        
end 