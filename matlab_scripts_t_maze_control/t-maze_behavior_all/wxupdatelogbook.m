%% wxupdatelogbook
 % Updates logbook with relevant event timestamps
 % JL Alatorre Warren

function mainLogbook = wxupdatelogbook(mainLogbook, eventName, timestampVector)

% Find the next empty element in the cell array mainLogbook
newEventIndex = find(cellfun(@isempty,mainLogbook),1);

% Column 1: Event name
% Column 2: Timestamp as a date vector
% Column 3: Timestamp as a serial date number
mainLogbook{newEventIndex,1} = eventName;
mainLogbook{newEventIndex,2} = timestampVector;
mainLogbook{newEventIndex,3} = datenum(timestampVector);

% Column 4: Elapsed time since the previous event
% Column 5: Elapsed time since the beginning of... time
if newEventIndex > 1
  mainLogbook{newEventIndex,4} = etime(mainLogbook{newEventIndex  ,2}, ...
                                       mainLogbook{newEventIndex-1,2}); 
  mainLogbook{newEventIndex,5} = etime(mainLogbook{newEventIndex  ,2}, ...
                                       mainLogbook{1              ,2});
else
  mainLogbook{newEventIndex,4} = 0;
  mainLogbook{newEventIndex,5} = 0;
end

