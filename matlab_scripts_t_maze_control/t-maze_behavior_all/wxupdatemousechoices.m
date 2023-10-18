%% wxupdatemousechoices
 % Updates the cell array with the sides chosen by the mouse
 % JL Alatorre Warren

function sequenceOfSides = wxupdatemousechoices(sequenceOfSides, ...
                                                counterOfRuns, ...
                                                detectedSide)

% Notes about sequenceOfSides
% Column 01: randomized sequence of sides
% Column 02: flags with forced changes to avoid 4 or more consecutives
% Column 03: sides chosen by the mouse
% Column 04: correct (1) or wrong (0) choice
% Column 05: cummulative correct choices
% Column 06: current success rate
% Column 07: timestamps: start forced run (date and time)
% Column 08: timestamps: stop forced run (date and time)
% Column 09: timestamps: start free (open) run (date and time)
% Column 10: timestamps: stop free (open) run (date and time)
% Column 11: elapsed time (seconds): forced run
% Column 12: elapsed time (seconds): free run
% Column 13: elapsed time (seconds): full trial (forced run + free run)
                                              
% Write in sequenceOfSides the side chosen by the mouse
if strcmp(detectedSide, 'right') == 1
  sequenceOfSides{counterOfRuns,3} = 'R';
elseif strcmp(detectedSide, 'left') == 1
  sequenceOfSides{counterOfRuns,3} = 'L';
end

% Write correct (1) or wrong (0) choice
% Wrong choice (0): the mouse chose the same gate
% Good choice (1): the mouse chose the opposite gate
if sequenceOfSides{counterOfRuns,1} == sequenceOfSides{counterOfRuns,3}
  sequenceOfSides{counterOfRuns,4} = 0;
else
  sequenceOfSides{counterOfRuns,4} = 1;
end

% Update cummulative correct choices
if counterOfRuns == 1
  sequenceOfSides{counterOfRuns,5} = sequenceOfSides{counterOfRuns,4};
else
  sequenceOfSides{counterOfRuns,5} = sequenceOfSides{counterOfRuns,4} + sequenceOfSides{counterOfRuns-1,5};
end

% Update current success rate
sequenceOfSides{counterOfRuns,6} = sequenceOfSides{counterOfRuns,5}/counterOfRuns;