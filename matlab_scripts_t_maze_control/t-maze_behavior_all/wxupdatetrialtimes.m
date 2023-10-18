%% wxupdatetrialtimes
 % Updates the trial times in the cell array with the sides chosen by the mouse
 % JL Alatorre Warren

function sequenceOfSides = wxupdatetrialtimes(sequenceOfSides, ...
                                              counterOfRuns, ...
                                              startTimeForcedRun, ...
                                              endTimeForcedRun, ...
                                              startTimeFreeRun, ...
                                              endTimeFreeRun)                                            

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

% Write timestamps in sequenceOfSides
sequenceOfSides{counterOfRuns, 7} = startTimeForcedRun;
sequenceOfSides{counterOfRuns, 8} = endTimeForcedRun;
sequenceOfSides{counterOfRuns, 9} = startTimeFreeRun;
sequenceOfSides{counterOfRuns,10} = endTimeFreeRun;
sequenceOfSides{counterOfRuns,11} = etime(endTimeForcedRun,startTimeForcedRun);
sequenceOfSides{counterOfRuns,12} = etime(endTimeFreeRun  ,startTimeFreeRun);
sequenceOfSides{counterOfRuns,13} = sequenceOfSides{counterOfRuns,11} + sequenceOfSides{counterOfRuns,12};