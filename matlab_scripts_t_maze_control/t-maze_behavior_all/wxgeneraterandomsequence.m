%% wxgeneraterandomsequence
 % Generates a random sequence for the learning period
 % JL Alatorre Warren

function sequenceOfSides = wxgeneraterandomsequence(numberOfDecisions)

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
% Column 14: randomized sequence of inhibition runs

% Create cell array sequenceOfSides
sequenceOfSides = cell(numberOfDecisions,14);

% Randomly assign the side (right or left) for each run
consecutiveRight = 0;
consecutiveLeft = 0;
for ii = 1:numberOfDecisions

  % Get a pseudorandom integer: 1 (right gate) or 2 (left gate)
	% Here, rng('shuffle') avoids repeting a result from a previous MATLAB session
  rng('shuffle');
  currentValueSide = randi(2,1);
  
  % Count number 1s or 2s in a row
  if currentValueSide == 1
    consecutiveRight = consecutiveRight+1;
    consecutiveLeft = 0;    
  elseif currentValueSide == 2
    consecutiveLeft = consecutiveLeft+1;
    consecutiveRight = 0;
  end
  
  % If the consecutive counter is greater than 3, force a change
  if consecutiveRight > 3 %3
    currentValueSide = 2;
    consecutiveRight = 0;
    consecutiveLeft = consecutiveLeft + 1;
    sequenceOfSides{ii,2} = 1;
    
  elseif consecutiveLeft > 3 %3
    currentValueSide = 1;
    consecutiveLeft = 0;
    consecutiveRight = consecutiveRight + 1;
    sequenceOfSides{ii,2} = 1;
  else
    sequenceOfSides{ii,2} = 0;
  end
  
  % Fill the vector with 'R's and 'L's
  if currentValueSide == 1
    sequenceOfSides{ii,1} = 'R';
  elseif currentValueSide == 2
    sequenceOfSides{ii,1} = 'L';
  end
  
end

% Randomly assign runs with and without inhibition
consecutiveOn  = 0;
consecutiveOff = 0;
for ii = 1:numberOfDecisions

  % Get a pseudorandom integer: 1 (inhibition on) or 2 (inhibition off)
	% Here, rng('shuffle') avoids repeting a result from a previous MATLAB session
  rng('shuffle');
  currentValueInhibition = randi(2,1);
  
  % Count number 1s or 2s in a row
  if currentValueInhibition == 1
    consecutiveOn = consecutiveOn+1;
    consecutiveOff = 0;    
  elseif currentValueInhibition == 2
    consecutiveOff = consecutiveOff+1;
    consecutiveOn = 0;
  end
  
  % If the consecutive counter is greater than 3, force a change
  if consecutiveOn > 3
    currentValueInhibition = 2;
    consecutiveOn = 0;
    consecutiveOff = consecutiveOff + 1;
    sequenceOfSides{ii,2} = 1;
    
  elseif consecutiveOff > 3
    currentValueInhibition = 1;
    consecutiveOff = 0;
    consecutiveOn = consecutiveOn + 1;
    sequenceOfSides{ii,2} = 1;
  else
    sequenceOfSides{ii,2} = 0;
  end  
  
  % Fill the vector with '1's and '0' for inhibition
  if currentValueInhibition == 1
    sequenceOfSides{ii,14} = 1;
  elseif currentValueInhibition == 2
    sequenceOfSides{ii,14} = 0;
  end
  
end