%% wxvisualizealternationandtrialtimes
 % Updates the plot of the mouse alternation task
 % JL Alatorre-Warren

function wxvisualizealternationandtrialtimes(sequenceOfAlternation, ...
                                             counterOfRuns, ...
                                             maxNumberOfRunsShapingPeriod)
                          
% Preallocate space in memory for listOfAlternation and listOfTrialTimes
listOfAlternation = zeros(maxNumberOfRunsShapingPeriod,1);
listOfTrialTimes  = zeros(maxNumberOfRunsShapingPeriod,1);
                          
% Fill listOfAlternation with the choices stored in the cell array sequenceOfAlternation
for ii = 1:counterOfRuns
  listOfAlternation(ii) = sequenceOfAlternation{ii,1};
  if listOfAlternation(ii) == 'R'
    listOfAlternation(ii) = 1;    
  elseif listOfAlternation(ii) == 'L'
    listOfAlternation(ii) = -1;
  end
end

% Fill listOfTrialTimes with the times stored in the cell array sequenceOfAlternation
for ii = 1:counterOfRuns
  
  % Forced run time
  listOfTrialTimes(ii,1) = sequenceOfAlternation{ii,11};
  
end
  
% Create new figure
% Place it in the desired position
figureChoices = figure(2);
figureChoices.OuterPosition = [1 741 1280 700];

% Create a bar subplot: alternation sides
subplot(2,1,1)
bar(listOfAlternation(1:counterOfRuns))
colormap(cool)
title('Alternation task')
xlim([0 (maxNumberOfRunsShapingPeriod+1)])
ylim([-1.25 1.25])
set(gca,'YTickLabel',[]);
pause(0.200)

% Create a bar subplot: trial times
subplot(2,1,2)
bar(listOfTrialTimes(1:counterOfRuns,1))
colormap(cool)
title('Trial Times')
xlim([0 (maxNumberOfRunsShapingPeriod+1)])
ylim('auto')
pause(0.200)