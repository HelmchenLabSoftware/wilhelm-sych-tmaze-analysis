%% wxvisualizechoicesandtrialtimes
 % Updates the plot of the mouse choices
 % JL Alatorre-Warren

function wxvisualizechoicesandtrialtimes(sequenceOfSides, ...
                                         counterOfRuns, ...
                                         maxNumberOfRunsLearningPeriod)
                          
% Preallocate space in memory for listOfMouseChoices and listOfTrialTimes
listOfMouseChoices = zeros(maxNumberOfRunsLearningPeriod,1);
listOfTrialTimes   = zeros(maxNumberOfRunsLearningPeriod,2);
                          
% Fill listOfMouseChoices with the choices stored in the cell array sequenceOfSides
for ii = 1:counterOfRuns
  listOfMouseChoices(ii) = sequenceOfSides{ii,4};
  if listOfMouseChoices(ii) == 0
    listOfMouseChoices(ii) = -1;
  end
end

% Fill listOfTrialTimes with the times stored in the cell array sequenceOfSides
for ii = 1:counterOfRuns
  
  % Forced run time
  listOfTrialTimes(ii,1) = sequenceOfSides{ii,11};
  
  % Free run time  
  listOfTrialTimes(ii,2) = sequenceOfSides{ii,12};
  
end
  
% Create new figure
% Place it in the desired position
figureChoices = figure(2);
figureChoices.OuterPosition = [1 741 1280 700];

% Create a bar subplot: mouse choices
subplot(2,1,1)

for ii=1:counterOfRuns
if sequenceOfSides{ii,1}=='L'&&sequenceOfSides{ii,14}==0
bar(ii,listOfMouseChoices(ii),'b')
hold on
elseif sequenceOfSides{ii,1}=='L'&&sequenceOfSides{ii,14}==1
bar(ii,listOfMouseChoices(ii),'c')
hold on
elseif sequenceOfSides{ii,1}=='R'&&sequenceOfSides{ii,14}==0
bar(ii,listOfMouseChoices(ii),'r')
hold on
elseif sequenceOfSides{ii,1}=='R'&&sequenceOfSides{ii,14}==1
bar(ii,listOfMouseChoices(ii),'m')
hold on
end
end
% colormap(cool)
title({'Correct Versus Incorrect Mouse Choices';'Blue-Left Sample, Red-Right Sample'})
xlim([0 (maxNumberOfRunsLearningPeriod+1)])
ylim([-1.25 1.25])
set(gca,'YTickLabel',[]);
pause(0.200)

% Create a bar subplot: trial times
subplot(2,1,2)
if counterOfRuns == 1
  bar(listOfTrialTimes(1:counterOfRuns,:))
elseif counterOfRuns > 1
  bar(listOfTrialTimes(1:counterOfRuns,:),'stacked')
end
colormap(cool)
title('Trial Times')
xlim([0 (maxNumberOfRunsLearningPeriod+1)])
ylim('auto')
pause(0.200)

% Notify the current performance of the mouse
disp(['Current performance: ' num2str(sequenceOfSides{counterOfRuns,5}) ' '...
      'out of ' num2str(counterOfRuns) ' trials '...
      '(' num2str(100*sequenceOfSides{counterOfRuns,6}) '%)'])