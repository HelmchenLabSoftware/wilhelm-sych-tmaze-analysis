%% wxinhibitionsetup
 % Decides when to turn on and off the inhibition laser based on a decision tree
 % JL Alatorre-Warren
 
function mainLogbook = wxinhibitionsetup(laserSelection, ...
                                         experimentStage, ...
                                         pathFlags, ...
                                         currentInhibitionStatus, ...
                                         mainLogbook)

% Relevant stages:
% 1 InhibitionStageStart 
%   Just before the start box opens during the forced run
% 2 InhibitionStageVertex
%   When the mice is detected at the vertex ROI after licking water
% 3 InhibitionStageDelay
%   Just before the start box opens during the free run
% 4 InhibitionStageEnd
%   When the mice licks water during the free run

% Display current inhibition status
disp(['Inhibition laser status: ' currentInhibitionStatus])

% Select the relevant case
if strcmp(currentInhibitionStatus, 'on') == 1
  switch laserSelection
    case 'Recording'
      switch experimentStage
        case 'InhibitionStageStart'
        case 'InhibitionStageVertex'
        case 'InhibitionStageDelay'
        case 'InhibitionStageEnd'
      end
    case 'Tonic'
      switch experimentStage
        case 'InhibitionStageStart'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');
        case 'InhibitionStageVertex'
        case 'InhibitionStageDelay'
        case 'InhibitionStageEnd'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');                  
      end      
    case 'Encoding'
      switch experimentStage
        case 'InhibitionStageStart'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');
        case 'InhibitionStageVertex'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');        
        case 'InhibitionStageDelay'
        case 'InhibitionStageEnd'
      end
    case 'Delay'
      switch experimentStage
        case 'InhibitionStageStart'
        case 'InhibitionStageVertex'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');        
        case 'InhibitionStageDelay'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');        
        case 'InhibitionStageEnd'
      end
    case 'Retrieval'
      switch experimentStage
        case 'InhibitionStageStart'
        case 'InhibitionStageVertex'
        case 'InhibitionStageDelay'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');        
        case 'InhibitionStageEnd'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');
      end
    case 'E&R'
      switch experimentStage
        case 'InhibitionStageStart'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');
        case 'InhibitionStageVertex'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');        
        case 'InhibitionStageDelay'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'on');        
        case 'InhibitionStageEnd'
          mainLogbook = wxinhibitionlaser(pathFlags,mainLogbook,'off');
      end
  end
end