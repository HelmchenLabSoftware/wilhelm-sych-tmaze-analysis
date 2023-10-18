function varargout = wxgui(varargin)
% WXGUI MATLAB code for wxgui.fig
%      WXGUI, by itself, creates a new WXGUI or raises the existing
%      singleton*.
%
%      H = WXGUI returns the handle to a new WXGUI or the handle to
%      the existing singleton*.
%
%      WXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WXGUI.M with the given input arguments.
%
%      WXGUI('Property','Value',...) creates a new WXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wxgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wxgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wxgui

% Last Modified by GUIDE v2.5 24-Feb-2018 23:59:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wxgui_OpeningFcn, ...
                   'gui_OutputFcn',  @wxgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before wxgui is made visible.
function wxgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wxgui (see VARARGIN)

% Choose default command line output for wxgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wxgui wait for user response (see UIRESUME)
% uiwait(handles.MainGui);

% --- Outputs from this function are returned to the command line.
function varargout = wxgui_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in PushButtonLaunch.
function PushButtonLaunch_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonLaunch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This directory will be used to save the "smooth stop" flag
pathSmoothStop = 'I:\2017_m001-m006_backup/smoothstop';

% Delete the smooth stop flag if found
flagSmoothStop = exist([pathSmoothStop '/' 'flagsmoothstop.mat'],'file');
if flagSmoothStop == 2
  delete([pathSmoothStop '/' 'flagsmoothstop.mat']);
end

% Get user input
mouseName                           = get(handles.EditMouseName,'String');
experimentName                      = get(handles.EditExperimentName,'String'); %#ok<NASGU>
maxNumberOfRunsShapingPeriod        = str2double(get(handles.EditNumberOfShapingRuns,'String'));
maxNumberOfRunsLearningPeriod       = str2double(get(handles.EditNumberOfLearningRuns,'String'));
secondsInsideStartBoxAfterForcedRun = str2double(get(handles.EditSecondsInsideStartBoxAfterForcedRun,'String'));
secondsInsideStartBoxAfterFreeRun   = str2double(get(handles.EditSecondsInsideStartBoxAfterFreeRun,'String'));

% Get experiment type (shaping vs learning)
getTypeOfExperiment = get(handles.UiButtonGroupTypeOfExperiment,'SelectedObject');
experimentName = get(getTypeOfExperiment,'String');
if experimentName == 'Learning Sequence (Sample and Choice Runs)'
  experimentName = 'learning';
else
  experimentName = 'shaping';
end

% Get inhibition option (recording, encoding, delay, retrieval and E&R)
getLaserOptions = get(handles.UiButtonGroupLaserOptions,'SelectedObject');
laserSelection = get(getLaserOptions,'String');

% Get image format (MATLAB, JPEG, PNG, TIFF)
getImageFormat = get(handles.UiButtonGroupImageFormat,'SelectedObject');
imageFormat = get(getImageFormat,'String');

% Get additional user options
optionWhiteNoiseInsideStartBoxAfterForcedRun = get(handles.CheckBoxOptionWhiteNoiseInsideStartBoxAfterForcedRun,'Value');
optionWhiteNoiseInsideStartBoxAfterFreeRun   = get(handles.CheckBoxOptionWhiteNoiseInsideStartBoxAfterFreeRun,'Value');
optionUseStartGateAfterFreeingTheMouse       = get(handles.CheckBoxUseStartGate,'Value');
optionLickToOpenValve                        = get(handles.CheckboxLickToOpenValve,'Value');
optionSaveFrames                             = get(handles.CheckboxSaveFrames,'Value');
optionSaveVideo                              = get(handles.CheckboxSaveVideo,'Value');

% Call wxguitmaze, the main program
[~, ~] = wxguitmaze(mouseName, ...
                    experimentName, ...
                    maxNumberOfRunsShapingPeriod, ...
                    maxNumberOfRunsLearningPeriod, ...
                    secondsInsideStartBoxAfterForcedRun, ...
                    secondsInsideStartBoxAfterFreeRun, ... 
                    optionWhiteNoiseInsideStartBoxAfterForcedRun, ...
                    optionWhiteNoiseInsideStartBoxAfterFreeRun, ...  
                    optionUseStartGateAfterFreeingTheMouse, ...
                    optionLickToOpenValve, ...
                    optionSaveFrames, ...
                    optionSaveVideo, ...
                    laserSelection, ...
                    imageFormat);

function EditMouseName_Callback(hObject, eventdata, handles)
% hObject    handle to EditMouseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditMouseName as text
%        str2double(get(hObject,'String')) returns contents of EditMouseName as a double

% --- Executes during object creation, after setting all properties.
function EditMouseName_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to EditMouseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditExperimentName_Callback(hObject, eventdata, handles)
% hObject    handle to EditExperimentName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of EditExperimentName as text
%        str2double(get(hObject,'String')) returns contents of EditExperimentName as a double

% --- Executes during object creation, after setting all properties.
function EditExperimentName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditExperimentName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditNumberOfShapingRuns_Callback(hObject, eventdata, handles)
% hObject    handle to EditNumberOfShapingRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of EditNumberOfShapingRuns as text
%        str2double(get(hObject,'String')) returns contents of EditNumberOfShapingRuns as a double

% --- Executes during object creation, after setting all properties.
function EditNumberOfShapingRuns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditNumberOfShapingRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditNumberOfLearningRuns_Callback(hObject, eventdata, handles)
% hObject    handle to EditNumberOfLearningRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of EditNumberOfLearningRuns as text
%        str2double(get(hObject,'String')) returns contents of EditNumberOfLearningRuns as a double

% --- Executes during object creation, after setting all properties.
function EditNumberOfLearningRuns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditNumberOfLearningRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditSecondsInsideStartBoxAfterForcedRun_Callback(hObject, eventdata, handles)
% hObject    handle to EditSecondsInsideStartBoxAfterForcedRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of EditSecondsInsideStartBoxAfterForcedRun as text
%        str2double(get(hObject,'String')) returns contents of EditSecondsInsideStartBoxAfterForcedRun as a double

% --- Executes during object creation, after setting all properties.
function EditSecondsInsideStartBoxAfterForcedRun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSecondsInsideStartBoxAfterForcedRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EditSecondsInsideStartBoxAfterFreeRun_Callback(hObject, eventdata, handles)
% hObject    handle to EditSecondsInsideStartBoxAfterFreeRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of EditSecondsInsideStartBoxAfterFreeRun as text
%        str2double(get(hObject,'String')) returns contents of EditSecondsInsideStartBoxAfterFreeRun as a double

% --- Executes during object creation, after setting all properties.
function EditSecondsInsideStartBoxAfterFreeRun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSecondsInsideStartBoxAfterFreeRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CheckBoxOptionWhiteNoiseInsideStartBoxAfterForcedRun.
function CheckBoxOptionWhiteNoiseInsideStartBoxAfterForcedRun_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxOptionWhiteNoiseInsideStartBoxAfterForcedRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxOptionWhiteNoiseInsideStartBoxAfterForcedRun

% --- Executes on button press in CheckBoxOptionWhiteNoiseInsideStartBoxAfterFreeRun.
function CheckBoxOptionWhiteNoiseInsideStartBoxAfterFreeRun_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxOptionWhiteNoiseInsideStartBoxAfterFreeRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxOptionWhiteNoiseInsideStartBoxAfterFreeRun

% --- Executes on button press in CheckBoxUseStartGate.
function CheckBoxUseStartGate_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxUseStartGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxUseStartGate

% --- Executes on button press in CheckboxLickToOpenValve.
function CheckboxLickToOpenValve_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxLickToOpenValve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckboxLickToOpenValve

% --- Executes on button press in RadioButtonShapingSequence.
function RadioButtonShapingSequence_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonShapingSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonShapingSequence

% --- Executes on button press in RadioButtonLearningSequence.
function RadioButtonLearningSequence_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonLearningSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonLearningSequence

% --- Executes during object creation, after setting all properties.
function UiButtonGroupTypeOfExperiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UiButtonGroupTypeOfExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in PushButtonSmoothStop.
function PushButtonSmoothStop_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonSmoothStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Send flag to wxtmaze
pathSmoothStop = 'I:\2017_m001-m006_backup/smoothstop';
flagSmoothStop = true; %#ok<NASGU>
save([pathSmoothStop '/' 'flagsmoothstop.mat'], 'flagSmoothStop');

% --- Executes on button press in PushButtonCloseAllGates.
function PushButtonCloseAllGates_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonCloseAllGates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ccc

% --- Executes on button press in PushButtonOpenAllGates.
function PushButtonOpenAllGates_Callback(hObject, eventdata, handles)
% hObject    handle to PushButtonOpenAllGates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ooo

% --- Executes on button press in RadioButtonRecording.
function RadioButtonRecording_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonRecording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonRecording

% --- Executes on button press in RadioButtonEncoding.
function RadioButtonEncoding_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonEncoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonEncoding

% --- Executes on button press in RadioButtonDelay.
function RadioButtonDelay_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonDelay

% --- Executes on button press in RadioButtonRetrieval.
function RadioButtonRetrieval_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonRetrieval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonRetrieval

% --- Executes on button press in RadioButtonEncodingAndRetrieval.
function RadioButtonEncodingAndRetrieval_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonEncodingAndRetrieval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonEncodingAndRetrieval

% --- Executes on button press in RadioButtonTonic.
function RadioButtonTonic_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonTonic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonTonic

% --- Executes during object creation, after setting all properties.
function UiButtonGroupLaserOptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UiButtonGroupLaserOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function UiButtonGroupImageFormat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UiButtonGroupImageFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in RadioButtonMatlab.
function RadioButtonMatlab_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonMatlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonMatlab

% --- Executes on button press in RadioButtonJpeg.
function RadioButtonJpeg_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonJpeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonJpeg

% --- Executes on button press in RadioButtonPng.
function RadioButtonPng_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonPng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonPng

% --- Executes on button press in RadioButtonTiff.
function RadioButtonTiff_Callback(hObject, eventdata, handles)
% hObject    handle to RadioButtonTiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of RadioButtonTiff

% --- Executes on button press in CheckboxSaveFrames.
function CheckboxSaveFrames_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxSaveFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckboxSaveFrames

% --- Executes on button press in CheckboxSaveVideo.
function CheckboxSaveVideo_Callback(hObject, eventdata, handles)
% hObject    handle to CheckboxSaveVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckboxSaveVideo
