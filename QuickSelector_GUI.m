function varargout = QuickSelector_GUI(varargin)
% QUICKSELECTOR_GUI MATLAB code for QuickSelector_GUI.fig
%      QUICKSELECTOR_GUI, by itself, creates a new QUICKSELECTOR_GUI or raises the existing
%      singleton*.
%
%      H = QUICKSELECTOR_GUI returns the handle to a new QUICKSELECTOR_GUI or the handle to
%      the existing singleton*.
%
%      QUICKSELECTOR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICKSELECTOR_GUI.M with the given input arguments.
%
%      QUICKSELECTOR_GUI('Property','Value',...) creates a new QUICKSELECTOR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QuickSelector_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QuickSelector_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QuickSelector_GUI

% Last Modified by GUIDE v2.5 19-Aug-2025 12:27:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QuickSelector_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @QuickSelector_GUI_OutputFcn, ...
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


% --- Executes just before QuickSelector_GUI is made visible.
function QuickSelector_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QuickSelector_GUI (see VARARGIN)

% Choose default command line output for QuickSelector_GUI
% Choose default command line output for MicrostrainOptions_GUI
handles.output = hObject;

% Import initial values
handles.Labels=varargin{1}.Labels;
% find unique labels
UniqueLabels=unique(handles.Labels);
handles.LabelString_pop.String=[UniqueLabels; 'my custom string'];
handles.LabelString_pop.Value=1;

% Default output is input
handles.OutputData.NewSelection=varargin{1}.CurrentSelection;

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataManualSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QuickSelector_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutputData;

% The GUI can be deleted now
delete(handles.figure1);


% --- Executes on selection change in LabelOption_pop.
function LabelOption_pop_Callback(hObject, eventdata, handles)
% hObject    handle to LabelOption_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LabelOption_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LabelOption_pop


% --- Executes during object creation, after setting all properties.
function LabelOption_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelOption_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MyString_ed_Callback(hObject, eventdata, handles)
% hObject    handle to MyString_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MyString_ed as text
%        str2double(get(hObject,'String')) returns contents of MyString_ed as a double


% --- Executes during object creation, after setting all properties.
function MyString_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MyString_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Process_btn.
function Process_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Process_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.LabelString_pop.String{handles.LabelString_pop.Value},'my custom string')
    MyString=handles.MyString_ed.String;
else
    MyString=handles.LabelString_pop.String{handles.LabelString_pop.Value};
end

switch handles.LabelOption_pop.String{handles.LabelOption_pop.Value}
    case 'containing'
        % contains(handles.Labels,MyString)
        handles.OutputData.NewSelection=contains(handles.Labels,MyString);
    case 'being exactly'
        % strcmp(handles.Labels,MyString)
        handles.OutputData.NewSelection=strcmp(handles.Labels,MyString);
end

guidata(hObject, handles);
uiresume;


% --- Executes on button press in Cancel_btn.
function Cancel_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;


% --- Executes on selection change in LabelString_pop.
function LabelString_pop_Callback(hObject, eventdata, handles)
% hObject    handle to LabelString_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LabelString_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LabelString_pop
if strcmp(hObject.String{hObject.Value},'my custom string')
    handles.MyString_ed.Enable='on';
else
    handles.MyString_ed.Enable='off';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function LabelString_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelString_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1,'waitstatus'),'waiting')
    uiresume;
else
    delete(gcf);
end
