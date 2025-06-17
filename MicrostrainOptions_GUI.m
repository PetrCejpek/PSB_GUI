function varargout = MicrostrainOptions_GUI(varargin)
% MICROSTRAINOPTIONS_GUI MATLAB code for MicrostrainOptions.fig
%      MICROSTRAINOPTIONS_GUI, by itself, creates a new MICROSTRAINOPTIONS_GUI or raises the existing
%      singleton*.
%
%      H = MICROSTRAINOPTIONS_GUI returns the handle to a new MICROSTRAINOPTIONS_GUI or the handle to
%      the existing singleton*.
%
%      MICROSTRAINOPTIONS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MICROSTRAINOPTIONS_GUI.M with the given input arguments.
%
%      MICROSTRAINOPTIONS_GUI('Property','Value',...) creates a new MICROSTRAINOPTIONS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MicrostrainOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MicrostrainOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MicrostrainOptions_GUI

% Last Modified by GUIDE v2.5 09-May-2025 16:12:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MicrostrainOptions_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MicrostrainOptions_GUI_OutputFcn, ...
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


% --- Executes just before MicrostrainOptions_GUI is made visible.
function MicrostrainOptions_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MicrostrainOptions_GUI (see VARARGIN)

% Choose default command line output for MicrostrainOptions_GUI
handles.output = hObject;

% Import initial values
MicrostrainOptions=varargin{1};
switch MicrostrainOptions.MicrostrainModel
    case 'Isotropic'
        handles.MicrostrainModel_pop.Value=1;
        DATA=[{'e'} {MicrostrainOptions.ParGuess} {MicrostrainOptions.ParFix} {MicrostrainOptions.ParLB} {MicrostrainOptions.ParUB} {'Unitless'}];
    case 'Dislocations (polycrystall)'
        handles.MicrostrainModel_pop.Value=2;
        DATA=[{'e'} {MicrostrainOptions.ParGuess(1)} {MicrostrainOptions.ParFix(1)} {MicrostrainOptions.ParLB(1)} {MicrostrainOptions.ParUB(1)} {'Unitless'};
              {'q'} {MicrostrainOptions.ParGuess(2)} {MicrostrainOptions.ParFix(2)} {MicrostrainOptions.ParLB(2)} {MicrostrainOptions.ParUB(2)} {'Unitless'};];
    case 'Dislocations (specific Burgers vector)'
        handles.MicrostrainModel_pop.Value=3;
        DATA=[{['e']} {MicrostrainOptions.ParGuess(1)} {MicrostrainOptions.ParFix(1)} {MicrostrainOptions.ParLB(1)} {MicrostrainOptions.ParUB(1)} {'Unitless'};
              {['PsiB']} {MicrostrainOptions.ParGuess(2)} {MicrostrainOptions.ParFix(2)} {MicrostrainOptions.ParLB(2)} {MicrostrainOptions.ParUB(2)} {'deg'};
              {['PhiB']} {MicrostrainOptions.ParGuess(3)} {MicrostrainOptions.ParFix(3)} {MicrostrainOptions.ParLB(3)} {MicrostrainOptions.ParUB(3)} {'deg'};];
end

handles.Parameters_tab.Data=DATA; 

% Default output is input
handles.OutputData=varargin{1};

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataManualSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MicrostrainOptions_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutputData;

% The GUI can be deleted now
delete(handles.figure1);


% --- Executes on selection change in MicrostrainModel_pop.
function MicrostrainModel_pop_Callback(hObject, eventdata, handles)
% hObject    handle to MicrostrainModel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MicrostrainModel_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MicrostrainModel_pop
% --- Executes on selection change in MicrostrainModel_pop.
switch hObject.Value
    case 1
        % spherical
        DATA=[{['e']} {NaN} {false} {NaN} {NaN} {'Unitless'}];
    case 2
        DATA=[{'e'} {NaN} {false} {NaN} {NaN} {'Unitless'};
              {'q'} {NaN} {false} {NaN} {NaN} {'Unitless'};];
    case 3
        DATA=[{'e'} {NaN} {false} {NaN} {NaN} {'Unitless'};
              {'PsiB'} {NaN} {false} {NaN} {NaN} {'deg'};
              {'PhiB'} {NaN} {false} {NaN} {NaN} {'deg'};];
end
handles.Parameters_tab.Data=DATA;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MicrostrainModel_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MicrostrainModel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OK_btn.
function OK_btn_Callback(hObject, eventdata, handles)
% hObject    handle to OK_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.OutputData.MicrostrainModel=handles.MicrostrainModel_pop.String{handles.MicrostrainModel_pop.Value};
DATA=handles.Parameters_tab.Data;
handles.OutputData.ParGuess=cell2mat(DATA(:,2))';
handles.OutputData.ParFix=cell2mat(DATA(:,3))';
handles.OutputData.ParLB=cell2mat(DATA(:,4))';
handles.OutputData.ParUB=cell2mat(DATA(:,5))';
guidata(hObject, handles);
uiresume;


% --- Executes on button press in Cancel_btn.
function Cancel_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;


% --- Executes on button press in Help_btn.
function Help_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Help_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg({'You can select parameters for the choosen model:',...
    '',...
    'Guess - guess value for the parameter',...
    'Fix - true/false value, if the parameter is fixed or not',...
    'lb - lower boundary for the parameter',...
    'ub - upper boundary for the parameter',...
    '',...
    'If you use NaN, the value will be estimated directly before the fit. It might be useful for example as the first guess and if the broadening has relatively simple behaviour.'},...
    'Help for model parameters');


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
