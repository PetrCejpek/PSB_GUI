function varargout = InstrumentalProfile(varargin)
% INSTRUMENTALPROFILE MATLAB code for InstrumentalProfile.fig
%      INSTRUMENTALPROFILE, by itself, creates a new INSTRUMENTALPROFILE or raises the existing
%      singleton*.
%
%      H = INSTRUMENTALPROFILE returns the handle to a new INSTRUMENTALPROFILE or the handle to
%      the existing singleton*.
%
%      INSTRUMENTALPROFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSTRUMENTALPROFILE.M with the given input arguments.
%
%      INSTRUMENTALPROFILE('Property','Value',...) creates a new INSTRUMENTALPROFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InstrumentalProfile_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InstrumentalProfile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InstrumentalProfile

% Last Modified by GUIDE v2.5 04-Apr-2024 18:42:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InstrumentalProfile_OpeningFcn, ...
                   'gui_OutputFcn',  @InstrumentalProfile_OutputFcn, ...
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


% --- Executes just before InstrumentalProfile is made visible.
function InstrumentalProfile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InstrumentalProfile (see VARARGIN)

% Choose default command line output for InstrumentalProfile
handles.output = hObject;

% initial values
if isempty(varargin)
    data=[{0} {'+'} {0} {'*'} {'1'};
        {0} {'+'} {0} {'*'} {'1'};
        {0.003} {'+'} {0} {'*'} {'1'};
        {0} {'+'} {0} {'*'} {'1'};
        {0} {'+'} {0} {'*'} {'1'};
        {0} {'+'} {0} {'*'} {'1'};
        {0} {'+'} {0} {'*'} {'1'};];
    Asym2ThetaMax=60;
else
    % Default Output is the Input
    data0=varargin{1};
    data=[{data0.u0} {'+'} {data0.u1} {'*'} {data0.UonChi};
        {data0.v0} {'+'} {data0.v1} {'*'} {data0.VonChi};
        {data0.w0} {'+'} {data0.w1} {'*'} {data0.WonChi};
        {data0.eta00} {'+'} {data0.eta01} {'*'} {data0.Eta0onChi};
        {data0.eta10} {'+'} {data0.eta11} {'*'} {data0.Eta1onChi};
        {data0.asym00} {'+'} {data0.asym01} {'*'} {data0.Asym0onChi};
        {data0.asym10} {'+'} {data0.asym11} {'*'} {data0.Asym1onChi};];
    Asym2ThetaMax=data0.Asym2ThetaMax;
end
handles.ProfileParameters_tab.Data=data;
handles.Asym2ThetaMax_ed.String=num2str(Asym2ThetaMax);

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

handles.OutputData=GetOutput(handles);
guidata(hObject, handles);

% UIWAIT makes InstrumentalProfile wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InstrumentalProfile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.OutputData;

% The GUI can be deleted now
delete(handles.figure1);


% --- Executes on button press in Apply_btn.
function Apply_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.OutputData=GetOutput(handles);
guidata(hObject, handles);
uiresume;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1,'waitstatus'),'waiting')
    Apply_btn_Callback(handles.Apply_btn, eventdata, handles)
%     uiresume;
else
    delete(gcf);
end

function OutputData=GetOutput(handles)
data=handles.ProfileParameters_tab.Data;
OutputData.u0=data{1,1};
OutputData.u1=data{1,3};
OutputData.UonChi=data{1,5};
OutputData.v0=data{2,1};
OutputData.v1=data{2,3};
OutputData.VonChi=data{2,5};
OutputData.w0=data{3,1};
OutputData.w1=data{3,3};
OutputData.WonChi=data{3,5};
OutputData.eta00=data{4,1};
OutputData.eta01=data{4,3};
OutputData.Eta0onChi=data{4,5};
OutputData.eta10=data{5,1};
OutputData.eta11=data{5,3};
OutputData.Eta1onChi=data{5,5};
OutputData.asym00=data{6,1};
OutputData.asym01=data{6,3};
OutputData.Asym0onChi=data{6,5};
OutputData.asym10=data{7,1};
OutputData.asym11=data{7,3};
OutputData.Asym1onChi=data{7,5};
OutputData.Asym2ThetaMax=str2num(handles.Asym2ThetaMax_ed.String);


function Asym2ThetaMax_ed_Callback(hObject, eventdata, handles)
% hObject    handle to Asym2ThetaMax_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Asym2ThetaMax_ed as text
%        str2double(get(hObject,'String')) returns contents of Asym2ThetaMax_ed as a double


% --- Executes during object creation, after setting all properties.
function Asym2ThetaMax_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Asym2ThetaMax_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
