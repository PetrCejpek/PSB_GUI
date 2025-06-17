function varargout = SizeOptions_GUI(varargin)
% SIZEOPTIONS_GUI MATLAB code for SizeOptions_GUI.fig
%      SIZEOPTIONS_GUI, by itself, creates a new SIZEOPTIONS_GUI or raises the existing
%      singleton*.
%
%      H = SIZEOPTIONS_GUI returns the handle to a new SIZEOPTIONS_GUI or the handle to
%      the existing singleton*.
%
%      SIZEOPTIONS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIZEOPTIONS_GUI.M with the given input arguments.
%
%      SIZEOPTIONS_GUI('Property','Value',...) creates a new SIZEOPTIONS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SizeOptions_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SizeOptions_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SizeOptions_GUI

% Last Modified by GUIDE v2.5 09-May-2025 13:14:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SizeOptions_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SizeOptions_GUI_OutputFcn, ...
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


% --- Executes just before SizeOptions_GUI is made visible.
function SizeOptions_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SizeOptions_GUI (see VARARGIN)

% Choose default command line output for SizeOptions_GUI
handles.output = hObject;

% Import initial values
SizeOptions=varargin{1};
switch SizeOptions.GrainsShape
    case 'Sphere'
        handles.GrainsShape_pop.Value=1;
        DATA=[{'D'} {SizeOptions.ParGuess} {SizeOptions.ParFix} {SizeOptions.ParLB} {SizeOptions.ParUB} {'Angstroem'}];
    case 'Cylinder'
        handles.GrainsShape_pop.Value=2;
        DATA=[{['D' char(8869)]} {SizeOptions.ParGuess(1)} {SizeOptions.ParFix(1)} {SizeOptions.ParLB(1)} {SizeOptions.ParUB(1)} {'Angstroem'};
              {['D' char(8741)]} {SizeOptions.ParGuess(2)} {SizeOptions.ParFix(2)} {SizeOptions.ParLB(2)} {SizeOptions.ParUB(2)} {'Angstroem'};];
    case 'Cylinder with elliptical base'
        handles.GrainsShape_pop.Value=3;
        DATA=[{['D1']} {SizeOptions.ParGuess(1)} {SizeOptions.ParFix(1)} {SizeOptions.ParLB(1)} {SizeOptions.ParUB(1)} {'Angstroem'};
              {['D2']} {SizeOptions.ParGuess(2)} {SizeOptions.ParFix(2)} {SizeOptions.ParLB(2)} {SizeOptions.ParUB(2)} {'Angstroem'};
              {['D3']} {SizeOptions.ParGuess(3)} {SizeOptions.ParFix(3)} {SizeOptions.ParLB(3)} {SizeOptions.ParUB(3)} {'Angstroem'};];
    case 'Rotational ellipsoid'
        handles.GrainsShape_pop.Value=4;
        DATA=[{['D' char(8869)]} {SizeOptions.ParGuess(1)} {SizeOptions.ParFix(1)} {SizeOptions.ParLB(1)} {SizeOptions.ParUB(1)} {'Angstroem'};
              {['D' char(8741)]} {SizeOptions.ParGuess(2)} {SizeOptions.ParFix(2)} {SizeOptions.ParLB(2)} {SizeOptions.ParUB(2)} {'Angstroem'};];
    case 'General ellipsoid'
        handles.GrainsShape_pop.Value=5;
        DATA=[{['D1']} {SizeOptions.ParGuess(1)} {SizeOptions.ParFix(1)} {SizeOptions.ParLB(1)} {SizeOptions.ParUB(1)} {'Angstroem'};
              {['D2']} {SizeOptions.ParGuess(2)} {SizeOptions.ParFix(2)} {SizeOptions.ParLB(2)} {SizeOptions.ParUB(2)} {'Angstroem'};
              {['D3']} {SizeOptions.ParGuess(3)} {SizeOptions.ParFix(3)} {SizeOptions.ParLB(3)} {SizeOptions.ParUB(3)} {'Angstroem'};];
end
switch SizeOptions.GrownDirectionOpt
    case 'Sample normal'
        handles.Direction_pop.Value=1;
        handles.PerpendicularDirection_txt2.String='Direction in the sample surface specified with azimuth';
    case 'Specific crystallographic direction'
        handles.Direction_pop.Value=2;
        handles.PerpendicularDirection_txt2.String='Specific crystallographic direction';
end
handles.Parameters_tab.Data=DATA;
handles.DirectionH_ed.String=num2str(SizeOptions.GrownDirectionHKL(1));
handles.DirectionK_ed.String=num2str(SizeOptions.GrownDirectionHKL(2));
handles.DirectionL_ed.String=num2str(SizeOptions.GrownDirectionHKL(3));
handles.DirectionH2_ed.String=num2str(SizeOptions.GrownDirectionHKL2(1));
handles.DirectionK2_ed.String=num2str(SizeOptions.GrownDirectionHKL2(2));
handles.DirectionL2_ed.String=num2str(SizeOptions.GrownDirectionHKL2(3));
handles.Azimuth_ed.String=num2str(SizeOptions.GrownDirectionPerpAzimuth);
guidata(hObject, handles);
%%%% input parameters for the specific directions
switch SizeOptions.GrainsShape
    case 'Sphere'
        handles.Direction_pop.Enable='off';
        handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
        handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
        handles.PerpendicularDirection_txt2.Enable='off';
        handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
        handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        handles.Azimuth_ed.Enable='off';
        handles.Azimuth_txt.Enable='off';
    case {'Cylinder'; 'Rotational ellipsoid'}
        handles.Direction_pop.Enable='on';
        if handles.Direction_pop.Value==1
            handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
            handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
        else
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
        end
        handles.PerpendicularDirection_txt2.Enable='off';
        handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
        handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        handles.Azimuth_ed.Enable='off';
        handles.Azimuth_txt.Enable='off';
    case {'Cylinder with elliptical base'; 'General ellipsoid'}
        handles.Direction_pop.Enable='on';
        if handles.Direction_pop.Value==1
            handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
            handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
            handles.PerpendicularDirection_txt2.Enable='on';
            handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
            handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
            handles.Azimuth_ed.Enable='on';
            handles.Azimuth_txt.Enable='on';        
        else
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
            handles.PerpendicularDirection_txt2.Enable='on';
            handles.DirectionH2_txt.Enable='on'; handles.DirectionK2_txt.Enable='on'; handles.DirectionL2_txt.Enable='on';
            handles.DirectionH2_ed.Enable='on'; handles.DirectionK2_ed.Enable='on'; handles.DirectionL2_ed.Enable='on';
            handles.Azimuth_ed.Enable='off';
            handles.Azimuth_txt.Enable='off'; 
        end
end

% Default output is input
handles.OutputData=varargin{1};

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataManualSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SizeOptions_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutputData;

% The GUI can be deleted now
delete(handles.figure1);


% --- Executes on selection change in GrainsShape_pop.
function GrainsShape_pop_Callback(hObject, eventdata, handles)
% hObject    handle to GrainsShape_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GrainsShape_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GrainsShape_pop
switch hObject.Value
    case 1
        % spherical
        DATA=[{['D']} {NaN} {false} {NaN} {NaN} {'Angstroem'}];
        handles.Direction_pop.Enable='off';
        handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
        handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
        handles.PerpendicularDirection_txt2.Enable='off';
        handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
        handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        handles.Azimuth_ed.Enable='off';
        handles.Azimuth_txt.Enable='off';
    case {2; 4}
        DATA=[{['D' char(8869)]} {NaN} {false} {NaN} {NaN} {'Angstroem'};
              {['D' char(8741)]} {NaN} {false} {NaN} {NaN} {'Angstroem'};];
        handles.Direction_pop.Enable='on';
        if handles.Direction_pop.Value==1
            handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
            handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
        else
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
        end
        handles.PerpendicularDirection_txt2.Enable='off';
        handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
        handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        handles.Azimuth_ed.Enable='off';
        handles.Azimuth_txt.Enable='off';
    case {3; 5}
        DATA=[{['D1']} {NaN} {false} {NaN} {NaN} {'Angstroem'};
              {['D2']} {NaN} {false} {NaN} {NaN} {'Angstroem'};
              {['D3']} {NaN} {false} {NaN} {NaN} {'Angstroem'}];
        handles.Direction_pop.Enable='on';
        if handles.Direction_pop.Value==1
            handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
            handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
            handles.PerpendicularDirection_txt2.Enable='on';
            handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
            handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
            handles.Azimuth_ed.Enable='on';
            handles.Azimuth_txt.Enable='on';        
        else
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
            handles.PerpendicularDirection_txt2.Enable='on';
            handles.DirectionH2_txt.Enable='on'; handles.DirectionK2_txt.Enable='on'; handles.DirectionL2_txt.Enable='on';
            handles.DirectionH2_ed.Enable='on'; handles.DirectionK2_ed.Enable='on'; handles.DirectionL2_ed.Enable='on';
            handles.Azimuth_ed.Enable='off';
            handles.Azimuth_txt.Enable='off'; 
        end
end
handles.Parameters_tab.Data=DATA;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function GrainsShape_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GrainsShape_pop (see GCBO)
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
handles.OutputData.GrainsShape=handles.GrainsShape_pop.String{handles.GrainsShape_pop.Value};
DATA=handles.Parameters_tab.Data;
handles.OutputData.ParGuess=cell2mat(DATA(:,2))';
handles.OutputData.ParFix=cell2mat(DATA(:,3))';
handles.OutputData.ParLB=cell2mat(DATA(:,4))';
handles.OutputData.ParUB=cell2mat(DATA(:,5))';
handles.OutputData.GrownDirectionOpt=handles.Direction_pop.String{handles.Direction_pop.Value};
handles.OutputData.GrownDirectionHKL=[str2num(handles.DirectionH_ed.String) str2num(handles.DirectionK_ed.String) str2num(handles.DirectionL_ed.String)];
% handles.OutputData.GrownDirectionOpt2=handles.PerpendicularDirection_txt2.String{handles.PerpendicularDirection_txt2.Value};
handles.OutputData.GrownDirectionHKL2=[str2num(handles.DirectionH2_ed.String) str2num(handles.DirectionK2_ed.String) str2num(handles.DirectionL2_ed.String)];
handles.OutputData.GrownDirectionPerpAzimuth=str2num(handles.Azimuth_ed.String);
guidata(hObject, handles);
uiresume;


% --- Executes on button press in Cancel_btn.
function Cancel_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Here, 'uiresume' should be sufficient, because Input data are in Output
% as default
uiresume;


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



function DirectionH_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionH_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionH_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionH_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionH_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionH_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DirectionK_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionK_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionK_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionK_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionK_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionK_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DirectionL_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionL_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionL_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionL_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionL_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionL_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Direction_pop.
function Direction_pop_Callback(hObject, eventdata, handles)
% hObject    handle to Direction_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Direction_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Direction_pop
if handles.Direction_pop.Value==1
    handles.PerpendicularDirection_txt2.String='Direction in sample surface specified by azimuth';
    handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
    handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
    handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
    handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
    switch handles.GrainsShape_pop.Value
        case {1; 2; 4}
            handles.Azimuth_ed.Enable='off';
            handles.Azimuth_txt.Enable='off';
            handles.PerpendicularDirection_txt2.Enable='off';
        case {3; 5}
            handles.Azimuth_ed.Enable='on';
            handles.Azimuth_txt.Enable='on';
            handles.PerpendicularDirection_txt2.Enable='on';
    end
else
    handles.PerpendicularDirection_txt2.String='Specific crystallographic direction';
    handles.Azimuth_ed.Enable='off';
    handles.Azimuth_txt.Enable='off';
    handles.PerpendicularDirection_txt2.Enable='off';
    switch handles.GrainsShape_pop.Value
        case 1 % sphere (this in fact cannot happen
            handles.DirectionH_txt.Enable='off'; handles.DirectionK_txt.Enable='off'; handles.DirectionL_txt.Enable='off';
            handles.DirectionH_ed.Enable='off'; handles.DirectionK_ed.Enable='off'; handles.DirectionL_ed.Enable='off';
            handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
            handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        case {2; 4} % cylinder, rotational ellipsoid
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
            handles.DirectionH2_txt.Enable='off'; handles.DirectionK2_txt.Enable='off'; handles.DirectionL2_txt.Enable='off';
            handles.DirectionH2_ed.Enable='off'; handles.DirectionK2_ed.Enable='off'; handles.DirectionL2_ed.Enable='off';
        case {3; 5} % cylinder with elliptical base, general ellipsoid
            handles.DirectionH_txt.Enable='on'; handles.DirectionK_txt.Enable='on'; handles.DirectionL_txt.Enable='on';
            handles.DirectionH_ed.Enable='on'; handles.DirectionK_ed.Enable='on'; handles.DirectionL_ed.Enable='on';
            handles.DirectionH2_txt.Enable='on'; handles.DirectionK2_txt.Enable='on'; handles.DirectionL2_txt.Enable='on';
            handles.DirectionH2_ed.Enable='on'; handles.DirectionK2_ed.Enable='on'; handles.DirectionL2_ed.Enable='on';
    end
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Direction_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Direction_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DirectionH2_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionH2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionH2_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionH2_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionH2_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionH2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DirectionK2_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionK2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionK2_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionK2_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionK2_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionK2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DirectionL2_ed_Callback(hObject, eventdata, handles)
% hObject    handle to DirectionL2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectionL2_ed as text
%        str2double(get(hObject,'String')) returns contents of DirectionL2_ed as a double


% --- Executes during object creation, after setting all properties.
function DirectionL2_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectionL2_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PerpendicularDirection_txt2.
function PerpendicularDirection_txt2_Callback(hObject, eventdata, handles)
% hObject    handle to PerpendicularDirection_txt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PerpendicularDirection_txt2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PerpendicularDirection_txt2


% --- Executes during object creation, after setting all properties.
function PerpendicularDirection_txt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PerpendicularDirection_txt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Azimuth_ed_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuth_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Azimuth_ed as text
%        str2double(get(hObject,'String')) returns contents of Azimuth_ed as a double


% --- Executes during object creation, after setting all properties.
function Azimuth_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimuth_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
