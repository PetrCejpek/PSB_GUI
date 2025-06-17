function varargout = DataManualSelection(varargin)
% DATAMANUALSELECTION MATLAB code for DataManualSelection.fig
%      DATAMANUALSELECTION, by itself, creates a new DATAMANUALSELECTION or raises the existing
%      singleton*.
%
%      H = DATAMANUALSELECTION returns the handle to a new DATAMANUALSELECTION or the handle to
%      the existing singleton*.
%
%      DATAMANUALSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATAMANUALSELECTION.M with the given input arguments.
%
%      DATAMANUALSELECTION('Property','Value',...) creates a new DATAMANUALSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataManualSelection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataManualSelection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataManualSelection

% Last Modified by GUIDE v2.5 08-May-2025 16:24:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataManualSelection_OpeningFcn, ...
                   'gui_OutputFcn',  @DataManualSelection_OutputFcn, ...
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


% --- Executes just before DataManualSelection is made visible.
function DataManualSelection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataManualSelection (see VARARGIN)

% Choose default command line output for DataManualSelection
handles.output = hObject;

% Set initial values for all components
handles.Filename_txt.String=['Filename: ' varargin{1}.Filename];
handles.Wavelength_ed.String=num2str(varargin{1}.Wavelength);

handles.psi_panel.Title=char(968);
handles.tth_panel.Title=['2' char(952)];
handles.beta_panel.Title=char(946);

if varargin{1}.DataAvailability(1)==0
    handles.h_NotAvailable_rad.Value=1;
    handles.h_columnN_ed.Enable='off';
    handles.h_columnNplus_btn.Enable='off';
    handles.h_columnNminus_btn.Enable='off';
    handles.h_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(1)~=0
    handles.h_columnN_rad.Value=1;
    handles.h_columnN_ed.String=num2str(varargin{1}.DataColumnN(1));
    handles.h_columnN_ed.Enable='on';
    handles.h_columnNplus_btn.Enable='on';
    handles.h_columnNminus_btn.Enable='on';
    handles.h_formula_ed.Enable='off';
else
    handles.h_formula_rad.Value=1;
    handles.h_formula_ed.String=num2str(varargin{1}.DataFormula{1});
    handles.h_columnN_ed.Enable='off';
    handles.h_columnNplus_btn.Enable='off';
    handles.h_columnNminus_btn.Enable='off';
    handles.h_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(2)==0
    handles.k_NotAvailable_rad.Value=1;
    handles.k_columnN_ed.Enable='off';
    handles.k_columnNplus_btn.Enable='off';
    handles.k_columnNminus_btn.Enable='off';
    handles.k_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(2)~=0
    handles.k_columnN_rad.Value=1;
    handles.k_columnN_ed.String=num2str(varargin{1}.DataColumnN(2));
    handles.k_columnN_ed.Enable='on';
    handles.k_columnNplus_btn.Enable='on';
    handles.k_columnNminus_btn.Enable='on';
    handles.k_formula_ed.Enable='off';
else
    handles.k_formula_rad.Value=1;
    handles.k_formula_ed.String=num2str(varargin{1}.DataFormula{2});
    handles.k_columnN_ed.Enable='off';
    handles.k_columnNplus_btn.Enable='off';
    handles.k_columnNminus_btn.Enable='off';
    handles.k_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(3)==0
    handles.l_NotAvailable_rad.Value=1;
    handles.l_columnN_ed.Enable='off';
    handles.l_columnNplus_btn.Enable='off';
    handles.l_columnNminus_btn.Enable='off';
    handles.l_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(3)~=0
    handles.l_columnN_rad.Value=1;
    handles.l_columnN_ed.String=num2str(varargin{1}.DataColumnN(3));
    handles.l_columnN_ed.Enable='on';
    handles.l_columnNplus_btn.Enable='on';
    handles.l_columnNminus_btn.Enable='on';
    handles.l_formula_ed.Enable='off';
else
    handles.l_formula_rad.Value=1;
    handles.l_formula_ed.String=num2str(varargin{1}.DataFormula{3});
    handles.l_columnN_ed.Enable='off';
    handles.l_columnNplus_btn.Enable='off';
    handles.l_columnNminus_btn.Enable='off';
    handles.l_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(4)==0
    handles.psi_NotAvailable_rad.Value=1;
    handles.psi_columnN_ed.Enable='off';
    handles.psi_columnNplus_btn.Enable='off';
    handles.psi_columnNminus_btn.Enable='off';
    handles.psi_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(4)~=0
    handles.psi_columnN_rad.Value=1;
    handles.psi_columnN_ed.String=num2str(varargin{1}.DataColumnN(4));
    handles.psi_columnN_ed.Enable='on';
    handles.psi_columnNplus_btn.Enable='on';
    handles.psi_columnNminus_btn.Enable='on';
    handles.psi_formula_ed.Enable='off';
    handles.psi_Unit_txt.Enable='on';
    handles.psi_Unit_pop.Enable='on';
else
    handles.psi_formula_rad.Value=1;
    handles.psi_formula_ed.String=num2str(varargin{1}.DataFormula{4});
    handles.psi_columnN_ed.Enable='off';
    handles.psi_columnNplus_btn.Enable='off';
    handles.psi_columnNminus_btn.Enable='off';
    handles.psi_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(5)==0
    handles.phi_NotAvailable_rad.Value=1;
    handles.phi_columnN_ed.Enable='off';
    handles.phi_columnNplus_btn.Enable='off';
    handles.phi_columnNminus_btn.Enable='off';
    handles.phi_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(5)~=0
    handles.phi_columnN_rad.Value=1;
    handles.phi_columnN_ed.String=num2str(varargin{1}.DataColumnN(5));
    handles.phi_columnN_ed.Enable='on';
    handles.phi_columnNplus_btn.Enable='on';
    handles.phi_columnNminus_btn.Enable='on';
    handles.phi_formula_ed.Enable='off';
    handles.phi_Unit_txt.Enable='on';
    handles.phi_Unit_pop.Enable='on';
else
    handles.phi_formula_rad.Value=1;
    handles.phi_formula_ed.String=num2str(varargin{1}.DataFormula{5});
    handles.phi_columnN_ed.Enable='off';
    handles.phi_columnNplus_btn.Enable='off';
    handles.phi_columnNminus_btn.Enable='off';
    handles.phi_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(6)==0
    handles.tth_NotAvailable_rad.Value=1;
    handles.tth_columnN_ed.Enable='off';
    handles.tth_columnNplus_btn.Enable='off';
    handles.tth_columnNminus_btn.Enable='off';
    handles.tth_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(6)~=0
    handles.tth_columnN_rad.Value=1;
    handles.tth_columnN_ed.String=num2str(varargin{1}.DataColumnN(6));
    handles.tth_columnN_ed.Enable='on';
    handles.tth_columnNplus_btn.Enable='on';
    handles.tth_columnNminus_btn.Enable='on';
    handles.tth_formula_ed.Enable='off';
    handles.tth_Unit_txt.Enable='on';
    handles.tth_Unit_pop.Enable='on';
else
    handles.tth_formula_rad.Value=1;
    handles.tth_formula_ed.String=num2str(varargin{1}.DataFormula{6});
    handles.tth_columnN_ed.Enable='off';
    handles.tth_columnNplus_btn.Enable='off';
    handles.tth_columnNminus_btn.Enable='off';
    handles.tth_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(7)==0
    handles.fwhm_NotAvailable_rad.Value=1;
    handles.fwhm_columnN_ed.Enable='off';
    handles.fwhm.Enable='off';
    handles.fwhm_columnNplus_btn.Enable='off';
    handles.fwhm_columnNminus_btn.Enable='off';
    handles.fwhm_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(7)~=0
    handles.fwhm_columnN_rad.Value=1;
    handles.fwhm_columnN_ed.String=num2str(varargin{1}.DataColumnN(7));
    handles.fwhm_columnN_ed.Enable='on';
    handles.fwhm_columnNplus_btn.Enable='on';
    handles.fwhm_columnNminus_btn.Enable='on';
    handles.fwhm_formula_ed.Enable='off';
    handles.fwhm_Unit_txt.Enable='on';
    handles.fwhm_Unit_pop.Enable='on';
else
    handles.fwhm_formula_rad.Value=1;
    handles.fwhm_formula_ed.String=num2str(varargin{1}.DataFormula{7});
    handles.fwhm_columnN_ed.Enable='off';
    handles.fwhm_columnNplus_btn.Enable='off';
    handles.fwhm_columnNminus_btn.Enable='off';
    handles.fwhm_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(8)==0
    handles.beta_NotAvailable_rad.Value=1;
    handles.beta_columnN_ed.Enable='off';
    handles.beta.Enable='off';
    handles.beta_columnNplus_btn.Enable='off';
    handles.beta_columnNminus_btn.Enable='off';
    handles.beta_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(8)~=0
    handles.beta_columnN_rad.Value=1;
    handles.beta_columnN_ed.String=num2str(varargin{1}.DataColumnN(8));
    handles.beta_columnN_ed.Enable='on';
    handles.beta_columnNplus_btn.Enable='on';
    handles.beta_columnNminus_btn.Enable='on';
    handles.beta_formula_ed.Enable='off';
    handles.beta_Unit_txt.Enable='on';
    handles.beta_Unit_pop.Enable='on';
else
    handles.beta_formula_rad.Value=1;
    handles.beta_formula_ed.String=num2str(varargin{1}.DataFormula{8});
    handles.beta_columnN_ed.Enable='off';
    handles.beta_columnNplus_btn.Enable='off';
    handles.beta_columnNminus_btn.Enable='off';
    handles.beta_formula_ed.Enable='on';
end
if varargin{1}.DataAvailability(9)==0
    handles.label_NotAvailable_rad.Value=1;
    handles.label_columnN_ed.Enable='off';
    handles.label.Enable='off';
    handles.label_columnNplus_btn.Enable='off';
    handles.label_columnNminus_btn.Enable='off';
    handles.label_formula_ed.Enable='off';
elseif varargin{1}.DataColumnN(9)~=0
    handles.label_columnN_rad.Value=1;
    handles.label_columnN_ed.String=num2str(varargin{1}.DataColumnN(9));
    handles.label_columnN_ed.Enable='on';
    handles.label_columnNplus_btn.Enable='on';
    handles.label_columnNminus_btn.Enable='on';
    handles.label_formula_ed.Enable='off';
else
    handles.label_formula_rad.Value=1;
    handles.label_formula_ed.String=num2str(varargin{1}.DataFormula{9});
    handles.label_columnN_ed.Enable='off';
    handles.label_columnNplus_btn.Enable='off';
    handles.label_columnNminus_btn.Enable='off';
    handles.label_formula_ed.Enable='on';
end

% load the table inside uitable
InputTable=varargin{1}.InputTable;
ColumnNames=InputTable.Properties.VariableNames; % original
ColumnNames1=InputTable.Properties.VariableNames; % actual label shown in table
data=[];
for c=1:numel(ColumnNames)
    if iscell(InputTable.(ColumnNames{c}))
        data=[data InputTable.(ColumnNames{c})];
    else
        data=[data mat2cell(InputTable.(ColumnNames{c}),ones(size(InputTable.(ColumnNames{c}))))];
    end

    ColumnNames1{c}=['var' num2str(c) ': ' ColumnNames{c}];
end
set(handles.Data_tab,'Data',data,'ColumnName',ColumnNames1,'ColumnWidth','auto');
handles.InputTable=InputTable;
handles.ColumnNames=ColumnNames;
handles.InputDataFormula=varargin{1}.DataFormula;

% load units options
if isempty(varargin{1}.DataUnitsOpt{4})
    handles.psi_Unit_pop.Value=1;
else
    switch varargin{1}.DataUnitsOpt{4} % psi
        case 'uknown'
            handles.psi_Unit_pop.Value=1;
        case 'deg'
            handles.psi_Unit_pop.Value=2;
        case 'rad'
            handles.psi_Unit_pop.Value=3;
    end
end
if isempty(varargin{1}.DataUnitsOpt{5})
    handles.phi_Unit_pop.Value=1;
else
    switch varargin{1}.DataUnitsOpt{5} % phi
        case 'uknown'
            handles.phi_Unit_pop.Value=1;
        case 'deg'
            handles.phi_Unit_pop.Value=2;
        case 'rad'
            handles.phi_Unit_pop.Value=3;
    end
end
if isempty(varargin{1}.DataUnitsOpt{6})
    handles.tth_Unit_pop.Value=1;
else
    switch varargin{1}.DataUnitsOpt{6} % tth
        case 'uknown'
            handles.tth_Unit_pop.Value=1;
        case 'deg'
            handles.tth_Unit_pop.Value=2;
        case 'rad'
            handles.tth_Unit_pop.Value=3;
        case 'A-1'
            handles.tth_Unit_pop.Value=4;
    end
end
if isempty(varargin{1}.DataUnitsOpt{7})
    handles.fwhm_Unit_pop.Value=1;
else
    switch varargin{1}.DataUnitsOpt{7} % fwhm
        case 'uknown'
            handles.fwhm_Unit_pop.Value=1;
        case 'deg'
            handles.fwhm_Unit_pop.Value=2;
        case 'rad'
            handles.fwhm_Unit_pop.Value=3;
        case 'A-1'
            handles.fwhm_Unit_pop.Value=4;
    end
end
if isempty(varargin{1}.DataUnitsOpt{8})
    handles.beta_Unit_pop.Value=1;
else
    switch varargin{1}.DataUnitsOpt{8} % beta
        case 'uknown'
            handles.beta_Unit_pop.Value=1;
        case 'deg'
            handles.beta_Unit_pop.Value=2;
        case 'rad'
            handles.beta_Unit_pop.Value=3;
        case 'A-1'
            handles.beta_Unit_pop.Value=4;
    end
end

% Default Output is the Input
handles.OutputData=varargin{1}.InputData;
handles.Wavelength=varargin{1}.Wavelength;
handles.DataFormula=varargin{1}.DataFormula;
handles.DataColumnN=varargin{1}.DataColumnN;
handles.DataAvailability=varargin{1}.DataAvailability;
handles.DataUnitsOpt=varargin{1}.DataUnitsOpt;

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataManualSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataManualSelection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.OutputData;
varargout{2} = handles.Wavelength;
varargout{3} = handles.DataAvailability;
varargout{4} = handles.DataColumnN;
varargout{5} = handles.DataFormula;
varargout{6} = handles.DataUnitsOpt;

% The GUI can be deleted now
delete(handles.figure1);


% --- Executes on button press in Apply_btn.
function Apply_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script will go through every panel and will try to select the data, or
% 'process' the data with the formula. Direct selection of columns should
% be safe, because the edit function of edit boxes will fix numbers to
% integer 1<=n<=ColumnNumber. If there will be a mistake in any formula,
% errordialog is created. If everything is ok, new dataset is export to the
% main GUI.
% Check wavelength
handles.Wavelength=str2num(handles.Wavelength_ed.String);
if isnan(handles.Wavelength)
    errordlg('Wavelength should be a real number greater than zero!','Wrong input');
else
    if handles.Wavelength<=0
        errordlg('Wavelength should be a real number greater than zero!','Wrong input');
    end
end

% Check, if all units are known
isserUnits=0;
handles.DataUnitsOpt{4}=handles.psi_Unit_pop.String{handles.psi_Unit_pop.Value};
if handles.psi_columnN_rad.Value==1 && strcmp(handles.psi_Unit_pop.String{handles.psi_Unit_pop.Value},'uknown')
    isserUnits=1;
end
handles.DataUnitsOpt{5}=handles.phi_Unit_pop.String{handles.phi_Unit_pop.Value};
if handles.phi_columnN_rad.Value==1 && strcmp(handles.phi_Unit_pop.String{handles.phi_Unit_pop.Value},'uknown')
    isserUnits=1;
end
handles.DataUnitsOpt{6}=handles.tth_Unit_pop.String{handles.tth_Unit_pop.Value};
if handles.tth_columnN_rad.Value==1 && strcmp(handles.tth_Unit_pop.String{handles.tth_Unit_pop.Value},'uknown')
    isserUnits=1;
end
handles.DataUnitsOpt{7}=handles.fwhm_Unit_pop.String{handles.fwhm_Unit_pop.Value};
if handles.fwhm_columnN_rad.Value==1 && strcmp(handles.fwhm_Unit_pop.String{handles.fwhm_Unit_pop.Value},'uknown')
    isserUnits=1;
end
handles.DataUnitsOpt{8}=handles.beta_Unit_pop.String{handles.beta_Unit_pop.Value};
if handles.beta_columnN_rad.Value==1 && strcmp(handles.beta_Unit_pop.String{handles.beta_Unit_pop.Value},'uknown')
    isserUnits=1;
end
if isserUnits==1
    errordlg('Units of all data should be known!','Wrong input');
end
guidata(hObject, handles);

handles.DataFormula=handles.InputDataFormula;
iserr=zeros(1,8);

[handles,iserr]=GetColumnOrFormula('h',handles,iserr);
[handles,iserr]=GetColumnOrFormula('k',handles,iserr);
[handles,iserr]=GetColumnOrFormula('l',handles,iserr);
[handles,iserr]=GetColumnOrFormula('psi',handles,iserr);
[handles,iserr]=GetColumnOrFormula('phi',handles,iserr);
[handles,iserr]=GetColumnOrFormula('tth',handles,iserr);
[handles,iserr]=GetColumnOrFormula('fwhm',handles,iserr);
[handles,iserr]=GetColumnOrFormula('beta',handles,iserr);
[handles,iserr]=GetColumnOrFormula('label',handles,iserr);

if sum(iserr)==0
    guidata(hObject, handles);
    uiresume;
end


function label_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to label_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of label_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function label_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in label_columnNplus_btn.
function label_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to label_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.label_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.label_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


% --- Executes on button press in label_columnNminus_btn.
function label_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to label_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.label_columnN_ed.String);
if n-1>0
    handles.label_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);



function label_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to label_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of label_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function label_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function beta_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to beta_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of beta_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function beta_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in beta_columnNplus_btn.
function beta_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to beta_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.beta_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.beta_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


% --- Executes on button press in beta_columnNminus_btn.
function beta_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to beta_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.beta_columnN_ed.String);
if n-1>0
    handles.beta_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);



function beta_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to beta_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of beta_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function beta_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fwhm_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fwhm_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of fwhm_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function fwhm_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fwhm_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fwhm_columnNplus_btn.
function fwhm_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.fwhm_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.fwhm_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


% --- Executes on button press in fwhm_columnNminus_btn.
function fwhm_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.fwhm_columnN_ed.String);
if n-1>0
    handles.fwhm_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);



function fwhm_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fwhm_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of fwhm_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function fwhm_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fwhm_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function tth_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to tth_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tth_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of tth_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function tth_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tth_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tth_columnNminus_btn.
function tth_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to tth_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.tth_columnN_ed.String);
if n-1>0
    handles.tth_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


% --- Executes on button press in tth_columnNplus_btn.
function tth_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to tth_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.tth_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.tth_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);



function tth_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to tth_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tth_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of tth_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function tth_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tth_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function psi_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to psi_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psi_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of psi_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function psi_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psi_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in psi_columnNminus_btn.
function psi_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to psi_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.psi_columnN_ed.String);
if n-1>0
    handles.psi_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


% --- Executes on button press in psi_columnNplus_btn.
function psi_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to psi_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.psi_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.psi_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);



function psi_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to psi_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psi_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of psi_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function psi_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psi_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function l_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to l_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of l_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of l_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function l_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in l_columnNminus_btn.
function l_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to l_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.l_columnN_ed.String);
if n-1>0
    handles.l_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


% --- Executes on button press in l_columnNplus_btn.
function l_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to l_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.l_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.l_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);



function l_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to l_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of l_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of l_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function l_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to k_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of k_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function k_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k_columnNminus_btn.
function k_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to k_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.k_columnN_ed.String);
if n-1>0
    handles.k_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


% --- Executes on button press in k_columnNplus_btn.
function k_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to k_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.k_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.k_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


function k_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to k_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of k_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function k_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to h_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of h_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function h_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in h_columnNplus_btn.
function h_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to h_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.h_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.h_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


% --- Executes on button press in h_columnNminus_btn.
function h_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to h_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.h_columnN_ed.String);
if n-1>0
    handles.h_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


function h_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to h_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of h_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function h_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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


% --- Executes on button press in Help_btn.
function Help_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Help_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg({'In this dialog window you can choose, in which columns are your data located in the file/table you have loaded. You have three options: 1) say that the data are not available, 2) specify the columns where the data are, 3) specify the exact formula, how to compute your data points from the available columns.',...
'','Concerning the third option, the formula should be written in Matlab syntax. The names of the variables representing columns should be ''var#'' (# is a number). These names could be seen in the table on the left together with the column names automaticaly detected by Matlab. However, be catious - the automaticaly detected names might not look very nice, because Matlab replaced probably every character excluding letter and number with an underscore character.',...
'','During entering the formula you can also call the wavelength variable by the expression ''lam''.'},...
'Help for manual column selection');


% --- Executes on button press in h_NotAvailable_rad.
function h_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to h_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of h_NotAvailable_rad
handles.h_columnN_ed.Enable='off';
handles.h_columnNplus_btn.Enable='off';
handles.h_columnNminus_btn.Enable='off';
handles.h_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in h_columnN_rad.
function h_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to h_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of h_columnN_rad
handles.h_columnN_ed.Enable='on';
handles.h_columnNplus_btn.Enable='on';
handles.h_columnNminus_btn.Enable='on';
handles.h_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in h_formula_rad.
function h_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to h_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of h_formula_rad
handles.h_columnN_ed.Enable='off';
handles.h_columnNplus_btn.Enable='off';
handles.h_columnNminus_btn.Enable='off';
handles.h_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in k_NotAvailable_rad.
function k_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to k_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k_NotAvailable_rad
handles.k_columnN_ed.Enable='off';
handles.k_columnNplus_btn.Enable='off';
handles.k_columnNminus_btn.Enable='off';
handles.k_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in k_columnN_rad.
function k_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to k_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k_columnN_rad
handles.k_columnN_ed.Enable='on';
handles.k_columnNplus_btn.Enable='on';
handles.k_columnNminus_btn.Enable='on';
handles.k_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in k_formula_rad.
function k_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to k_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k_formula_rad
handles.k_columnN_ed.Enable='off';
handles.k_columnNplus_btn.Enable='off';
handles.k_columnNminus_btn.Enable='off';
handles.k_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in l_NotAvailable_rad.
function l_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to l_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l_NotAvailable_rad
handles.l_columnN_ed.Enable='off';
handles.l_columnNplus_btn.Enable='off';
handles.l_columnNminus_btn.Enable='off';
handles.l_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in l_columnN_rad.
function l_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to l_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l_columnN_rad
handles.l_columnN_ed.Enable='on';
handles.l_columnNplus_btn.Enable='on';
handles.l_columnNminus_btn.Enable='on';
handles.l_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in l_formula_rad.
function l_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to l_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of l_formula_rad
handles.l_columnN_ed.Enable='off';
handles.l_columnNplus_btn.Enable='off';
handles.l_columnNminus_btn.Enable='off';
handles.l_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in psi_NotAvailable_rad.
function psi_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to psi_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of psi_NotAvailable_rad
handles.psi_columnN_ed.Enable='off';
handles.psi_columnNplus_btn.Enable='off';
handles.psi_columnNminus_btn.Enable='off';
handles.psi_Unit_txt.Enable='off';
handles.psi_Unit_pop.Enable='off';
handles.psi_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in psi_columnN_rad.
function psi_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to psi_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of psi_columnN_rad
handles.psi_columnN_ed.Enable='on';
handles.psi_columnNplus_btn.Enable='on';
handles.psi_columnNminus_btn.Enable='on';
handles.psi_Unit_txt.Enable='on';
handles.psi_Unit_pop.Enable='on';
handles.psi_columnNminus_btn.Enable='on';
handles.psi_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in psi_formula_rad.
function psi_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to psi_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of psi_formula_rad
handles.psi_columnN_ed.Enable='off';
handles.psi_columnNplus_btn.Enable='off';
handles.psi_columnNminus_btn.Enable='off';
handles.psi_Unit_txt.Enable='off';
handles.psi_Unit_pop.Enable='off';
handles.psi_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in tth_NotAvailable_rad.
function tth_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to tth_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tth_NotAvailable_rad
handles.tth_columnN_ed.Enable='off';
handles.tth_columnNplus_btn.Enable='off';
handles.tth_columnNminus_btn.Enable='off';
handles.tth_Unit_txt.Enable='off';
handles.tth_Unit_pop.Enable='off';
handles.tth_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in tth_columnN_rad.
function tth_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to tth_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tth_columnN_rad
handles.tth_columnN_ed.Enable='on';
handles.tth_columnNplus_btn.Enable='on';
handles.tth_columnNminus_btn.Enable='on';
handles.tth_Unit_txt.Enable='on';
handles.tth_Unit_pop.Enable='on';
handles.tth_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in tth_formula_rad.
function tth_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to tth_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tth_formula_rad
handles.tth_columnN_ed.Enable='off';
handles.tth_columnNplus_btn.Enable='off';
handles.tth_columnNminus_btn.Enable='off';
handles.tth_Unit_txt.Enable='off';
handles.tth_Unit_pop.Enable='off';
handles.tth_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in fwhm_NotAvailable_rad.
function fwhm_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fwhm_NotAvailable_rad
handles.fwhm_columnN_ed.Enable='off';
handles.fwhm_columnNplus_btn.Enable='off';
handles.fwhm_columnNminus_btn.Enable='off';
handles.fwhm_Unit_txt.Enable='off';
handles.fwhm_Unit_pop.Enable='off';
handles.fwhm_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in fwhm_columnN_rad.
function fwhm_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fwhm_columnN_rad
handles.fwhm_columnN_ed.Enable='on';
handles.fwhm_columnNplus_btn.Enable='on';
handles.fwhm_columnNminus_btn.Enable='on';
handles.fwhm_Unit_txt.Enable='on';
handles.fwhm_Unit_pop.Enable='on';
handles.fwhm_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in fwhm_formula_rad.
function fwhm_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fwhm_formula_rad
handles.fwhm_columnN_ed.Enable='off';
handles.fwhm_columnNplus_btn.Enable='off';
handles.fwhm_columnNminus_btn.Enable='off';
handles.fwhm_Unit_txt.Enable='off';
handles.fwhm_Unit_pop.Enable='off';
handles.fwhm_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in beta_NotAvailable_rad.
function beta_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to beta_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beta_NotAvailable_rad
handles.beta_columnN_ed.Enable='off';
handles.beta_columnNplus_btn.Enable='off';
handles.beta_columnNminus_btn.Enable='off';
handles.beta_Unit_txt.Enable='off';
handles.beta_Unit_pop.Enable='off';
handles.beta_formula_ed.Enable='off';
guidata(hObject, handles);


% % --- Executes during object creation, after setting all properties.
% function beta_columnN_rad_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to beta_columnN_rad (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in beta_columnN_rad.
function beta_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to beta_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beta_columnN_rad
handles.beta_columnN_ed.Enable='on';
handles.beta_columnNplus_btn.Enable='on';
handles.beta_columnNminus_btn.Enable='on';
handles.beta_Unit_txt.Enable='on';
handles.beta_Unit_pop.Enable='on';
handles.beta_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in beta_formula_rad.
function beta_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to beta_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beta_formula_rad
handles.beta_columnN_ed.Enable='off';
handles.beta_columnNplus_btn.Enable='off';
handles.beta_columnNminus_btn.Enable='off';
handles.beta_Unit_txt.Enable='off';
handles.beta_Unit_pop.Enable='off';
handles.beta_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in label_NotAvailable_rad.
function label_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to label_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of label_NotAvailable_rad
handles.label_columnN_ed.Enable='off';
handles.label_columnNplus_btn.Enable='off';
handles.label_columnNminus_btn.Enable='off';
handles.label_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in label_columnN_rad.
function label_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to label_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of label_columnN_rad
handles.label_columnN_ed.Enable='on';
handles.label_columnNplus_btn.Enable='on';
handles.label_columnNminus_btn.Enable='on';
handles.label_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in label_formula_rad.
function label_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to label_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of label_formula_rad
handles.label_columnN_ed.Enable='off';
handles.label_columnNplus_btn.Enable='off';
handles.label_columnNminus_btn.Enable='off';
handles.label_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on button press in Cancel_btn.
function Cancel_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Here, 'uiresume' should be sufficient, because Input data are in Output
% as default
uiresume;


function [handles,iserr]=GetColumnOrFormula(var,handles,iserr)
N=numel(handles.InputTable.(handles.ColumnNames{1}));
indVar=find(strcmp(var,{'h', 'k', 'l', 'psi', 'phi', 'tth', 'fwhm', 'beta', 'label'}));
lam=handles.Wavelength;
%%%%% var = 'h', 'k', 'l', 'psi', 'phi', 'tth', 'fwhm', 'beta', 'label' %%%%%
if handles.([var '_NotAvailable_rad']).Value==1
    if indVar==9
        handles.OutputData.(var)=cell(N,1);
    else
        handles.OutputData.(var)=NaN*ones(N,1);
    end
    handles.DataAvailability(indVar)=0;
    handles.DataColumnN(indVar)=0;
elseif handles.([var '_columnN_rad']).Value==1
    ind=str2num(handles.([var '_columnN_ed']).String);
    %%% Here, the unit option should be correctly treated and recompute
    %%% everything into deg (this unit is used for computation and fit)
    hlp=handles.InputTable.(handles.ColumnNames{ind});
    switch var
        case 'psi'
            switch handles.psi_Unit_pop.String{handles.psi_Unit_pop.Value}
                case 'deg'
                    handles.OutputData.psi=hlp;
                case 'rad'
                    handles.OutputData.psi=hlp*180/pi;
            end
        case 'phi'
            switch handles.phi_Unit_pop.String{handles.phi_Unit_pop.Value}
                case 'deg'
                    handles.OutputData.phi=hlp;
                case 'rad'
                    handles.OutputData.phi=hlp*180/pi;
            end
        case 'tth'
            switch handles.tth_Unit_pop.String{handles.tth_Unit_pop.Value}
                case 'deg'
                    handles.OutputData.tth=hlp;
                case 'rad'
                    handles.OutputData.tth=hlp*180/pi;
                case 'A-1'
                    handles.OutputData.tth=2*asind(hlp*lam/(4*pi));
            end
        case 'fwhm'
            switch handles.fwhm_Unit_pop.String{handles.fwhm_Unit_pop.Value}
               case 'deg'
                   handles.OutputData.fwhm=hlp;
               case 'rad'
                   handles.OutputData.fwhm=hlp*180/pi;
               case 'A-1'
                   handles.OutputData.fwhm=2*hlp*lam/(4*pi)./cosd(handles.OutputData.tth/2); % tth should be already known
            end
        case 'beta'
            switch handles.beta_Unit_pop.String{handles.beta_Unit_pop.Value}
               case 'deg'
                   handles.OutputData.beta=hlp;
               case 'rad'
                   handles.OutputData.beta=hlp*180/pi;
               case 'A-1'
                   handles.OutputData.beta=2*hlp*lam/(4*pi)./cosd(handles.OutputData.tth/2); % tth should be already known
            end
        otherwise
            handles.OutputData.(var)=handles.InputTable.(handles.ColumnNames{ind});
    end
    handles.DataColumnN(indVar)=ind;
    handles.DataAvailability(indVar)=1;
else
    formula=handles.([var '_formula_ed']).String;
    if contains(formula,'formula')
        iserr(indVar)=1;
    else
        % I will create variables var# from
        % handles.InputTable.(ColumnNames{c}). Then, if the formula contain
        % var1, var2, ..., it can be directly evaluated
        for c=1:numel(handles.ColumnNames)
            eval(['var' num2str(c) '=handles.InputTable.' handles.ColumnNames{c} ';']);
        end

        try
            handles.OutputData.(var)=eval([formula ';']);
            handles.DataFormula{indVar}=formula; % I must save the original formula, or the Manual Selection dialog will not work properly on the next run
            handles.DataColumnN(indVar)=0;
            handles.DataAvailability(indVar)=1;

            if indVar==8
                % label should be string
                if ~(ischar(handles.OutputData.(var)) | isstring(handles.OutputData.(var)))
                    try 
                        handles.OutputData.(var)=num2str(handles.OutputData.(var));
                    catch
                        iserr(indVar)=1;
                    end
                end
            else
                if ~isnumeric(handles.OutputData.(var))
                    handles.OutputData.(var)=NaN*ones(N,1);
                    iserr(indVar)=1;
                end
            end
        catch
            iserr(indVar)=1;
        end
    end
    % If there were some error, show error dialog
    if iserr(indVar)
        errordlg(['The formula for the parameter ''' var ''' could not be evaluated. Please, check the syntax.'],'Error during formula evaluation');
    end
    % Finally, one have to test, if the resulting array/matrix has the
    % correct size. (missing dot in the formula could cause troubles)
    sz=size(handles.OutputData.(var));
    if ~(sz(1)==N && sz(2)==1)
        iserr(indVar)=1;
        errordlg(['The size of the resulting ''' var ''' variable is wrong. Check the formula. When working with matrices/arrays, the problem is caused often by a missing dot.'],'Error during formula evaluation');
    end
end


% --- Executes on button press in phi_columnNminus_btn.
function phi_columnNminus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to phi_columnNminus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.phi_columnN_ed.String);
if n-1>0
    handles.phi_columnN_ed.String=num2str(n-1);
end
guidata(hObject,handles);


% --- Executes on button press in phi_columnNplus_btn.
function phi_columnNplus_btn_Callback(hObject, eventdata, handles)
% hObject    handle to phi_columnNplus_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=str2num(handles.phi_columnN_ed.String);
if n+1<=numel(handles.ColumnNames)
    handles.phi_columnN_ed.String=num2str(n+1);
end
guidata(hObject,handles);


function phi_columnN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to phi_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi_columnN_ed as text
%        str2double(get(hObject,'String')) returns contents of phi_columnN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if n~=round(n)
        iserr=1;
    else
        if n<1 | n>numel(handles.ColumnNames)
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg(['Column number should an integer from 1 to ' num2str(numel(handles.ColumnNames)) '.']);
    hObject.String='1';
end

% --- Executes during object creation, after setting all properties.
function phi_columnN_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi_columnN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phi_formula_ed_Callback(hObject, eventdata, handles)
% hObject    handle to phi_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi_formula_ed as text
%        str2double(get(hObject,'String')) returns contents of phi_formula_ed as a double


% --- Executes during object creation, after setting all properties.
function phi_formula_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi_formula_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in phi_NotAvailable_rad.
function phi_NotAvailable_rad_Callback(hObject, eventdata, handles)
% hObject    handle to phi_NotAvailable_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phi_NotAvailable_rad
handles.phi_columnN_ed.Enable='off';
handles.phi_columnNplus_btn.Enable='off';
handles.phi_columnNminus_btn.Enable='off';
handles.phi_Unit_txt.Enable='off';
handles.phi_Unit_pop.Enable='off';
handles.phi_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in phi_columnN_rad.
function phi_columnN_rad_Callback(hObject, eventdata, handles)
% hObject    handle to phi_columnN_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phi_columnN_rad
handles.phi_columnN_ed.Enable='on';
handles.phi_columnNplus_btn.Enable='on';
handles.phi_columnNminus_btn.Enable='on';
handles.phi_Unit_txt.Enable='on';
handles.phi_Unit_pop.Enable='on';
handles.phi_formula_ed.Enable='off';
guidata(hObject, handles);


% --- Executes on button press in phi_formula_rad.
function phi_formula_rad_Callback(hObject, eventdata, handles)
% hObject    handle to phi_formula_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phi_formula_rad
handles.phi_columnN_ed.Enable='off';
handles.phi_columnNplus_btn.Enable='off';
handles.phi_columnNminus_btn.Enable='off';
handles.phi_Unit_txt.Enable='off';
handles.phi_Unit_pop.Enable='off';
handles.phi_formula_ed.Enable='on';
guidata(hObject, handles);


% --- Executes on selection change in tth_Unit_pop.
function tth_Unit_pop_Callback(hObject, eventdata, handles)
% hObject    handle to tth_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tth_Unit_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tth_Unit_pop


% --- Executes during object creation, after setting all properties.
function tth_Unit_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tth_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in phi_Unit_pop.
function phi_Unit_pop_Callback(hObject, eventdata, handles)
% hObject    handle to phi_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns phi_Unit_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from phi_Unit_pop


% --- Executes during object creation, after setting all properties.
function phi_Unit_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in beta_Unit_pop.
function beta_Unit_pop_Callback(hObject, eventdata, handles)
% hObject    handle to beta_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns beta_Unit_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from beta_Unit_pop


% --- Executes during object creation, after setting all properties.
function beta_Unit_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fwhm_Unit_pop.
function fwhm_Unit_pop_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fwhm_Unit_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fwhm_Unit_pop


% --- Executes during object creation, after setting all properties.
function fwhm_Unit_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fwhm_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in psi_Unit_pop.
function psi_Unit_pop_Callback(hObject, eventdata, handles)
% hObject    handle to psi_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns psi_Unit_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from psi_Unit_pop


% --- Executes during object creation, after setting all properties.
function psi_Unit_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psi_Unit_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wavelength_ed_Callback(hObject, eventdata, handles)
% hObject    handle to Wavelength_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wavelength_ed as text
%        str2double(get(hObject,'String')) returns contents of Wavelength_ed as a double


% --- Executes during object creation, after setting all properties.
function Wavelength_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wavelength_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Wavelength_pop.
function Wavelength_pop_Callback(hObject, eventdata, handles)
% hObject    handle to Wavelength_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Wavelength_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Wavelength_pop
switch hObject.String{hObject.Value}
    case 'CuKalpha'
        handles.Wavelength_ed.String='1.5406';
    case 'CoKalpha'
        handles.Wavelength_ed.String='1.7902';
    case 'MoKalpha'
        handles.Wavelength_ed.String='0.7107';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Wavelength_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wavelength_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
