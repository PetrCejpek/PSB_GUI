function varargout = AdvancedShiftOpt_GUI(varargin)
% ADVANCEDSHIFTOPT_GUI MATLAB code for AdvancedShiftOpt_GUI.fig
%      ADVANCEDSHIFTOPT_GUI, by itself, creates a new ADVANCEDSHIFTOPT_GUI or raises the existing
%      singleton*.
%
%      H = ADVANCEDSHIFTOPT_GUI returns the handle to a new ADVANCEDSHIFTOPT_GUI or the handle to
%      the existing singleton*.
%
%      ADVANCEDSHIFTOPT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDSHIFTOPT_GUI.M with the given input arguments.
%
%      ADVANCEDSHIFTOPT_GUI('Property','Value',...) creates a new ADVANCEDSHIFTOPT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdvancedShiftOpt_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdvancedShiftOpt_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdvancedShiftOpt_GUI

% Last Modified by GUIDE v2.5 06-May-2025 16:44:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdvancedShiftOpt_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AdvancedShiftOpt_GUI_OutputFcn, ...
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


% --- Executes just before AdvancedShiftOpt_GUI is made visible.
function AdvancedShiftOpt_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AdvancedShiftOpt_GUI (see VARARGIN)

% Fill the form elements
% Material info
switch varargin{1}.MyOptions.MaterialInfo
    case 'SpecificOrientation'
        handles.SpecificOrientation_rad.Value=1;

        handles.StressFactorsModel_pop.Enable='off';
        handles.IdealPolycrystal_rad.Enable='off';
        handles.PolycrystalWithTexture_rad.Enable='off';
        handles.LoadODFfile_btn.Enable='off';
        handles.ODFfilename_ed.Enable='off';
    case 'XEC'
        handles.XEC_rad.Value=1;

        handles.StressFactorsModel_pop.Enable='off';
        handles.IdealPolycrystal_rad.Enable='off';
        handles.PolycrystalWithTexture_rad.Enable='off';
        handles.LoadODFfile_btn.Enable='off';
        handles.ODFfilename_ed.Enable='off';

        handles.Orientation_tab.Enable='off';
        handles.IndividualDatagroupOrientation_chbx.Enable='off';
    case 'StressFactors'
        handles.StressFactorsFromODF_rad.Value=1;
        
        handles.Orientation_tab.Enable='off';
        handles.IndividualDatagroupOrientation_chbx.Enable='off';
end
handles.KroenerCoeff_ed.String=num2str(varargin{1}.MyOptions.KroenerCoeff);
switch varargin{1}.MyOptions.XECmodel
    case 'Voigt'
        handles.XECmodel_pop.Value=1;
        handles.KroenerCoeff_txt.Enable='off';
        handles.KroenerCoeff_ed.Enable='off';
    case 'Reuss'
        handles.XECmodel_pop.Value=2;
        handles.KroenerCoeff_txt.Enable='off';
        handles.KroenerCoeff_ed.Enable='off';
    case 'Kroener'
        handles.XECmodel_pop.Value=3;
        handles.KroenerCoeff_txt.Enable='on';
        handles.KroenerCoeff_ed.Enable='on';
end

% Number of data groups
handles.DataGroupN=varargin{1}.DataGroupN;
% Orientation angles
DATA_or=varargin{1}.MyOptions.OrientationAngles;
if isempty(DATA_or)
    DATA={0 0 0 'z' 'x' 'z' 'extrinsic'};
    handles.Orientation_tab.Data=DATA;
else
    if handles.DataGroupN~=size(DATA_or,1)
        DATA=cell(handles.DataGroupN,7);
        for n=1:handles.DataGroupN
            DATA(n,:)=DATA_or(1,:);
        end
        handles.Orientation_tab.Data=DATA;
    else
        handles.Orientation_tab.Data=DATA_or;
    end
end
% Orientation of the individual datagroups
handles.IndividualDatagroupOrientation_chbx.Value=varargin{1}.MyOptions.DatagroupsIndividual;
% Computation method for the stress factors
switch varargin{1}.MyOptions.StressFactorsMethod
    case 'Reuss'
        handles.StressFactorsModel_pop.Value=1;
    case 'Modified Voigt'
        handles.StressFactorsModel_pop.Value=2;
    case 'Kroener'
        handles.StressFactorsModel_pop.Value=3;
    case 'Inverse Kroener'
        handles.StressFactorsModel_pop.Value=4;
end
% ODF type
switch varargin{1}.MyOptions.ODFtype
    case 'IdealPolycrystal'
        handles.IdealPolycrystal_rad.Value=1;
        
        handles.LoadODFfile_btn.Enable='off';
        handles.ODFfilename_ed.Enable='off';
    case 'PolycrystalWithTexture'
        handles.PolycrystalWithTexture_rad.Value=1;
end
% ODF filename
handles.ODFfilename_ed.String=varargin{1}.MyOptions.ODFfilename;
% Stress component to fit
handles.StressComponents_tab.Data=varargin{1}.MyOptions.StressComponents;
% Number of fitted parameters
handles.FitN_ed.String=num2str(varargin{1}.MyOptions.FitN);
% Measurement geometry
switch varargin{1}.MyOptions.Geometry
    case 'GAXRD'
        handles.MeasurementGeometry_pop.Value=1;
    case 'Pole figure'
        handles.MeasurementGeometry_pop.Value=2;
end

% Special characters in tables etc.
handles.Orientation_tab.ColumnName={[char(934) char(8321) '(' char(176) ')'],[char(936) '(' char(176) ')'],[char(934) char(8322) '(' char(176) ')'],'Axes 1','Axes 2','Axes 3','Rotation type'};
handles.Orientation_tab.ColumnEditable=true(1,7);
handles.StressComponents_tab.RowName={[char(963) char(8321) char(8321)];[char(963) char(8322) char(8322)];[char(963) char(8323) char(8323)];[char(963) char(8323) char(8322)];[char(963) char(8323) char(8321)];[char(963) char(8322) char(8321)]};

% Default Output is the Input
handles.OutputData=varargin{1}.MyOptions;

% Make the GUI modal
set(gcf,'WindowStyle','modal')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataManualSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AdvancedShiftOpt_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutputData;

% The GUI can be deleted now
delete(handles.figure1);

% --- Executes on button press in Help_btn.
function Help_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Help_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg({'In the table you can input formulae for individual stress components you want to fit. In the syntax, you should use components of vector p=[p(1); ... p(N)].';...
    '';...
    'For example, if you want to fit all of 6 stress parameters, you will put p(1), p(2) ... p(6) into corresponding rows.';...
    '';...
    ['If you want to fit only isotropic "in-plane" residual stress, you will set ' char(963) char(8321) char(8321) '=p(1), ' char(963) char(8322) char(8322) '=p(1), ' char(963) char(8323) char(8323) '=0, ' char(963) char(8323) char(8322) '=0, ' char(963) char(8323) char(8321) '=0, ' char(963) char(8322) char(8321) '=0.'];...
    },...
    'Help for formulae input');


% --- Executes on selection change in StressFactorsModel_pop.
function StressFactorsModel_pop_Callback(hObject, eventdata, handles)
% hObject    handle to StressFactorsModel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StressFactorsModel_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StressFactorsModel_pop


% --- Executes during object creation, after setting all properties.
function StressFactorsModel_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StressFactorsModel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in IndividualDatagroupOrientation_chbx.
function IndividualDatagroupOrientation_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to IndividualDatagroupOrientation_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IndividualDatagroupOrientation_chbx


% --- Executes on button press in LoadODFfile_btn.
function LoadODFfile_btn_Callback(hObject, eventdata, handles)
% hObject    handle to LoadODFfile_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastFilePath')
    [FileName,PathName,FilterIndex] = uigetfile({'*.*',  'All Files (*.*)'},'Select File to Open',handles.LastFilePath);
    handles.LastFilePath=PathName;
else
    [FileName,PathName,FilterIndex] = uigetfile({'*.*',  'All Files (*.*)'},'Select File to Open');
    handles.LastFilePath=PathName;
end
if FileName~=0
    handles.ODFfilename_ed.String=[PathName FileName];
end


function ODFfilename_ed_Callback(hObject, eventdata, handles)
% hObject    handle to ODFfilename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ODFfilename_ed as text
%        str2double(get(hObject,'String')) returns contents of ODFfilename_ed as a double


% --- Executes during object creation, after setting all properties.
function ODFfilename_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ODFfilename_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Apply_btn.
function Apply_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iserr=zeros(1,3);

[handles,iserr]=GetOptions(handles,iserr);

if sum(iserr)==0
    guidata(hObject, handles);
    uiresume;
end


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


% --- Executes on button press in SpecificOrientation_rad.
function SpecificOrientation_rad_Callback(hObject, eventdata, handles)
% hObject    handle to SpecificOrientation_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SpecificOrientation_rad
handles.Orientation_tab.Enable='on';
handles.IndividualDatagroupOrientation_chbx.Enable='on';

handles.StressFactorsModel_pop.Enable='off';
handles.IdealPolycrystal_rad.Enable='off';
handles.PolycrystalWithTexture_rad.Enable='off';
handles.LoadODFfile_btn.Enable='off';
handles.ODFfilename_ed.Enable='off';

guidata(hObject, handles);


% --- Executes on button press in StressFactorsFromODF_rad.
function StressFactorsFromODF_rad_Callback(hObject, eventdata, handles)
% hObject    handle to StressFactorsFromODF_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StressFactorsFromODF_rad
handles.Orientation_tab.Enable='off';
handles.IndividualDatagroupOrientation_chbx.Enable='off';

handles.StressFactorsModel_pop.Enable='on';
handles.IdealPolycrystal_rad.Enable='on';
handles.PolycrystalWithTexture_rad.Enable='on';
if handles.IdealPolycrystal_rad.Value==1
    handles.LoadODFfile_btn.Enable='off';
    handles.ODFfilename_ed.Enable='off';
else
    handles.LoadODFfile_btn.Enable='on';
    handles.ODFfilename_ed.Enable='on';
end

guidata(hObject, handles);

function [handles,iserr]=GetOptions(handles,iserr)
if handles.SpecificOrientation_rad.Value==1
    handles.OutputData.MaterialInfo='SpecificOrientation';
elseif handles.XEC_rad.Value==1
    handles.OutputData.MaterialInfo='XEC';
else
    handles.OutputData.MaterialInfo='StressFactors';
end

handles.OutputData.KroenerCoeff=str2num(handles.KroenerCoeff_ed.String);
handles.OutputData.XECmodel=handles.XECmodel_pop.String{handles.XECmodel_pop.Value};

% this should be ok, because MatLab should allow only numeric values in the
% table
handles.OutputData.OrientationAngles=handles.Orientation_tab.Data;

handles.OutputData.DatagroupsIndividual=handles.IndividualDatagroupOrientation_chbx.Value;

handles.OutputData.StressFactorsMethod=handles.StressFactorsModel_pop.String{handles.StressFactorsModel_pop.Value};

if handles.IdealPolycrystal_rad.Value==1
    handles.OutputData.ODFtype='IdealPolycrystal';
else
    handles.OutputData.ODFtype='PolycrystalWithTexture';
end

if handles.PolycrystalWithTexture_rad.Value==1
    try
        data=importdata(handles.ODFfilename_ed.String);
        if size(data.data,1)<10 | size(data.data,2)~=4
            iserr(1)=1;
        else
            handles.OutputData.ODFfilename=handles.ODFfilename_ed.String;
        end
    catch
        iserr(1)=1;
    end
    if iserr(1)
        errordlg('The data inside the provided ODF file do not have the desired format. They should contain the scattered data in 4 columns (3 orientation angles and 1 ODF value).','Wrong file format');
    end
end

FitN=str2num(handles.FitN_ed.String);
handles.OutputData.FitN=FitN;
formulae=handles.StressComponents_tab.Data;
p=ones(1,handles.OutputData.FitN); % this is for evaluation attempt of formulas
PinF=zeros(1,handles.OutputData.FitN);
for fn=1:6
    for pn=1:FitN
        if contains(formulae{fn},['p(' num2str(pn) ')'])
            PinF(pn)=PinF(pn)+1;
        end
    end
end

if sum(PinF<1)>0
    iserr(2)=1;
else
    try
        for fn=1:6
            eval([formulae{fn} ';']);
        end
    catch
        iserr(2)=1;
    end
end
if iserr(2)
    errordlg({'The provided formulae for stress components could not be evaluated from the following possible reasons:';'';'1) They contain a syntax mistake';'2) You did not refer to all parameters which you want to fit';'3) You refered to more parameters than is the number you have set'},'Wrong parameters to fit');
else
    handles.OutputData.StressComponents=formulae;
end

handles.OutputData.Geometry=handles.MeasurementGeometry_pop.String{handles.MeasurementGeometry_pop.Value};


% --- Executes on button press in IdealPolycrystal_rad.
function IdealPolycrystal_rad_Callback(hObject, eventdata, handles)
% hObject    handle to IdealPolycrystal_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IdealPolycrystal_rad
handles.LoadODFfile_btn.Enable='off';
handles.ODFfilename_ed.Enable='off';

guidata(hObject, handles);


% --- Executes on button press in PolycrystalWithTexture_rad.
function PolycrystalWithTexture_rad_Callback(hObject, eventdata, handles)
% hObject    handle to PolycrystalWithTexture_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PolycrystalWithTexture_rad
handles.LoadODFfile_btn.Enable='on';
handles.ODFfilename_ed.Enable='on';

guidata(hObject, handles);



function FitN_ed_Callback(hObject, eventdata, handles)
% hObject    handle to FitN_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FitN_ed as text
%        str2double(get(hObject,'String')) returns contents of FitN_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if ~(isnumeric(n) && n==round(n))
        iserr=1;
    else
        if n<1
            iserr=1;
        end
    end
catch
    iserr=1;
end
if iserr==1
    errordlg('Number of parameters to fit should be an integer higher than 0.','Wrong number format');
    hObject.String='1';
end


% --- Executes on selection change in MeasurementGeometry_pop.
function MeasurementGeometry_pop_Callback(hObject, eventdata, handles)
% hObject    handle to MeasurementGeometry_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MeasurementGeometry_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MeasurementGeometry_pop


% --- Executes during object creation, after setting all properties.
function MeasurementGeometry_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasurementGeometry_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in XEC_rad.
function XEC_rad_Callback(hObject, eventdata, handles)
% hObject    handle to XEC_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XEC_rad
handles.Orientation_tab.Enable='off';
handles.IndividualDatagroupOrientation_chbx.Enable='off';

handles.StressFactorsModel_pop.Enable='off';
handles.IdealPolycrystal_rad.Enable='off';
handles.PolycrystalWithTexture_rad.Enable='off';
handles.LoadODFfile_btn.Enable='off';
handles.ODFfilename_ed.Enable='off';

guidata(hObject, handles);


% --- Executes on selection change in XECmodel_pop.
function XECmodel_pop_Callback(hObject, eventdata, handles)
% hObject    handle to XECmodel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XECmodel_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XECmodel_pop
switch hObject.String{hObject.Value}
    case 'Kroener'
        handles.KroenerCoeff_txt.Enable='on';
        handles.KroenerCoeff_ed.Enable='on';
    otherwise
        handles.KroenerCoeff_txt.Enable='off';
        handles.KroenerCoeff_ed.Enable='off';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function XECmodel_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XECmodel_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function KroenerCoeff_ed_Callback(hObject, eventdata, handles)
% hObject    handle to KroenerCoeff_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of KroenerCoeff_ed as text
%        str2double(get(hObject,'String')) returns contents of KroenerCoeff_ed as a double
iserr=0;
try
    n=str2num(hObject.String);
    if isempty(n) | (n<0 | n>1)
        iserr=1;
    end
catch
    iserr=1;
end
if iserr==1
    errordlg('The weigth factor should be a number between 0 and 1 or equal.');
    hObject.String='0.5';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function KroenerCoeff_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KroenerCoeff_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
