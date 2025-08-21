function varargout = PSB_GUI(varargin)
% PSB_GUI MATLAB code for PSB_GUI.fig
%      PSB_GUI, by itself, creates a new PSB_GUI or raises the existing
%      singleton*.
%
%      H = PSB_GUI returns the handle to a new PSB_GUI or the handle to
%      the existing singleton*.
%
%      PSB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSB_GUI.M with the given input arguments.
%
%      PSB_GUI('Property','Value',...) creates a new PSB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PSB_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PSB_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PSB_GUI

% Last Modified by GUIDE v2.5 21-Aug-2025 15:25:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PSB_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PSB_GUI_OutputFcn, ...
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


% --- Executes just before PSB_GUI is made visible.
function PSB_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PSB_GUI (see VARARGIN)

% Choose default command line output for PSB_GUI
handles.output = hObject;

% initial values for the some components
% Material constants for Fe
handles.MaterialConstants_tab.ColumnName={char(957),'E (GPa)','G (GPa)',[char(957) char(8321) char(8320) char(8320)],['E' char(8321) char(8320) char(8320) ' (GPa)'],['G' char(8321) char(8320) char(8320) ' (GPa)']};
handles.DataTable_tab.ColumnName={'h','k','l',['2' char(952) ' (deg)'],[char(968) ' (deg)'],[char(966) ' (deg)'],'FWHM (deg)',[char(946) ' (deg)'],'label','In group'};
handles.YoungModulus_rad.String=['E, ' char(957) ', G'];
% handles.SijMat_rad.String=['S' char(7522) char(11388)];
% handles.CijMat_rad.String=['C' char(7522) char(11388)];
handles.MaterialConstants_tab.Data={0.37 85 30 NaN NaN NaN};
handles.Stiffness_tab.Data={NaN NaN NaN 0 0 0; NaN NaN NaN 0 0 0; NaN NaN NaN 0 0 0; 0 0 0 NaN 0 0; 0 0 0 0 NaN 0; 0 0 0 0 0 NaN};
% handles.Wavelength_txt2.String=char(8491);

opts = detectImportOptions('Constants_database.xlsx');
opts.VariableNames={'Label','Material','Reference','System','nu','E','G'};
handles.MaterialConstants_database=readtable('Constants_database.xlsx',opts);
opts = detectImportOptions('Stiffness_database.xlsx');
opts.VariableNames={'Label','Material','Reference','System','s11','s12','s13','s33','s44','c11','c12','c13','c33','c44'};
handles.StiffnessMatrices_database=readtable('Stiffness_database.xlsx',opts);
handles.MaterialDatabase_list.String=handles.MaterialConstants_database.Label;
handles.MaterialDatabase_list.Value=1;
switch handles.MaterialConstants_database.System{1}
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'hexagonal'
        handles.System_pop.Value=3;
end

% Default options for information about the material
handles.GroupsModels(1).sys=handles.MaterialConstants_database.System{1};
handles.GroupsModels(1).Material.av=0;

% Default options for fitting models
handles.GroupsModels(1).PeakShiftModel='IsotropicStress';
handles.GroupsModels(1).SizeOptions.GrainsShape='Sphere';
handles.GroupsModels(1).SizeOptions.Fitted=1;
handles.GroupsModels(1).SizeOptions.ParGuess=NaN;
handles.GroupsModels(1).SizeOptions.ParFix=false;
handles.GroupsModels(1).SizeOptions.ParLB=NaN;
handles.GroupsModels(1).SizeOptions.ParUB=NaN;
handles.GroupsModels(1).SizeOptions.GrownDirectionOpt='Sample normal';
handles.GroupsModels(1).SizeOptions.GrownDirectionHKL=[0 0 1];
% handles.GroupsModels(1).SizeOptions.GrownDirectionOpt2='Direction in sample surface specified by azimuth';
handles.GroupsModels(1).SizeOptions.GrownDirectionHKL2=[1 0 0];
handles.GroupsModels(1).SizeOptions.GrownDirectionPerpAzimuth=0;
% handles.GroupsModels(1).MicrostrainOptions.PeakBroadeningModel='Isotropic';
handles.GroupsModels(1).MicrostrainOptions.MicrostrainModel='Isotropic';
handles.GroupsModels(1).MicrostrainOptions.ParGuess=NaN;
handles.GroupsModels(1).MicrostrainOptions.ParFix=false;
handles.GroupsModels(1).MicrostrainOptions.ParLB=NaN;
handles.GroupsModels(1).MicrostrainOptions.ParUB=NaN;
handles.GroupsModels(1).MicrostrainOptions.Fitted=1;
handles.GroupsModels(1).FitSF=0;
handles.GroupsModels(1).AdvancedShiftOptions.MaterialInfo='XEC';
handles.GroupsModels(1).AdvancedShiftOptions.XECmodel='Voigt';
handles.GroupsModels(1).AdvancedShiftOptions.KroenerCoeff=0.5;
handles.GroupsModels(1).AdvancedShiftOptions.OrientationAngles={0 0 0 'z' 'x' 'z' 'extrinsic'};
handles.GroupsModels(1).AdvancedShiftOptions.DatagroupsIndividual=0;
handles.GroupsModels(1).AdvancedShiftOptions.StressFactorsMethod='Reuss';
handles.GroupsModels(1).AdvancedShiftOptions.ODFtype='IdealPolycrystal';
handles.GroupsModels(1).AdvancedShiftOptions.ODFfilename='';
handles.GroupsModels(1).AdvancedShiftOptions.FitN=1;
handles.GroupsModels(1).AdvancedShiftOptions.StressComponents={'p(1)'; 'p(1)'; '0'; '0'; '0'; '0'};
handles.GroupsModels(1).AdvancedShiftOptions.Geometry='GAXRD';

handles.GroupSelected=1;

handles.Wavelength=NaN;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PSB_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PSB_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in System_pop.
function System_pop_Callback(hObject, eventdata, handles)
% hObject    handle to System_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns System_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from System_pop
switch hObject.String{hObject.Value}
    case {'cubic_fcc';'cubic_bcc'}
        handles.FitHcpLP_chbx.Enable='off';
    case 'hexagonal'
        handles.FitHcpLP_chbx.Enable='on';
end
I=handles.GroupSelected;
handles.GroupsModels(I).sys=hObject.String{hObject.Value};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function System_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to System_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MaterialDatabase_list.
function MaterialDatabase_list_Callback(hObject, eventdata, handles)
% hObject    handle to MaterialDatabase_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MaterialDatabase_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MaterialDatabase_list
ind=hObject.Value;
I=handles.GroupSelected;
handles.GroupsModels(I).Material.DatabaseRow=ind;
switch handles.MaterialConstants_database.System{ind}
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'hexagonal'
        handles.System_pop.Value=3;
end

if handles.YoungModulus_rad.Value==1
    nu=handles.MaterialConstants_database.nu(ind);
    E=handles.MaterialConstants_database.E(ind);
    G=handles.MaterialConstants_database.G(ind);
    handles.MaterialConstants_tab.Data={nu E G NaN NaN NaN};
    handles.Reference_txt.Enable='on';
    handles.Reference_txt.String=['Reference: ' handles.MaterialConstants_database.Reference{ind}];
    switch handles.MaterialConstants_database.System{ind}
        case 'cubic_fcc'
            handles.System_pop.Value=1;
        case 'cubic_bcc'
            handles.System_pop.Value=2;
        case 'hexagonal'
            handles.System_pop.Value=3;
    end

    % Update material info for selected dataset
    tab=handles.MaterialConstants_tab.Data;
    handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
    handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
    handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
    handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
    handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
    handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)
elseif handles.SijMat_rad.Value==1
    switch handles.StiffnessMatrices_database.System{ind}
        case 'cubic_fcc'
            handles.System_pop.Value=1;
        case 'cubic_bcc'
            handles.System_pop.Value=2;
        case 'hexagonal'
            handles.System_pop.Value=3;
    end
    s11=handles.StiffnessMatrices_database.s11(ind);
    s12=handles.StiffnessMatrices_database.s12(ind);
    s13=handles.StiffnessMatrices_database.s13(ind);
    s33=handles.StiffnessMatrices_database.s33(ind);
    s44=handles.StiffnessMatrices_database.s44(ind);
    switch handles.StiffnessMatrices_database.System{ind}
        case {'cubic_fcc','cubic_bcc'}
            S=[s11 s12 s12 0 0 0; s12 s11 s12 0 0 0; s12 s12 s11 0 0 0; 0 0 0 s44 0 0; 0 0 0 0 s44 0; 0 0 0 0 0 s44];
            handles.Stiffness_tab.Data={s11 s12 s12 0 0 0; s12 s11 s12 0 0 0; s12 s12 s11 0 0 0; 0 0 0 s44 0 0; 0 0 0 0 s44 0; 0 0 0 0 0 s44};
        case 'hexagonal'
            S=[s11 s12 s13 0 0 0; s12 s11 s13 0 0 0; s13 s13 s33 0 0 0; 0 0 0 s44 0 0; 0 0 0 0 s44 0; 0 0 0 0 0 2*(s11-s12)];
            handles.Stiffness_tab.Data={s11 s12 s13 0 0 0; s12 s11 s13 0 0 0; s13 s13 s33 0 0 0; 0 0 0 s44 0 0; 0 0 0 0 s44 0; 0 0 0 0 0 2*(s11-s12)};
    end
    C=inv(S);
    nu100=-S(1,2)/S(1,1);
    E100=1/S(1,1);
    G100=1/S(4,4);
    [nu,E,G,~]=ElasticConstants_VoigtReussHill(C,handles.StiffnessMatrices_database.System{ind});
    handles.MaterialConstants_tab.Data={nu E G nu100 E100 G100};
    handles.Reference_txt.Enable='on';
    handles.Reference_txt.String=['Reference: ' handles.StiffnessMatrices_database.Reference{ind}];

    % Update material info for selected dataset
    handles.GroupsModels(I).Material.Sij=cell2mat(handles.Stiffness_tab.Data);
    handles.GroupsModels(I).Material.Cij=inv(handles.GroupsModels(I).Material.Sij);
else % Cij
    switch handles.StiffnessMatrices_database.System{ind}
        case 'cubic_fcc'
            handles.System_pop.Value=1;
        case 'cubic_bcc'
            handles.System_pop.Value=2;
        case 'hexagonal'
            handles.System_pop.Value=3;
    end
    c11=handles.StiffnessMatrices_database.c11(ind);
    c12=handles.StiffnessMatrices_database.c12(ind);
    c13=handles.StiffnessMatrices_database.c13(ind);
    c33=handles.StiffnessMatrices_database.c33(ind);
    c44=handles.StiffnessMatrices_database.c44(ind);
    switch handles.StiffnessMatrices_database.System{ind}
        case {'cubic_fcc','cubic_bcc'}
            C=[c11 c12 c12 0 0 0; c12 c11 c12 0 0 0; c12 c12 c11 0 0 0; 0 0 0 c44 0 0; 0 0 0 0 c44 0; 0 0 0 0 0 c44];
            handles.Stiffness_tab.Data={c11 c12 c12 0 0 0; c12 c11 c12 0 0 0; c12 c12 c11 0 0 0; 0 0 0 c44 0 0; 0 0 0 0 c44 0; 0 0 0 0 0 c44};
        case 'hexagonal'
            C=[c11 c12 c13 0 0 0; c12 c11 c13 0 0 0; c13 c13 c33 0 0 0; 0 0 0 c44 0 0; 0 0 0 0 c44 0; 0 0 0 0 0 0.5*(c11-c12)];
            handles.Stiffness_tab.Data={c11 c12 c13 0 0 0; c12 c11 c13 0 0 0; c13 c13 c33 0 0 0; 0 0 0 c44 0 0; 0 0 0 0 c44 0; 0 0 0 0 0 0.5*(c11-c12)};
    end
    S=inv(C);
    nu100=-S(1,2)/S(1,1);
    E100=1/S(1,1);
    G100=1/S(4,4);
    [nu,E,G,~]=ElasticConstants_VoigtReussHill(C,handles.StiffnessMatrices_database.System{ind});
    handles.MaterialConstants_tab.Data={nu E G nu100 E100 G100};
    handles.Reference_txt.Enable='on';
    handles.Reference_txt.String=['Reference: ' handles.StiffnessMatrices_database.Reference{ind}];

    % Update material info for selected dataset
    handles.GroupsModels(I).Material.Cij=cell2mat(handles.Stiffness_tab.Data);
    handles.GroupsModels(I).Material.Sij=inv(handles.GroupsModels(I).Material.Cij);
end
System_pop_Callback(handles.System_pop,eventdata,handles)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MaterialDatabase_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaterialDatabase_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FitAll_btn.
function FitAll_btn_Callback(hObject, eventdata, handles)
% hObject    handle to FitAll_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GlobOpt=GetGlobalOptions(handles);
handles.GlobOpt=GlobOpt; % it is used for plotting results afterwards
% I have to erase field with results from handles variable, because they
% could have variable sizes and it could be problematic
if isfield(handles,'PeakShiftResults')
    handles=rmfield(handles,'PeakShiftResults');
end
if isfield(handles,'PeakBroadeningResults_fwhm')
    handles=rmfield(handles,'PeakBroadeningResults_fwhm');
end
if isfield(handles,'PeakBroadeningResults_beta')
    handles=rmfield(handles,'PeakBroadeningResults_beta');
end
% close all warnings and error dialogs
hf=findall(groot,'type','figure');
for n=1:numel(hf)
    if contains(hf(n).Tag,'Msgbox_')
        close(hf(n));
    end
end

% Check wavelength input
if isnan(handles.Wavelength)
    errordlg('The known wavelength is neccessary for the fit!','Wrong input')
    iserrWL=1;
else
    iserrWL=0;
end
 
if iserrWL==0
    % Peak shift evaluation
    % Check availability of data
    iserrPS=zeros(1,size(handles.DataGroup,2));
    if sum(handles.DataAvailability(1:4)+handles.DataAvailability(6))<5 % h k l psi theta
        errordlg('There is not enough available data for the evaluation of peak shift.','Not enough data');
        for g=1:size(handles.DataGroup,2)
            iserrPS(g)=3;
            handles.PeakShiftResults(g).av=0;
        end
    end

for g=1:size(handles.DataGroup,2)
    if ~iserrPS(g)
        if (GlobOpt.FitGroups==0 & g>1) | (GlobOpt.FitGroups==1 & g==1)
            handles.PeakShiftResults(g).av=0;
        else
            FitOpt=GetFitOptions(handles,g); % TAHLE FUNKCE NASBÍRÁ MODELY, MOŽNÁ TO PŘEPSAT, VNOŘIT DO FOR CYKLU
            if strcmp(FitOpt.PeakShiftModel,'AdvancedShiftOpt') & FitOpt.Material.av==0
                handles.PeakShiftResults(g).av=0;
                iserrPS(g)=4;
                errordlg(['Peak shift (G#' num2str(g-1) '): It is neccessary to have known material constants if you want to fit with ''Advanced options''.'],'Material constants not available');
            end
            if ~iserrPS(g)
                if FitOpt.Material.av==1
                    switch FitOpt.Material.Option
                        case 'E_nu_G'
                            warndlg(['Peak shift (G#' num2str(g-1) '): Not a complete stiffness tensor is known, it will be recomputed as for isotropic material.','Material constants not fully available']);
                            s11=1/FitOpt.Material.E;
                            s12=FitOpt.Material.nu/FitOpt.Material.E;
                            s44=1/FitOpt.Material.G;
                            Sij=[s11 s12 s12  0   0    0;
                                 s12 s11 s12  0   0    0;
                                 s12 s12 s11  0   0    0;
                                 0    0   0  s44  0    0;
                                 0    0   0   0  s44   0;
                                 0    0   0   0   0   s44];
                            Cij=inv(Sij);
                        case 'Sij'
                            Sij=FitOpt.Material.Sij;
                            Cij=inv(Sij);
                        case 'Cij'
                            Cij=FitOpt.Material.Cij;
                            Sij=inv(Cij);
                    end
                end
                ind=find(handles.DataGroup(:,g)==1);
                myData.h=handles.Data.h(ind);
                myData.k=handles.Data.k(ind);
                myData.l=handles.Data.l(ind);
                myData.tth=handles.Data.tth(ind);
                myData.psi=handles.Data.psi(ind);
                myData.phi=handles.Data.phi(ind);
%                 disp([myData.h myData.k myData.l myData.tth myData.psi myData.phi])
                if strcmp(FitOpt.PeakShiftModel,'AdvancedShiftOpt')
                    disp(['Dataset #' num2str(g-1) ': Fitting peak shift with advanced options...']);
                    if FitOpt.AdvancedShiftOptions.DatagroupsIndividual==0
                        [results,ahkl,ahkl_T,iserrPS(g)]=EvaluateShiftAdvanced(myData,GlobOpt,FitOpt,Sij,Cij,1);
                    else
                        [results,ahkl,ahkl_T,iserrPS(g)]=EvaluateShiftAdvanced(myData,GlobOpt,FitOpt,Sij,Cij,g-1);
                    end
                else
                    disp(['Dataset #' num2str(g-1) ': Fitting peak shift...']);
                    disp(['Model: ' FitOpt.PeakShiftModel]);
                    [results,ahkl,ahkl_T,iserrPS(g)]=EvaluateShift(myData,GlobOpt,FitOpt);
                end
                
                switch iserrPS(g)
                    case -1
                        handles.PeakShiftResults(g).results=results;
                        handles.PeakShiftResults(g).ahkl=ahkl;
                        handles.PeakShiftResults(g).ahkl_T=ahkl_T;
                        handles.PeakShiftResults(g).av=1;
                        warndlg(['Peak shift (G#' num2str(g-1) '): The parameters resulting from the fit are complex or they have large errors. They probably correlate.'],'Large errors or complex values of resulting parameters');
                        guidata(hObject,handles)
                    case 0
                        handles.PeakShiftResults(g).results=results;
                        handles.PeakShiftResults(g).ahkl=ahkl;
                        handles.PeakShiftResults(g).ahkl_T=ahkl_T;
                        handles.PeakShiftResults(g).av=1;
                        guidata(hObject,handles);
                    case 1
                        errordlg(['Peak shift (G#' num2str(g-1) '): The model of the fit has more parameters than available data points.'],'Not enough data points');
                        handles.PeakShiftResults(g).av=0;
                    case 2
                        errordlg(['Peak shift (G#' num2str(g-1) '): Initial guess of the fitting parameters has complex values. Please, check your data.'],'Complex values of parameters at initial guess');
                        handles.PeakShiftResults(g).av=0;
                end
            else
                handles.PeakShiftResults(g).av=0;
            end
        end
    end
end

% Peak broadening from FWHM
% Check availability of data
iserrPB_fwhm=zeros(1,size(handles.DataGroup,2));
if sum(handles.DataAvailability(1:3))+sum(handles.DataAvailability(6:7))<5 % h k l theta fwhm
    errordlg('There is not enough available data for the evaluation of peak broadening from FWHM parameter.','Not enough data');
    for g=1:size(handles.DataGroup,2)
        iserrPB_fwhm(g)=3;
        handles.PeakBroadeningResults_fwhm(g).av=0;
    end
end

iserrPB_beta=zeros(1,size(handles.DataGroup,2));
if sum(handles.DataAvailability(1:3))+sum(handles.DataAvailability([6 8]))<5 % h k l theta beta
    errordlg('There is not enough available data for the evaluation of peak broadening from beta parameter.','Not enough data');
    for g=1:size(handles.DataGroup,2)
        iserrPB_beta(g)=3;
        handles.PeakBroadeningResults_beta(g).av=0;
    end
end

for g=1:size(handles.DataGroup,2)
    if (GlobOpt.FitGroups==0 & g>1) | (GlobOpt.FitGroups==1 & g==1)
        handles.PeakBroadeningResults_fwhm(g).av=0;
        handles.PeakBroadeningResults_beta(g).av=0;
    else
        FitOpt=GetFitOptions(handles,g); % TAHLE FUNKCE NASBÍRÁ MODELY, MOŽNÁ TO PŘEPSAT, VNOŘIT DO FOR CYKLU
        if FitOpt.FitSize==0 & FitOpt.FitMicrostrain==0 & FitOpt.FitSF==0
            errordlg('At least one effect contributing to the peak broadening should be selected.','Wrong fitting option');
            handles.PeakBroadeningResults_fwhm(g).av=0;
            handles.PeakBroadeningResults_beta(g).av=0;
        else
            ind=find(handles.DataGroup(:,g)==1);
            if handles.SubstractInstrumental_chbx.Value==1
                [fwhmIP,betaIP]=ComputeInstrumentalBroadening(handles.IP_Parameters,handles.Data.tth,handles.Data.psi);
                N=FitOpt.BroadeningNFitOrder;
                broaddata_fwhm=(handles.Data.fwhm(ind).^N-fwhmIP(ind).^N).^(1/N);
                broaddata_beta=(handles.Data.beta(ind).^N-betaIP(ind).^N).^(1/N);
            else
                broaddata_fwhm=handles.Data.fwhm(ind);
                broaddata_beta=handles.Data.beta(ind);
            end

            myData.h=handles.Data.h(ind);
            myData.k=handles.Data.k(ind);
            myData.l=handles.Data.l(ind);
            myData.tth=handles.Data.tth(ind);
            myData.psi=handles.Data.psi(ind);
            myData.phi=handles.Data.phi(ind);
%             disp([myData.h myData.k myData.l myData.tth myData.psi myData.phi])
            % Fit FWHM
            if ~iserrPB_fwhm(g)
                [results,broad_0,broad_T,iserrPB_fwhm(g)]=EvaluateBroadening(broaddata_fwhm,myData,GlobOpt,FitOpt);
                Nth=EstimateStrainSizeEffect(myData.tth,broaddata_fwhm,GlobOpt.lam);
                if abs(Nth-FitOpt.BroadeningNFitOrder)>0.1
                    warndlg(['Peak broadening (FWHM:G#' num2str(g-1) '): The estimated ratio between the strain/size influence is different from the value corresponding to ''Broadening fitting order'' parameter. Consider changing it. Theoreticaly estimated value: ' num2str(Nth)],'Broadening fitting order');
                end

                switch iserrPB_fwhm(g)
                    case -1
                        handles.PeakBroadeningResults_fwhm(g).results=results;
                        handles.PeakBroadeningResults_fwhm(g).broad=broad_0;
                        handles.PeakBroadeningResults_fwhm(g).broad_T=broad_T;
                        handles.PeakBroadeningResults_fwhm(g).av=1;
                        warndlg(['Peak broadening (FWHM:G#' num2str(g-1) '): The parameters resulting from the fit are complex or they have large errors. They probably correlate.'],'Large errors or complex values of resulting parameters');
                        guidata(hObject,handles)
                    case 0
                        handles.PeakBroadeningResults_fwhm(g).results=results;
                        handles.PeakBroadeningResults_fwhm(g).broad=broad_0;
                        handles.PeakBroadeningResults_fwhm(g).broad_T=broad_T;
                        handles.PeakBroadeningResults_fwhm(g).av=1;
                        guidata(hObject,handles);
                    case 1
                        errordlg(['Peak broadening (FWHM:G#' num2str(g-1) '): The model of the fit has more parameters than available data points.'],'Not enough data points');
                        handles.PeakBroadeningResults_fwhm(g).av=0;
                    case 2
                        errordlg(['Peak broadening (FWHM:G#' num2str(g-1) '): Initial guess of the fitting parameters has complex values. Please, check your data.'],'Complex values of parameters at initial guess');
                        handles.PeakBroadeningResults_fwhm(g).av=0;
                end
            else
                handles.PeakBroadeningResults_fwhm(g).av=0;
            end
            
            % Fit beta
            if ~iserrPB_beta(g)
                [results,broad_0,broad_T,iserrPB_beta(g)]=EvaluateBroadening(broaddata_beta,myData,GlobOpt,FitOpt);
                Nth=EstimateStrainSizeEffect(myData.tth,broaddata_beta,GlobOpt.lam);
                if abs(Nth-FitOpt.BroadeningNFitOrder)>0.1
                    warndlg(['Peak broadening (beta:G#' num2str(g-1) '): The estimated ratio between the strain/size influence is different from the value corresponding to ''Broadening fitting order'' parameter. Consider changing it. Theoreticaly estimated value: ' num2str(Nth)],'Broadening fitting order');
                end

                switch iserrPB_beta(g)
                    case -1
                        handles.PeakBroadeningResults_beta(g).results=results;
                        handles.PeakBroadeningResults_beta(g).broad=broad_0;
                        handles.PeakBroadeningResults_beta(g).broad_T=broad_T;
                        handles.PeakBroadeningResults_beta(g).av=1;
                        warndlg(['Peak broadening (beta:G#' num2str(g-1) '): The parameters resulting from the fit are complex or they have large errors. They probably correlate.'],'Large errors or complex values of resulting parameters');
                        guidata(hObject,handles)
                    case 0
                        handles.PeakBroadeningResults_beta(g).results=results;
                        handles.PeakBroadeningResults_beta(g).broad=broad_0;
                        handles.PeakBroadeningResults_beta(g).broad_T=broad_T;
                        handles.PeakBroadeningResults_beta(g).av=1;
                        guidata(hObject,handles);
                    case 1
                        errordlg(['Peak broadening (beta:G#' num2str(g-1) '): The model of the fit has more parameters than available data points.'],'Not enough data points');
                        handles.PeakBroadeningResults_beta(g).av=0;
                    case 2
                        errordlg(['Peak broadening (beta:G#' num2str(g-1) '): Initial guess of the fitting parameters has complex values. Please, check your data.'],'Complex values of parameters at initial guess');
                        handles.PeakBroadeningResults_beta(g).av=0;
                end
            else
                handles.PeakBroadeningResults_beta(g).av=0;
            end
        end
    end
end

PlotData(handles);
str=DispResults(handles,GlobOpt);
handles.Status_list.String=str;
guidata(hObject,handles)
end
%%% END OF FIT SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in FitGroups_chbx.
function FitGroups_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to FitGroups_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SF_chbx.
function SF_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to SF_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SF_chbx
I=handles.GroupSelected;
handles.GroupsModels(I).FitSF=handles.SF_chbx.Value;
guidata(hObject, handles);


% --- Executes on button press in Size_chbx.
function Size_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to Size_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Size_chbx
I=handles.GroupSelected;
if hObject.Value==0
    handles.GroupsModels(I).SizeOptions.Fitted=0;
else
    handles.GroupsModels(I).SizeOptions.Fitted=1;
end
guidata(hObject, handles);


% --- Executes on button press in ImportData_btn.
function ImportData_btn_Callback(hObject, eventdata, handles)
% hObject    handle to ImportData_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastFilePath')
    [FileName,PathName,FilterIndex] = uigetfile({'*.xls,*.xlsx','Excel tables (*.xls, *.xlsx)'; '*.*',  'All Files (*.*)'},'Select File to Open',handles.LastFilePath);
    handles.LastFilePath=PathName;
else
    [FileName,PathName,FilterIndex] = uigetfile({'*.xls,*.xlsx','Excel tables (*.xls, *.xlsx)'; '*.*',  'All Files (*.*)'},'Select File to Open');
    handles.LastFilePath=PathName;
end
if FileName~=0
    if FilterIndex==1
        % excel spreadsheet
        mytable=readtable([PathName '\' FileName],'FileType','spreadsheet');
    else
        % textfile, default option is not always able to read the file,
        % header and "table cosmetics" (such as row of _ or -) are
        % problematic. Some extra code is neccessary
        nhMax=100;
        for nh=0:nhMax
            mytable=readtable([PathName '\' FileName],'FileType','delimitedtext','numheaderlines',nh,'importerrorrule','omitrow');
            ColumnNames=mytable.Properties.VariableNames;
            if size(mytable,1)>1 && size(mytable,2)>1 && ~strcmp(ColumnNames{1},'Var1') && ~strcmp(ColumnNames{2},'Var2')
                break;
            end
        end
        % Now, let's try to scan the header, if there is information about
        % the wavelength
        lam=ScanHeaderForWavelength([PathName '\' FileName],nh);
        handles.Wavelength=lam;
        if isnan(lam)
            warndlg('Information about the wavelength has not been automaticaly found. Check the data import with "Manual selection" tool.');
        end
    end

    b=zeros(1,9);
    col=zeros(1,9);
    formula=cell(1,9);
    DataUnitsOpt=cell(1,9);
    ind=find(ismember(mytable.Properties.VariableDescriptions,{'h','h_','h,','h;';'H','H_','H,','H;'})==1);
    if ~isempty(ind)
        b(1)=1;
        col(1)=ind;
        handles.Data.h=mytable.(ColumnNames{ind});
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'k','k_','k,','k;';'K','K_','K,','K;'})==1);
    if ~isempty(ind)
        b(2)=1;
        col(2)=ind;
        handles.Data.k=mytable.(ColumnNames{ind});
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'l','l_','l,','l;';'L','L_','L,','L;'})==1);
    if ~isempty(ind)
        b(3)=1;
        col(3)=ind;
        handles.Data.l=mytable.(ColumnNames{ind});
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Chi','Chi_','Chi,','Chi;',...
        'chi','chi_','chi,','chi;',...
        'CHI','CHI_','CHI,','CHI;',...
        'Psi','Psi_','Psi,','Psi;',...
        'psi','psi_','psi,','psi;',...
        'PSI','PSI_','PSI,','PSI;',...
        'Chi (deg)','Chi (deg)_','Chi (deg),','Chi (deg);',...
        'chi (deg)','chi (deg)_','chi (deg),','chi (deg);',...
        'CHI (deg)','CHI (deg)_','CHI (deg),','CHI (deg);',...
        'CHI (DEG)','CHI (DEG)_','CHI (DEG),','CHI (DEG);',...
        'Psi (deg)','Psi (deg)_','Psi (deg),','Psi (deg);',...
        'psi (deg)','psi (deg)_','psi (deg),','psi (deg);',...
        'PSI (deg)','PSI (deg)_','PSI (deg),','PSI (deg);',...
        'PSI (DEG)','PSI (DEG)_','PSI (DEG),','PSI (DEG);'})==1);
    if ~isempty(ind)
        b(4)=1;
        col(4)=ind;
        
        hlp=lower(mytable.Properties.VariableDescriptions{ind});
        if contains(hlp,'deg')
            DataUnitsOpt{4}='deg';
            handles.Data.psi=mytable.(ColumnNames{ind});
        elseif contains(hlp,'rad')
            DataUnitsOpt{4}='rad';
            handles.Data.psi=mytable.(ColumnNames{ind})*180/pi;
        else
            DataUnitsOpt{4}='uknown';
            handles.Data.psi=mytable.(ColumnNames{ind})*NaN;
            b(4)=0;
        end
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Phi','Phi_','Phi,','Phi;',...
        'phi','phi_','phi,','phi;',...
        'PHI','PHI_','PHI,','PHI;',...
        'Phi (deg)','Phi (deg)_','Phi (deg),','Phi (deg);',...
        'phi (deg)','phi (deg)_','phi (deg),','phi (deg);',...
        'PHI (deg)','PHI (deg)_','PHI (deg),','PHI (deg);',...
        'PHI (DEG)','PHI (DEG)_','PHI (DEG),','PHI (DEG);'})==1);
    if ~isempty(ind)
        b(5)=1;
        col(5)=ind;
        
        hlp=lower(mytable.Properties.VariableDescriptions{ind});
        if contains(hlp,'deg')
            DataUnitsOpt{5}='deg';
            handles.Data.phi=mytable.(ColumnNames{ind});
        elseif contains(hlp,'rad')
            DataUnitsOpt{5}='rad';
            handles.Data.phi=mytable.(ColumnNames{ind})*180/pi;
        else
            DataUnitsOpt{5}='uknown';
            handles.Data.phi=mytable.(ColumnNames{ind})*NaN;
            b(5)=0;
        end
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Tth','Tth_','Tth,','Tth;',...
        'tth','tth_','tth,','tth;',...
        'TTH','TTH_','TTH,','TTH;',...
        '2theta','2theta_','2theta,','2theta;',...
        '2Theta','2Theta_','2Theta,','2Theta;',...
        '2THETA','2THETA_','2THETA,','2THETA;',...
        'Theta2','Theta2_','Theta2,','Theta2;',...
        'THETA2','THETA2_','THETA2,','THETA2;',...
        'Tth (deg)','Tth (deg)_','Tth (deg),','Tth (deg);',...
        'tth (deg)','tth (deg)_','tth (deg),','tth (deg);',...
        'TTH (deg)','TTH (deg)_','TTH (deg),','TTH (deg);',...
        'TTH (DEG)','TTH (DEG)_','TTH (DEG),','TTH (DEG);',...
        '2theta (deg)','2theta (deg)_','2theta (deg),','2theta (deg);',...
        '2Theta (deg)','2Theta (deg)_','2Theta (deg),','2Theta (deg);',...
        '2THETA (deg)','2THETA (deg)_','2THETA (deg),','2THETA (deg);',...
        '2THETA (DEG)','2THETA (DEG)_','2THETA (DEG),','2THETA (DEG);',...
        'Theta2 (deg)','Theta2_ (deg)','Theta2 (deg),','Theta2 (deg);',...
        'THETA2 (deg)','THETA2 (deg)_','THETA2 (deg),','THETA2 (deg);',...
        'THETA2 (DEG)','THETA2 (DEG)_','THETA2 (DEG),','THETA2 (DEG);'})==1);
    if ~isempty(ind)
        b(6)=1;
        col(6)=ind;

        hlp=lower(mytable.Properties.VariableDescriptions{ind});
        if contains(hlp,'deg')
            DataUnitsOpt{6}='deg';
            handles.Data.tth=mytable.(ColumnNames{ind});
        elseif contains(hlp,'rad')
            DataUnitsOpt{6}='rad';
            handles.Data.tth=mytable.(ColumnNames{ind})*180/pi;
        elseif contains(hlp,'A-1')
            DataUnitsOpt{6}='A-1';
            hlp0=mytable.(ColumnNames{ind});
            handles.Data.tth=2*asind(hlp0*lam/(4*pi));
        else
            DataUnitsOpt{6}='uknown';
            handles.Data.tth=mytable.(ColumnNames{ind})*NaN;
            b(6)=0;
        end
    else
        % try to find theta column instead
        ind=find(ismember(mytable.Properties.VariableDescriptions,{'Th','Th_','Th,','Th;',...
        'th','th_','th,','th;',...
        'TH','TH_','TH,','TH;',...
        'Theta','Theta_','Theta,','Theta;',...
        'theta','theta_','theta,','theta;',...
        'Th (deg)','Th_ (deg)','Th (deg),','Th (deg);',...
        'th (deg)','th (deg)_','th (deg),','th (deg);',...
        'TH (deg)','TH (deg)_','TH (deg),','TH (deg);',...
        'TH (DEG)','TH (DEG)_','TH (DEG),','TH (DEG);',...
        'Theta (deg)','Theta (deg)_','Theta (deg),','Theta (deg);',...
        'theta (deg)','theta (deg)_','theta (deg),','theta (deg);',...
        'THETA (deg)','THETA (deg)_','THETA (deg),','THETA (deg);',...
        'THETA (DEG)','THETA (DEG)_','THETA (DEG),','THETA (DEG);'})==1);
        if ~isempty(ind)
            b(6)=1;
            col(6)=0;
            formula{6}=['2*' ColumnNames{ind}];

            hlp=lower(mytable.Properties.VariableDescriptions{ind});
            if contains(hlp,'deg')
                DataUnitsOpt{6}='deg';
                handles.Data.tth=mytable.(ColumnNames{ind});
            elseif contains(hlp,'rad')
                DataUnitsOpt{6}='rad';
                handles.Data.tth=mytable.(ColumnNames{ind})*180/pi;
            elseif contains(hlp,'A-1')
                DataUnitsOpt{6}='A-1';
                hlp0=mytable.(ColumnNames{ind});
                handles.Data.tth=2*asind(hlp0*lam/(4*pi));
            else
                DataUnitsOpt{6}='uknown';
                handles.Data.tth=mytable.(ColumnNames{ind})*NaN;
                b(6)=0;
            end
        end
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Fwhm','Fwhm_','Fwhm,','Fwhm;',...
        'fwhm','fwhm_','fwhm,','fwhm;',...
        'FWHM','FWHM_','FWHM,','FWHM;',...
        'Fwhm (deg)','Fwhm (deg)_','Fwhm (deg),','Fwhm (deg);',...
        'fwhm (deg)','fwhm (deg)_','fwhm (deg),','fwhm (deg);',...
        'FWHM (deg)','FWHM (deg)_','FWHM (deg),','FWHM (deg);',...
        'FWHM (DEG)','FWHM (DEG)_','FWHM (DEG),','FWHM (DEG);'})==1);
    if ~isempty(ind)
        b(7)=1;
        col(7)=ind;
        
        hlp=lower(mytable.Properties.VariableDescriptions{ind});
        if contains(hlp,'deg')
            DataUnitsOpt{7}='deg';
            handles.Data.fwhm=mytable.(ColumnNames{ind});
        elseif contains(hlp,'rad')
            DataUnitsOpt{7}='rad';
            handles.Data.fwhm=mytable.(ColumnNames{ind})*180/pi;
        elseif contains(hlp,'A-1')
            DataUnitsOpt{7}='A-1';
            hlp0=mytable.(ColumnNames{ind});
            handles.Data.fwhm=2*hlp0*lam/(4*pi)./cosd(handles.Data.tth/2);
        else
            DataUnitsOpt{7}='uknown';
            handles.Data.fwhm=mytable.(ColumnNames{ind})*NaN;
            b(7)=0;
        end
    else
        % try to find hwhm column instead
        ind=find(ismember(mytable.Properties.VariableDescriptions,{'Hwhm','Hwhm_','Hwhm,','Hwhm;',...
        'hwhm','hwhm_','hwhm,','hwhm;',...
        'HWHM','HWHM_','HWHM,','HWHM;',...
        'Hwhm (deg)','Hwhm (deg)_','Hwhm (deg),','Hwhm (deg);',...
        'hwhm (deg)','hwhm (deg)_','hwhm (deg),','hwhm (deg);',...
        'HWHM (deg)','HWHM (deg)_','HWHM (deg),','HWHM (deg);',...
        'HWHM (DEG)','HWHM (DEG)_','HWHM (DEG),','HWHM (DEG);'})==1);
        if ~isempty(ind)
            b(7)=1;
            col(7)=0;
            formula{7}=['2*' ColumnNames{ind}];
            
            hlp=lower(mytable.Properties.VariableDescriptions{ind});
            if contains(hlp,'deg')
                DataUnitsOpt{7}='deg';
                handles.Data.fwhm=2*mytable.(ColumnNames{ind});
            elseif contains(hlp,'rad')
                DataUnitsOpt{7}='rad';
                handles.Data.fwhm=mytable.(ColumnNames{ind})*360/pi;
            elseif contains(hlp,'A-1')
                DataUnitsOpt{7}='A-1';
                hlp0=mytable.(ColumnNames{ind});
                handles.Data.fwhm=4*hlp0*lam/(4*pi)./cosd(handles.Data.tth/2);
            else
                DataUnitsOpt{7}='uknown';
                handles.Data.fwhm=mytable.(ColumnNames{ind})*NaN;
                b(7)=0;
            end
        end
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Beta','Beta_','Beta,','Beta;',...
        'beta','beta_','beta,','beta;',...
        'BETA','BETA_','BETA,','BETA;',...
        'Broadening','Broadening_','Broadening,','Broadening;',...
        'broadening','broadening_','broadening,','broadening;',...
        'BROADENING','BROADENING_','BROADENING,','BROADENING;',...
        'Beta (deg)','Beta (deg)_','Beta (deg),','Beta (deg);',...
        'beta (deg)','beta (deg)_','beta (deg),','beta (deg);',...
        'BETA (deg)','BETA (deg)_','BETA (deg),','BETA (deg);',...
        'BETA (DEG)','BETA (DEG)_','BETA (DEG),','BETA (DEG);',...
        'Broadening (deg)','Broadening (deg)_','Broadening (deg),','Broadening (deg);',...
        'broadening (deg)','broadening (deg)_','broadening (deg),','broadening (deg);',...
        'BROADENING (deg)','BROADENING (deg)_','BROADENING (deg),','BROADENING (deg);',...
        'BROADENING (DEG)','BROADENING (DEG)_','BROADENING (DEG),','BROADENING (DEG);'})==1);
    if ~isempty(ind)
        b(8)=1;
        col(8)=ind;
        handles.Data.beta=mytable.(ColumnNames{ind});

        hlp=lower(mytable.Properties.VariableDescriptions{ind});
        if contains(hlp,'deg')
            DataUnitsOpt{8}='deg';
            handles.Data.beta=mytable.(ColumnNames{ind});
        elseif contains(hlp,'rad')
            DataUnitsOpt{8}='rad';
            handles.Data.beta=mytable.(ColumnNames{ind})*180/pi;
        elseif contains(hlp,'A-1')
            DataUnitsOpt{8}='A-1';
            hlp0=mytable.(ColumnNames{ind});
            handles.Data.beta=2*hlp0*lam/(4*pi)./cosd(handles.Data.tth/2);
        else
            DataUnitsOpt{8}='uknown';
            handles.Data.beta=mytable.(ColumnNames{ind})*NaN;
            b(8)=0;
        end
    end

    ind=find(ismember(mytable.Properties.VariableDescriptions,{'Label','Label_','Label,','Label;',...
        'label','label_','label,','label;',...
        'LABEL','LABEL_','LABEL,','LABEL;',...
        'Comment','Comment_','Comment,','Comment;',...
        'comment','comment_','comment,','comment;',...
        'COMMENT','COMMENT_','COMMENT,','COMMENT;',...
        'Grain','Grain_','Grain,','Grain;',...
        'grain','grain_','grain,','grain;',...
        'GRAIN','GRAIN_','GRAIN,','GRAIN;' ...
        'Domain','Domain_','Domain,','Domain;',...
        'domain','domain_','domain,','domain;',...
        'DOMAIN','DOMAIN_','DOMAIN,','DOMAIN;'})==1);
    if ~isempty(ind)
        b(9)=1;
        col(9)=ind;
        handles.Data.label=mytable.(ColumnNames{ind});
    end

    handles.DataUnitsOpt=DataUnitsOpt;

    if sum(b(1:8))<7
        if sum(b(1:8))==0
            warndlg('No columns have been automaticaly found. Check the file format, file structure or try to import the data with "Manual selection" tool.');
        else
            warndlg('Some columns have not been automaticaly found. Try to re-import the data with "Manual selection" tool.');
        end
        N=size(mytable,1);

        if b(1)==0
            handles.Data.h=NaN*ones(N,1);
        end
        if b(2)==0
            handles.Data.k=NaN*ones(N,1);
        end
        if b(3)==0
            handles.Data.l=NaN*ones(N,1);
        end
        if b(4)==0
            handles.Data.psi=NaN*ones(N,1);
        end
        if b(5)==0
            handles.Data.phi=NaN*ones(N,1);
        end
        if b(6)==0
            handles.Data.tth=NaN*ones(N,1);
        end
        if b(7)==0
            handles.Data.fwhm=NaN*ones(N,1);
        end
        if b(8)==0
            handles.Data.beta=NaN*ones(N,1);
        end
        if b(9)==0
            handles.Data.label=cell(N,1);
        end

        handles.Filename=FileName;
        handles.InputTable=mytable;
        handles.DataPointsNumber=N;
        handles.DataAvailability=b;
        handles.DataColumnN=col;
        handles.DataFormula=formula;
        handles.DataGroup=ones(N,1);
        handles.DataGroupMaxNumber=0;
        handles.Groups_list.String={'All peaks'};
        handles.Groups_list.Value=1;
        handles.FitGroups_chbx.Value=0;
        handles.FitGroups_chbx.Enable='off';
    else
        N=numel(handles.Data.h);
        handles.Filename=FileName;
        handles.InputTable=mytable;
        handles.DataPointsNumber=N;
        handles.DataAvailability=b;
        handles.DataColumnN=col;
        handles.DataFormula=formula;
        handles.DataGroup=ones(N,1);
        handles.DataGroupMaxNumber=0;
        handles.Groups_list.String='All peaks';
        handles.Groups_list.Value=1;
        handles.FitGroups_chbx.Value=0;
        handles.FitGroups_chbx.Enable='off';
    end

    if sum(handles.DataAvailability(1:3))==3 % h, k, l are available
        indskip=find(handles.Data.h==0 & handles.Data.k==0 & handles.Data.l==0);
        if ~isempty(indskip)
            answer = questdlg('There are diffractions with the indices h k l = 0 0 0. Do you want to skip them?','Diffraction indices 0 0 0','Yes','No','Yes');
            if strcmp(answer,'Yes')
                indstay=setdiff(1:N,indskip);
                N=numel(indstay);

                mytable(indstay,:)
                handles.InputTable=mytable(indstay,:);
                handles.DataPointsNumber=N;
                handles.Data.h=handles.Data.h(indstay);
                handles.Data.k=handles.Data.k(indstay);
                handles.Data.l=handles.Data.l(indstay);
                handles.Data.tth=handles.Data.tth(indstay);
                handles.Data.psi=handles.Data.psi(indstay);
                handles.Data.phi=handles.Data.phi(indstay);
                handles.Data.fwhm=handles.Data.fwhm(indstay);
                handles.Data.beta=handles.Data.beta(indstay);
                handles.Data.label=handles.Data.label(indstay);
                handles.DataGroup=ones(N,1);
            end
        end
    end

    TAB=[mat2cell(handles.Data.h,ones(N,1)) mat2cell(handles.Data.k,ones(N,1)) mat2cell(handles.Data.l,ones(N,1)) mat2cell(handles.Data.tth,ones(N,1)) mat2cell(handles.Data.psi,ones(N,1)) mat2cell(handles.Data.phi,ones(N,1)) mat2cell(handles.Data.fwhm,ones(N,1)) mat2cell(handles.Data.beta,ones(N,1)) handles.Data.label mat2cell(true(N,1),ones(N,1))];
    set(handles.DataTable_tab,'Data',TAB);
    handles.ManualSelection_btn.Enable='on';
    handles.AdvancedShiftOpt_btn.Enable='on';
    title(handles.PeakShift_ax,[handles.Filename ': Peak shift evaluation'],'interpreter','none')
    cla(handles.PeakShift_ax);
    title(handles.PeakBroadening_ax,[handles.Filename ': Peak broadening evaluation'],'interpreter','none')
    cla(handles.PeakBroadening_ax);
    handles.Status_list.String={};

    guidata(hObject, handles);
end


% --- Executes on selection change in Groups_list.
function Groups_list_Callback(hObject, eventdata, handles)
% hObject    handle to Groups_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Groups_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Groups_list
I=hObject.Value;
N=handles.DataPointsNumber;
if I==1
    % set checkboxes in DataTable uneditable and all checked (it
    % corresponds to the fact that the first column in the
    % handles.DataGroups is allways equal to ones
    handles.DataTable_tab.ColumnEditable(end)=0;
    for n=1:N
        handles.DataTable_tab.Data(n,end)={true};
    end
    handles.QuickSelector_btn.Enable='off';
else
    handles.DataTable_tab.ColumnEditable(end)=1;
    for n=1:N
        handles.DataTable_tab.Data(n,end)={logical(handles.DataGroup(n,I))}; % there might be better option than doing a for cycle, but I haven't found it
    end
    handles.QuickSelector_btn.Enable='on';
end
% now we need to fill the information about the models of the current group
% Label
if I==1
    handles.GroupModelSelection_txt.String='Group #0: All peaks';
else
    handles.GroupModelSelection_txt.String=['Group #' num2str(I-1)];
end

% Material information
switch handles.GroupsModels(I).sys
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'cubic_hexagonal'
        handles.System_pop.Value=3;
end

if handles.GroupsModels(I).Material.av==0
    handles.Uknown_rad.Value=1;
    handles.MaterialDatabase_list.Enable='off';
    handles.MaterialConstants_tab.Enable='off';
    handles.Stiffness1_txt.Enable='off';
    handles.Stiffness_tab.Enable='off';
    handles.Stiffness2_txt.Enable='off';
else
    switch handles.GroupsModels(I).Material.Option
        case 'E_nu_G'
            handles.YoungModulus_rad.Value=1;
            handles.MaterialDatabase_list.Enable='on';
            handles.MaterialConstants_tab.Enable='on';
            handles.Stiffness1_txt.Enable='off';
            handles.Stiffness_tab.Enable='off';
            handles.Stiffness2_txt.Enable='off';
            
            % select correct database
            handles.MaterialDatabase_txt.String=['Material database (' char(957) ',E,G)'];
            handles.MaterialDatabase_list.String=handles.MaterialConstants_database.Label;
            handles.MaterialDatabase_list.Value=handles.GroupsModels(I).Material.DatabaseRow;

            %vložení dat do tabulky
            tab=handles.MaterialConstants_tab.Data; % to have the proper sizes
            tab{1,1}=handles.GroupsModels(I).Material.nu; % Young modulus (GPa)
            tab{1,2}=handles.GroupsModels(I).Material.E; % Poisson ration
            tab{1,3}=handles.GroupsModels(I).Material.G; % Shear modulus (GPa)
            tab{1,4}=handles.GroupsModels(I).Material.nu_100; % Young modulus (GPa)
            tab{1,5}=handles.GroupsModels(I).Material.E_100; % Poisson ration
            tab{1,6}=handles.GroupsModels(I).Material.G_100; % Shear modulus (GPa)
            handles.MaterialConstants_tab.Data=tab;
        case 'Sij'
            handles.SijMat_rad.Value=1;
            handles.MaterialDatabase_list.Enable='on';
            handles.MaterialConstants_tab.Enable='off';
            handles.Stiffness1_txt.Enable='on';
            handles.Stiffness_tab.Enable='on';
            handles.Stiffness2_txt.Enable='on';

            % select correct database
            handles.MaterialDatabase_txt.String=['Material database (' char(957) ',E,G)'];
            handles.MaterialDatabase_list.String=handles.StiffnessMatrices_database.Label;
            handles.MaterialDatabase_list.Value=handles.GroupsModels(I).Material.DatabaseRow;

            %vložení dat do tabulky
            tab=handles.Stiffness_tab.Data;
            for i=1:6
                for j=1:6
                    tab{i,j}=handles.GroupsModels(I).Material.Sij(i,j);
                end
            end
            handles.Stiffness_tab.Data=tab;
        case 'Cij'
            handles.CijMat_rad.Value=1;
            handles.MaterialDatabase_list.Enable='on';
            handles.MaterialConstants_tab.Enable='off';
            handles.Stiffness1_txt.Enable='on';
            handles.Stiffness_tab.Enable='on';
            handles.Stiffness2_txt.Enable='on';

            % select correct database
            handles.MaterialDatabase_txt.String=['Material database (' char(957) ',E,G)'];
            handles.MaterialDatabase_list.String=handles.StiffnessMatrices_database.Label;
            handles.MaterialDatabase_list.Value=handles.GroupsModels(I).Material.DatabaseRow;

            %vložení dat do tabulky
            tab=handles.Stiffness_tab.Data;
            for i=1:6
                for j=1:6
                    tab{i,j}=handles.GroupsModels(I).Material.Cij(i,j);
                end
            end
            handles.Stiffness_tab.Data=tab;
    end
end
guidata(hObject,handles)

% Peak shift model
switch handles.GroupsModels(I).PeakShiftModel
    case 'NoStress'
        handles.NoStress_rad.Value=1;
    case 'IsotropicStress'
        handles.IsotropicStress_rad.Value=1;
    case 'VookWitt'
        handles.VookWitt_rad.Value=1;
    case 'Reuss'
        handles.Reuss_rad.Value=1;
    case 'AdvancedShiftOpt'
        handles.AdvancedShiftOpt_rad.Value=1;
end

% Peak broadening model
if handles.GroupsModels(I).SizeOptions.Fitted==1
    handles.Size_chbx.Value=1;
else
    handles.Size_chbx.Value=0;
end

if handles.GroupsModels(I).MicrostrainOptions.Fitted==1
    handles.Microstrain_chbx.Value=1;
else
    handles.Microstrain_chbx.Value=0;
end

% switch handles.GroupsModels(I).MicrostrainOptions.PeakBroadeningModel
%     case 'Isotropic'
%         handles.IsotropicMicrostrain_rad.Value=1;
%     case 'Dislocations'
%         handles.Dislocations_rad.Value=1;
% end

if handles.GroupsModels(I).FitSF==1
    handles.SF_chbx.Value=1;
else
    handles.SF_chbx.Value=0;
end

handles.GroupSelected=I;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Groups_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Groups_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Shift_chbx.
function Shift_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to Shift_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Shift_chbx


% --- Executes on button press in Uknown_rad.
function Uknown_rad_Callback(hObject, eventdata, handles)
% hObject    handle to Uknown_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Uknown_rad
handles.MaterialConstants_tab.Enable='off';
handles.Stiffness_tab.Enable='off';
handles.Stiffness1_txt.Enable='off';
handles.Stiffness2_txt.Enable='off';
handles.MaterialDatabase_list.Enable='off';
handles.Reference_txt.Enable='off';

I=handles.GroupSelected;
handles.GroupsModels(I).Material.av=0;
guidata(hObject, handles);


% --- Executes on button press in YoungModulus_rad.
function YoungModulus_rad_Callback(hObject, eventdata, handles)
% hObject    handle to YoungModulus_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of YoungModulus_rad
handles.MaterialConstants_tab.Enable='on';
handles.Stiffness_tab.Enable='off';
handles.Stiffness1_txt.Enable='off';
handles.Stiffness2_txt.Enable='off';
handles.MaterialDatabase_list.Enable='on';
handles.MaterialDatabase_txt.String=['Material database (' char(957) ',E,G)'];
handles.MaterialDatabase_list.String=handles.MaterialConstants_database.Label;
handles.MaterialDatabase_list.Value=1;
switch handles.MaterialDatabase_list.String{1}
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'hexagonal'
        handles.System_pop.Value=3;
end
MaterialDatabase_list_Callback(handles.MaterialDatabase_list, eventdata, handles);

I=handles.GroupSelected;
handles.GroupsModels(I).Material.DatabaseRow=1;
handles.GroupsModels(I).Material.av=1;
handles.GroupsModels(I).Material.Option='E_nu_G';
tab=handles.MaterialConstants_tab.Data;
handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)

guidata(hObject, handles);


% --- Executes on button press in SijMat_rad.
function SijMat_rad_Callback(hObject, eventdata, handles)
% hObject    handle to SijMat_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SijMat_rad
handles.MaterialConstants_tab.Enable='off';
handles.Stiffness_tab.Enable='on';
handles.Stiffness1_txt.Enable='on';
handles.Stiffness1_txt.String='S=(';
handles.Stiffness2_txt.Enable='on';
handles.MaterialDatabase_list.Enable='on';
handles.MaterialDatabase_txt.String='Material database (Sij/Cij)';
handles.MaterialDatabase_list.String=handles.StiffnessMatrices_database.Label;
handles.MaterialDatabase_list.Value=1;
switch handles.MaterialDatabase_list.String{1}
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'hexagonal'
        handles.System_pop.Value=3;
end
MaterialDatabase_list_Callback(handles.MaterialDatabase_list, eventdata, handles);

I=handles.GroupSelected;
handles.GroupsModels(I).Material.DatabaseRow=1;
handles.GroupsModels(I).Material.av=1;
handles.GroupsModels(I).Material.Option='Sij';
tab=handles.MaterialConstants_tab.Data;
handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.Sij=cell2mat(handles.Stiffness_tab.Data);
handles.GroupsModels(I).Material.Cij=inv(handles.GroupsModels(I).Material.Sij);

guidata(hObject, handles);


% --- Executes on button press in ManualSelection_btn.
function ManualSelection_btn_Callback(hObject, eventdata, handles)
% hObject    handle to ManualSelection_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InputData.Filename=handles.Filename;
InputData.InputTable=handles.InputTable;
InputData.InputData=handles.Data;
InputData.Wavelength=handles.Wavelength;
InputData.DataPointsNumber=handles.DataPointsNumber;
InputData.DataAvailability=handles.DataAvailability;
InputData.DataColumnN=handles.DataColumnN;
InputData.DataFormula=handles.DataFormula;
InputData.DataUnitsOpt=handles.DataUnitsOpt;
[OutputData,Wavelength,DataAvailability,DataColumnN,DataFormula,DataUnitsOpt]=DataManualSelection(InputData);
% set new values in handles and in table
handles.Data=OutputData;
handles.Wavelength=Wavelength;
handles.DataAvailability=DataAvailability;
handles.DataColumnN=DataColumnN;
handles.DataFormula=DataFormula;
handles.DataUnitsOpt=DataUnitsOpt;

N=handles.DataPointsNumber;
TAB=[mat2cell(handles.Data.h,ones(N,1)) mat2cell(handles.Data.k,ones(N,1)) mat2cell(handles.Data.l,ones(N,1)) mat2cell(handles.Data.tth,ones(N,1)) mat2cell(handles.Data.psi,ones(N,1)) mat2cell(handles.Data.phi,ones(N,1)) mat2cell(handles.Data.fwhm,ones(N,1)) mat2cell(handles.Data.beta,ones(N,1)) handles.Data.label mat2cell(true(N,1),ones(N,1))];
set(handles.DataTable_tab,'Data',TAB);
guidata(hObject, handles);

Groups_list_Callback(handles.Groups_list, eventdata, handles)


% --- Executes on button press in NewGroup_btn.
function NewGroup_btn_Callback(hObject, eventdata, handles)
% hObject    handle to NewGroup_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataGroup=[handles.DataGroup ones(handles.DataPointsNumber,1)];
handles.DataGroupMaxNumber=handles.DataGroupMaxNumber+1;
handles.Groups_list.String{end+1}=['Group #' num2str(handles.DataGroupMaxNumber)];
% Default options for information about the material
handles.GroupsModels(handles.DataGroupMaxNumber+1).sys=handles.MaterialConstants_database.System{1};
handles.GroupsModels(handles.DataGroupMaxNumber+1).Material.av=0;
% set the default models for the new group datapoints evaluation
handles.GroupsModels(handles.DataGroupMaxNumber+1).PeakShiftModel='IsotropicStress';
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrainsShape='Sphere';
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.Fitted=1;
% handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.PeakBroadeningModel='Isotropic';
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.MicrostrainModel='Isotropic';
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.ParGuess=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.ParFix=false;
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.ParLB=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.ParUB=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).MicrostrainOptions.Fitted=1;
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.ParGuess=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.ParFix=false;
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.ParLB=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.ParUB=NaN;
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrownDirectionOpt='Sample normal';
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrownDirectionHKL=[0 0 1];
% handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrownDirectionOpt2='Direction in sample surface specified by azimuth';
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrownDirectionHKL2=[1 0 0];
handles.GroupsModels(handles.DataGroupMaxNumber+1).SizeOptions.GrownDirectionPerpAzimuth=0;
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.MaterialInfo='XEC';
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.XECmodel='Voigt';
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.KroenerCoeff=0.5;
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.OrientationAngles={0 0 0 'z' 'x' 'z' 'extrinsic'};
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.DatagroupsIndividual=0;
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.StressFactorsMethod='Reuss';
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.ODFtype='IdealPolycrystal';
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.ODFfilename='';
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.FitN=1;
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.StressComponents={'p(1)'; 'p(1)'; '0'; '0'; '0'; '0'};
handles.GroupsModels(handles.DataGroupMaxNumber+1).AdvancedShiftOptions.Geometry='GAXRD';
handles.GroupsModels(handles.DataGroupMaxNumber+1).FitSF=0;
handles.FitGroups_chbx.Enable='on';
guidata(hObject,handles);


% --- Executes on button press in RemoveGroup_btn.
function RemoveGroup_btn_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveGroup_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.Groups_list.Value;
if I==1
    % One cannot remove the default group. This condition solves also
    % removal of all groups from the list.
    errordlg('You cannot remove the default ''All peaks'' group.','Unallowed group removal');
else
    N0=size(handles.Groups_list.String,1); % original number of groups
    handles.DataGroup=[handles.DataGroup(:,1:I-1) handles.DataGroup(:,I+1:end)];
    handles.Groups_list.String=[handles.Groups_list.String(1:I-1); handles.Groups_list.String(I+1:end)];
    % The Value parameter of listbox remains the same, however, if the last
    % element is removed, it will cause error, because the "selection
    % index" will be higher than the number of elements.
    if I==N0
        handles.Groups_list.Value=I-1;
    end
end
if size(handles.DataGroup,2)==1
    handles.FitGroups_chbx.Value=0;
    handles.FitGroups_chbx.Enable='off';
end
guidata(hObject,handles);
Groups_list_Callback(handles.Groups_list, eventdata, handles)


function GlobOpt=GetGlobalOptions(handles)
% Vertical sample shift
if handles.VerticalSampleShift_chbx.Value==1
    GlobOpt.FitVerticalShift=1;
else
    GlobOpt.FitVerticalShift=0;
end
% Detector shift
if handles.DetectorShift_chbx.Value==1
    GlobOpt.FitDetectorShift=1;
else
    GlobOpt.FitDetectorShift=0;
end
GlobOpt.lam=handles.Wavelength;
GlobOpt.FitGroups=handles.FitGroups_chbx.Value;


function FitOpt=GetFitOptions(handles,I)
% Stress model
FitOpt.PeakShiftModel=handles.GroupsModels(I).PeakShiftModel;
% Stress model - advanced options
FitOpt.AdvancedShiftOptions=handles.GroupsModels(I).AdvancedShiftOptions;
% Broadening model
% Size
if handles.GroupsModels(I).SizeOptions.Fitted==1
    FitOpt.FitSize=1;
    FitOpt.GrainsShape=handles.GroupsModels(I).SizeOptions.GrainsShape;
    FitOpt.SizeOptions.Fitted=handles.GroupsModels(I).SizeOptions.Fitted;
    FitOpt.SizeOptions.ParGuess=handles.GroupsModels(I).SizeOptions.ParGuess;
    FitOpt.SizeOptions.ParFix=handles.GroupsModels(I).SizeOptions.ParFix;
    FitOpt.SizeOptions.ParLB=handles.GroupsModels(I).SizeOptions.ParLB;
    FitOpt.SizeOptions.ParUB=handles.GroupsModels(I).SizeOptions.ParUB;
    FitOpt.SizeOptions.GrownDirectionOpt=handles.GroupsModels(I).SizeOptions.GrownDirectionOpt;
    FitOpt.SizeOptions.GrownDirectionHKL=handles.GroupsModels(I).SizeOptions.GrownDirectionHKL;
%     FitOpt.SizeOptions.GrownDirectionOpt2=handles.GroupsModels(I).SizeOptions.GrownDirectionOpt2;
    FitOpt.SizeOptions.GrownDirectionHKL2=handles.GroupsModels(I).SizeOptions.GrownDirectionHKL2;
    FitOpt.SizeOptions.GrownDirectionPerpAzimuth=handles.GroupsModels(I).SizeOptions.GrownDirectionPerpAzimuth;
else
    FitOpt.FitSize=0;
end
% Microstrain
if handles.GroupsModels(I).MicrostrainOptions.Fitted==1
    FitOpt.FitMicrostrain=1;
%     FitOpt.MicrostrainOptions.PeakBroadeningModel=handles.GroupsModels(I).MicrostrainOptions.PeakBroadeningModel;
    FitOpt.MicrostrainOptions.MicrostrainModel=handles.GroupsModels(I).MicrostrainOptions.MicrostrainModel;
    FitOpt.MicrostrainOptions.Fitted=handles.GroupsModels(I).MicrostrainOptions.Fitted;
    FitOpt.MicrostrainOptions.ParGuess=handles.GroupsModels(I).MicrostrainOptions.ParGuess;
    FitOpt.MicrostrainOptions.ParFix=handles.GroupsModels(I).MicrostrainOptions.ParFix;
    FitOpt.MicrostrainOptions.ParLB=handles.GroupsModels(I).MicrostrainOptions.ParLB;
    FitOpt.MicrostrainOptions.ParUB=handles.GroupsModels(I).MicrostrainOptions.ParUB;
%     FitOpt.MicrostrainOptions.Fitted=handles.GroupsModels(I).MicrostrainOptions.Fitted;
else
    FitOpt.MicrostrainOptions.MicrostrainModel='none';
    FitOpt.FitMicrostrain=0;
end
% Stacking faults
if handles.GroupsModels(I).FitSF==1
    FitOpt.FitSF=1;
else
    FitOpt.FitSF=0;
end
% Crystallographic system XXX: společně s materiálem taky možno
% modifikovat, aby bylo různé pro jednotlivé datasety
FitOpt.sys=handles.System_pop.String{handles.System_pop.Value};
% Material information
if handles.GroupsModels(I).Material.av==0
    FitOpt.Material.av=0;
else
    % both second options should work, because the E, nu, G is computed
    % everytime the S_ij matrix is changed
    FitOpt.Material.av=1;
    FitOpt.Material.Option=handles.GroupsModels(I).Material.Option;
    % here, the data from the MaterialConstants_tab should be loaded
    switch handles.GroupsModels(I).Material.Option
        case 'E_nu_G'
            FitOpt.Material.nu=handles.GroupsModels(I).Material.nu; % Young modulus (GPa)
            FitOpt.Material.E=handles.GroupsModels(I).Material.E; % Poisson ration
            FitOpt.Material.G=handles.GroupsModels(I).Material.G; % Shear modulus (GPa)
            FitOpt.Material.nu_100=handles.GroupsModels(I).Material.nu_100; % Young modulus (GPa)
            FitOpt.Material.E_100=handles.GroupsModels(I).Material.E_100; % Poisson ration
            FitOpt.Material.G_100=handles.GroupsModels(I).Material.G_100; % Shear modulus (GPa)
        case 'Sij'
            FitOpt.Material.nu=handles.GroupsModels(I).Material.nu; % Young modulus (GPa)
            FitOpt.Material.E=handles.GroupsModels(I).Material.E; % Poisson ration
            FitOpt.Material.G=handles.GroupsModels(I).Material.G; % Shear modulus (GPa)
            FitOpt.Material.nu_100=handles.GroupsModels(I).Material.nu_100; % Young modulus (GPa)
            FitOpt.Material.E_100=handles.GroupsModels(I).Material.E_100; % Poisson ration
            FitOpt.Material.G_100=handles.GroupsModels(I).Material.G_100; % Shear modulus (GPa)
            
            FitOpt.Material.Sij=handles.GroupsModels(I).Material.Sij;
            FitOpt.Material.Cij=inv(FitOpt.Material.Sij);
        case 'Cij'
            FitOpt.Material.nu=handles.GroupsModels(I).Material.nu; % Young modulus (GPa)
            FitOpt.Material.E=handles.GroupsModels(I).Material.E; % Poisson ration
            FitOpt.Material.G=handles.GroupsModels(I).Material.G; % Shear modulus (GPa)
            FitOpt.Material.nu_100=handles.GroupsModels(I).Material.nu_100; % Young modulus (GPa)
            FitOpt.Material.E_100=handles.GroupsModels(I).Material.E_100; % Poisson ration
            FitOpt.Material.G_100=handles.GroupsModels(I).Material.G_100; % Shear modulus (GPa)

            FitOpt.Material.Cij=handles.GroupsModels(I).Material.Cij;
            FitOpt.Material.Sij=inv(FitOpt.Material.Cij);
    end
end
FitOpt.BroadeningNFitOrder=str2num(handles.BroadeningNFitOrder_ed.String); % XXX: může být různé, ale možná lepší stejné pro všechny kvůli kreslení
FitOpt.FitHcpLP=handles.FitHcpLP_chbx.Value; % XXX: V BUDOUCNU UPRAVIT


function BroadeningNFitOrder_ed_Callback(hObject, eventdata, handles)
% hObject    handle to BroadeningNFitOrder_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BroadeningNFitOrder_ed as text
%        str2double(get(hObject,'String')) returns contents of BroadeningNFitOrder_ed as a double


% --- Executes during object creation, after setting all properties.
function BroadeningNFitOrder_ed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BroadeningNFitOrder_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PlotData(handles)
hkl=[handles.Data.h handles.Data.k handles.Data.l];
tth=handles.Data.tth;
myclrlist='rgbcmyk';
% Peak shift
cla(handles.PeakShift_ax)
if handles.GlobOpt.FitGroups==0
    if handles.PeakShiftResults(1).av==1
        sin2psi=(sind(handles.Data.psi)).^2;
        [~,ind]=sort(sin2psi);
%         switch handles.GlobOpt.sys % XXX: should be edited in the future
%             case {'cubic_fcc','cubic_bcc'}
                ahkl=handles.PeakShiftResults(1).ahkl;
                ahkl_T=handles.PeakShiftResults(1).ahkl_T;
                hp=plot(handles.PeakShift_ax,sin2psi(ind),ahkl(ind),'bo','DisplayName','data');
                hp.Tag='Data_PS #1';
                hold(handles.PeakShift_ax,'on')
                plot(handles.PeakShift_ax,sin2psi(ind),ahkl_T(ind),'r','DisplayName','fit');
                hp.Tag='Fit_PS #1';
                hold(handles.PeakShift_ax,'off')
                yrng=max([ahkl;ahkl_T])-min([ahkl;ahkl_T]);
                for n=1:numel(ahkl)
                    text(handles.PeakShift_ax,sin2psi(ind(n)),ahkl(ind(n))+0.03*yrng,[num2str(hkl(ind(n),1)) ' ' num2str(hkl(ind(n),2)) ' ' num2str(hkl(ind(n),3))],'Rotation',90);
                end
                ylabel(handles.PeakShift_ax,'$$a_{hkl}$$ (\AA)','interpreter','latex');
%             case 'hexagonal' % XXX: should be edited in the future
%                 dhkl=handles.PeakShiftResults.ahkl;
%                 dhkl_T=handles.PeakShiftResults.ahkl_T;
%                 B=Reci(handles.PeakShiftResults.results.a0,handles.PeakShiftResults.results.a0,handles.PeakShiftResults.results.c0,90,90,120);
%                 dhkl0=zeros(size(tth));
%                 for n=1:numel(tth)
%                     dhkl0(n)=2*pi./norm(B*hkl(n,:)');
%                 end
%                 epshkl=(dhkl-dhkl0)./dhkl0;
%                 epshkl_T=(dhkl_T-dhkl0)./dhkl0;
%                 hold(handles.PeakShift_ax,'on')
%                 plot(handles.PeakShift_ax,sin2psi(ind),epshkl(ind),'bo','DisplayName','data');
%                 plot(handles.PeakShift_ax,sin2psi(ind),epshkl_T(ind),'r','DisplayName','fit');
%                 hold(handles.PeakShift_ax,'off')
%                 yrng=max([epshkl;epshkl_T])-min([epshkl;epshkl_T]);
%                 for n=1:numel(epshkl)
%                     text(handles.PeakShift_ax,sin2psi(ind(n)),epshkl(ind(n))+0.03*yrng,[num2str(hkl(ind(n),1)) ' ' num2str(hkl(ind(n),2)) ' ' num2str(hkl(ind(n),3))],'Rotation',90);
%                 end
%                 ylabel(handles.PeakShift_ax,'$$\epsilon_{hkl}$$','interpreter','latex');
%         end
    end
else
    for g=2:size(handles.DataGroup,2)
        if handles.PeakShiftResults(g).av==1
            indg=find(handles.DataGroup(:,g)==1);
            sin2psi=(sind(handles.Data.psi(indg))).^2;
            [~,ind]=sort(sin2psi);
            switch handles.GroupsModels(g).sys
                case {'cubic_fcc','cubic_bcc'}
                    ahkl=handles.PeakShiftResults(g).ahkl;
                    ahkl_T=handles.PeakShiftResults(g).ahkl_T;
                    hold(handles.PeakShift_ax,'on')
                    hp=plot(handles.PeakShift_ax,sin2psi(ind),ahkl(ind),'o','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['data G#' num2str(g-1)]);
                    hp.Tag=['Data_PS #' num2str(g)];
                    hp=plot(handles.PeakShift_ax,sin2psi(ind),ahkl_T(ind),'-','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['fit G#' num2str(g-1)]);
                    hp.Tag=['Fit_PS #' num2str(g)];
                    hold(handles.PeakShift_ax,'off')
                    yrng=max([ahkl;ahkl_T])-min([ahkl;ahkl_T]);
                    for n=1:numel(ahkl)
                        text(handles.PeakShift_ax,sin2psi(ind(n)),ahkl(ind(n))+0.03*yrng,[num2str(hkl(indg(ind(n)),1)) ' ' num2str(hkl(indg(ind(n)),2)) ' ' num2str(hkl(indg(ind(n)),3))],'Rotation',90,'Color',myclrlist(mod(g-2,7)+1));
                    end
                    ylabel(handles.PeakShift_ax,'$$a_{hkl}$$ (\AA)','interpreter','latex');
                case 'hexagonal'
                    dhkl=handles.PeakShiftResults(g).ahkl;
                    dhkl_T=handles.PeakShiftResults(g).ahkl_T;
                    B=Reci(handles.PeakShiftResults(g).results.a0,handles.PeakShiftResults(g).results.a0,handles.PeakShiftResults(g).results.c0,90,90,120);
                    dhkl0=zeros(size(tth(indg)));
                    for n=1:numel(tth(indg))
                        dhkl0(n)=2*pi./norm(B*hkl(n,:)');
                    end
                    epshkl=(dhkl-dhkl0)./dhkl0;
                    epshkl_T=(dhkl_T-dhkl0)./dhkl0;
                    hold(handles.PeakShift_ax,'on')
                    plot(handles.PeakShift_ax,sin2psi(ind),epshkl(ind),'o','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['data G#' num2str(g-1)]);
                    plot(handles.PeakShift_ax,sin2psi(ind),epshkl_T(ind),'-','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['fit G#' num2str(g-1)]);
                    hold(handles.PeakShift_ax,'off')
                    yrng=max([epshkl;epshkl_T])-min([epshkl;epshkl_T]);
                    for n=1:numel(epshkl)
                        text(handles.PeakShift_ax,sin2psi(ind(n)),epshkl(ind(n))+0.03*yrng,[num2str(hkl(indg(ind(n)),1)) ' ' num2str(hkl(indg(ind(n)),2)) ' ' num2str(hkl(indg(ind(n)),3))],'Rotation',90,'Color',myclrlist(mod(g-2,7)+1));
                    end
                    ylabel(handles.PeakShift_ax,'$$\epsilon_{hkl}$$','interpreter','latex');
            end
        end
    end
end

xlabel(handles.PeakShift_ax,'$$\sin^2\psi$$','interpreter','latex');
box(handles.PeakShift_ax,'on')
title(handles.PeakShift_ax,[handles.Filename ': Peak shift evaluation'],'interpreter','none')
hL=legend(handles.PeakShift_ax,'show');
if handles.PeakShiftResults(1).av==1 && isfield(handles.PeakShiftResults(1).results,'sig_ef') && handles.PeakShiftResults(1).results.sig_ef>0
    set(hL,'Location','southeast')
else
    set(hL,'Location','northeast')
end

% Peak broadening
sinth=sind(tth/2);
[~,ind]=sort(sinth);
% N=handles.FitOpt.BroadeningNFitOrder; % XXX: this should be set as one value, but in fact could differ for different datasets
N=str2num(handles.BroadeningNFitOrder_ed.String); % XXX: this should be set as one value, but in fact could differ for different datasets
cla(handles.PeakBroadening_ax)
if handles.GlobOpt.FitGroups==0
    % fwhm
    if handles.PeakBroadeningResults_fwhm(1).av==1
        broad=handles.PeakBroadeningResults_fwhm(1).broad;
        broad_T=handles.PeakBroadeningResults_fwhm(1).broad_T;
        hold(handles.PeakBroadening_ax,'on')
        hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad(ind).^N,'bo','DisplayName','fwhm - data');
        hp.Tag='Data_FWHM #1';
        hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad_T(ind).^N,'r','DisplayName','fwhm - fit');
        hp.Tag='Data_FWHM #1';
        hold(handles.PeakBroadening_ax,'off')
        yrng=max([broad(ind).^N;broad_T(ind).^N])-min([broad(ind).^N;broad_T(ind).^N]);
        for n=1:numel(sinth)
            text(handles.PeakBroadening_ax,sinth(ind(n)).^N,broad(ind(n)).^N+0.03*yrng,[num2str(hkl(ind(n),1)) ' ' num2str(hkl(ind(n),2)) ' ' num2str(hkl(ind(n),3))],'Rotation',90);
        end
    end
    % beta
    if handles.PeakBroadeningResults_beta(1).av==1
        broad=handles.PeakBroadeningResults_beta(1).broad;
        broad_T=handles.PeakBroadeningResults_beta(1).broad_T;
        hold(handles.PeakBroadening_ax,'on')
        hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad(ind).^N,'o','Color',[0 0.333 0],'DisplayName','\beta - data');
        hp.Tag='Data_beta #1';
        hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad_T(ind).^N,'g','DisplayName','\beta - fit');
        hp.Tag='Data_beta #1';
        hold(handles.PeakBroadening_ax,'off')
        yrng=max([broad(ind).^N;broad_T(ind).^N])-min([broad(ind).^N;broad_T(ind).^N]);
        for n=1:numel(sinth)
            text(handles.PeakBroadening_ax,sinth(ind(n)).^N,broad(ind(n)).^N+0.03*yrng,[num2str(hkl(ind(n),1)) ' ' num2str(hkl(ind(n),2)) ' ' num2str(hkl(ind(n),3))],'Rotation',90);
        end
    end
else
    for g=2:size(handles.DataGroup,2)
        if handles.PeakBroadeningResults_fwhm(g).av==1
            indg=find(handles.DataGroup(:,g)==1);
            sinth=sind(tth(indg)/2);
            [~,ind]=sort(sinth);
            broad=handles.PeakBroadeningResults_fwhm(g).broad;
            broad_T=handles.PeakBroadeningResults_fwhm(g).broad_T;
            hold(handles.PeakBroadening_ax,'on')
            hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad(ind).^N,'o','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['fwhm - data G#' num2str(g-1)]);
            hp.Tag=['Data_FWHM #' num2str(g)];
            hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad_T(ind).^N,'-','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['fwhm - fit G#' num2str(g-1)]);
            hp.Tag=['Data_FWHM #' num2str(g)];
            hold(handles.PeakBroadening_ax,'off')
            yrng=max([broad(ind).^N;broad_T(ind).^N])-min([broad(ind).^N;broad_T(ind).^N]);
            for n=1:numel(sinth)
                text(handles.PeakBroadening_ax,sinth(ind(n)).^N,broad(ind(n)).^N+0.03*yrng,[num2str(hkl(indg(ind(n)),1)) ' ' num2str(hkl(indg(ind(n)),2)) ' ' num2str(hkl(indg(ind(n)),3))],'Rotation',90,'Color',myclrlist(mod(g-2,7)+1));
            end
        end
        % beta
        if handles.PeakBroadeningResults_beta(g).av==1
            indg=find(handles.DataGroup(:,g)==1);
            sinth=sind(tth(indg)/2);
            [~,ind]=sort(sinth);
            broad=handles.PeakBroadeningResults_beta(g).broad;
            broad_T=handles.PeakBroadeningResults_beta(g).broad_T;
            hold(handles.PeakBroadening_ax,'on')
            hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad(ind).^N,'*','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['\beta - data G#' num2str(g-1)]);
            hp.Tag=['Data_beta #' num2str(g)];
            hp=plot(handles.PeakBroadening_ax,sinth(ind).^N,broad_T(ind).^N,':','Color',myclrlist(mod(g-2,7)+1),'DisplayName',['\beta - fit G#' num2str(g-1)]);
            hp.Tag=['Data_beta #' num2str(g)];
            hold(handles.PeakBroadening_ax,'off')
            yrng=max([broad(ind).^N;broad_T(ind).^N])-min([broad(ind).^N;broad_T(ind).^N]);
            for n=1:numel(sinth)
                text(handles.PeakBroadening_ax,sinth(ind(n)).^N,broad(ind(n)).^N+0.03*yrng,[num2str(hkl(indg(ind(n)),1)) ' ' num2str(hkl(indg(ind(n)),2)) ' ' num2str(hkl(indg(ind(n)),3))],'Rotation',90,'Color',myclrlist(mod(g-2,7)+1));
            end
        end
    end
end
if N==1
    xlabel(handles.PeakBroadening_ax,'$$\sin\theta$$','interpreter','latex');
    ylabel(handles.PeakBroadening_ax,'$$FWHM$$ $$\left(\textrm{\AA}^{-1}\right)$$, $$\beta$$ $$\left(\textrm{\AA}^{-1}\right)$$','interpreter','latex');
else
    xlabel(handles.PeakBroadening_ax,['$$\sin^{' num2str(N) '}\theta$$'],'interpreter','latex');
    ylabel(handles.PeakBroadening_ax,['$$FWHM^{' num2str(N) '}$$ $$\left(\textrm{\AA}^{' num2str(-N) '}\right)$$, $$\beta^{' num2str(N) '}$$ $$\left(\textrm{\AA}^{' num2str(-N) '}\right)$$'],'interpreter','latex');
end
box(handles.PeakBroadening_ax,'on')
title(handles.PeakBroadening_ax,[handles.Filename ': Peak broadening evaluation'],'interpreter','none')
hL=legend(handles.PeakBroadening_ax,'show');
set(hL,'Location','southeast')

dcm_obj = datacursormode(handles.figure1);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,handles.Data,handles.DataGroup});


function str=DispResults(handles,GlobOpt)
% Because each datapoints group can have different model, the results
% structure field will look differently. The help variables resultsPS,
% resultsPB_fwhm and resultsPB_beta have to be alocated directly during the
% for g=1:... cycles
str={};
fprintf('Results:\n');
str=[str; 'Peak shift:'];
fprintf('Peak shift:\n');

for g=1:size(handles.DataGroup,2)
    if (GlobOpt.FitGroups==0 & g==1) | (GlobOpt.FitGroups==1 & g>1)
        str=[str; '----------------------------------'];
        str=[str; 'GROUP #' num2str(g-1) ':'];
        str=[str; 'Selected peaks: ' sprintf('%d ',handles.DataGroup(:,g))];
        str=[str; 'Peak shift:'];
        if strcmp(handles.GroupsModels(g).PeakShiftModel,'AdvancedShiftOpt')
            if strcmp(handles.GroupsModels(g).AdvancedShiftOptions.MaterialInfo,'XEC')
                str=[str; 'Model: AdvancedShiftOpt - XECs (' handles.GroupsModels(g).AdvancedShiftOptions.XECmodel ' model)'];
            else
                str=[str; 'Model: AdvancedShiftOpt - ' handles.GroupsModels(g).AdvancedShiftOptions.MaterialInfo];
            end
        else
            str=[str; 'Model: ' handles.GroupsModels(g).PeakShiftModel];
        end
        fprintf('----------------------------------\n');
        fprintf(['GROUP #' num2str(g-1) ':\n']);
        fprintf('Selected peaks: ');
        fprintf('%d ',handles.DataGroup(:,g));
        fprintf('\n');
        fprintf('Peak shift:\n');
        if strcmp(handles.GroupsModels(g).PeakShiftModel,'AdvancedShiftOpt')
            if strcmp(handles.GroupsModels(g).AdvancedShiftOptions.MaterialInfo,'XEC')
                fprintf(['Model: AdvancedShiftOpt - XECs (' handles.GroupsModels(g).AdvancedShiftOptions.XECmodel ' model)\n']);
            else
                fprintf(['Model: AdvancedShiftOpt - ' handles.GroupsModels(g).AdvancedShiftOptions.MaterialInfo '\n']);
            end
        else
            fprintf(['Model: ' handles.GroupsModels(g).PeakShiftModel '\n']);
        end

        if handles.PeakShiftResults(g).av==1
            resultsPS=handles.PeakShiftResults(g).results;
            if strcmp(handles.GroupsModels(g).PeakShiftModel,'AdvancedShiftOpt')
                str=[str; 'Lattice parameter:'];
                str=[str; sprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
                
                str=[str; ' '];
                str=[str; 'Computed stress components:'];
                str=[str; [char(963) char(8321) char(8321) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{1}]];
                str=[str; [char(963) char(8322) char(8322) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{2}]];
                str=[str; [char(963) char(8323) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{3}]];
                str=[str; [char(963) char(8322) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{4}]];
                str=[str; [char(963) char(8321) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{5}]];
                str=[str; [char(963) char(8321) char(8322) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{6}]];
                for n=1:handles.GroupsModels(g).AdvancedShiftOptions.FitN
                    str=[str; sprintf(['p(' num2str(n) ') = %.4f ' char(177) ' %.4f'],resultsPS.sigP(n),resultsPS.sigP_err(n))];
                end
                str=[str; sprintf('Note: p(n) has such units which will produce [GPa] in stress components.')];

                str=[str; ' '];
                str=[str; 'Aberations:'];
                if GlobOpt.FitVerticalShift==1
                    str=[str; 'Sample vertical shift:'];
                    if handles.GroupsModels(g).Material.av==1
                        str=[str; sprintf(['sycos = %.4f ' char(177) ' %.4f'],resultsPS.sycos,resultsPS.sycos_err)];
                    else
                        str=[str; sprintf(['sycos = %.4f ' char(177) ' %.4f'],resultsPS.sycos_ef,resultsPS.sycos_ef_err)];
                    end
                else
                    str=[str; 'Sample vertical shift: not fitted'];
                end
                if GlobOpt.FitDetectorShift==1
                    str=[str; 'Detector zero shift:'];
                    if handles.GroupsModels(g).Material.av==1
                        str=[str; sprintf(['zero = %.4f ' char(177) ' %.4f deg'],resultsPS.dth*180/pi,resultsPS.dth_err*180/pi)];
                    else
                        str=[str; sprintf(['zero = %.4f ' char(177) ' %.4f'],resultsPS.dth_ef,resultsPS.dth_ef_err)];
                    end
                else
                    str=[str; 'Detector zero shift: not fitted'];
                end

                fprintf('Lattice parameter:\n');
                fprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);

                fprintf('Computed stress components:\n');
                fprintf([char(963) char(8321) char(8321) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{1} '\n']);
                fprintf([char(963) char(8322) char(8322) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{2} '\n']);
                fprintf([char(963) char(8323) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{3} '\n']);
                fprintf([char(963) char(8322) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{4} '\n']);
                fprintf([char(963) char(8321) char(8323) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{5} '\n']);
                fprintf([char(963) char(8321) char(8322) ' = ' handles.GroupsModels(g).AdvancedShiftOptions.StressComponents{6} '\n']);
                for n=1:handles.GroupsModels(g).AdvancedShiftOptions.FitN
                    fprintf(['p(' num2str(n) ') = %.4f ' char(177) ' %.4f\n'],resultsPS.sigP(n),resultsPS.sigP_err(n));
                end
                fprintf('Note: p(n) has such units which will produce [GPa] in stress components.\n');
                
                fprintf('Aberations:\n');
                if GlobOpt.FitVerticalShift==1
                    fprintf('Sample vertical shift:\n');
                    if handles.GroupsModels(g).Material.av==1
                        fprintf(['sycos = %.4f ' char(177) ' %.4f\n'],resultsPS.sycos,resultsPS.sycos_err);
                    else
                        fprintf(['sycos = %.4f ' char(177) ' %.4f\n'],resultsPS.sycos_ef,resultsPS.sycos_ef_err);
                    end
                else
                    fprintf('Sample vertical shift: not fitted\n');
                end
                if GlobOpt.FitDetectorShift==1
                    fprintf('Detector zero shift:\n');
                    if handles.GroupsModels(g).Material.av==1
                        fprintf(['zero = %.4f ' char(177) ' %.4f deg\n'],resultsPS.dth*180/pi,resultsPS.dth_err*180/pi);
                    else
                        fprintf(['zero = %.4f ' char(177) ' %.4f\n'],resultsPS.dth_ef,resultsPS.dth_ef_err);
                    end
                else
                    fprintf('Detector zero shift: not fitted\n');
                end
        else
            switch handles.GroupsModels(g).sys
                case {'cubic_fcc','cubic_bcc'}
                    if (strcmp(handles.GroupsModels(g).PeakShiftModel,'VookWitt') | strcmp(handles.GroupsModels(g).PeakShiftModel,'Reuss')) && handles.GroupsModels(g).Material.av==1 && isnan(handles.GroupsModels(g).Material.E_100)
                        str=[str; [char(957) char(8321) char(8320) char(8320) ' and E' char(8321) char(8320) char(8320) ' are not available, the tabled averages were used instead.']];
                    end
                    if handles.GroupsModels(g).Material.av==1 && ~strcmp(handles.GroupsModels(g).PeakShiftModel,'NoStress')
                        str=[str; sprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
                        str=[str; sprintf([char(963) ' = %.4f ' char(177) ' %.4f GPa'],resultsPS.sig,resultsPS.sig_err)];

                        fprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);
                        fprintf([char(963) ' = %.4f ' char(177) ' %.4f GPa\n'],resultsPS.sig,resultsPS.sig_err);
                    else
                        if strcmp(handles.GroupsModels(g).PeakShiftModel,'NoStress')
                            str=[str; sprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
   
                            fprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);
                        else
                            str=[str; sprintf(['a' char(8741) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a_par,resultsPS.a_par_err)];
                            str=[str; sprintf(['a' char(8869) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a_norm,resultsPS.a_norm_err)];
                            str=[str; sprintf(['a' char(8320) '*' char(963) '*(1+' char(957) ')/E = %.6f ' char(177) ' %.6f ' char(8491)],resultsPS.sig_ef,resultsPS.sig_ef_err)];

                            fprintf(['a' char(8741) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a_par,resultsPS.a_par_err);
                            fprintf(['a' char(8869) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a_norm,resultsPS.a_norm_err);
                            fprintf(['a' char(8320) '*' char(963) '*(1+' char(957) ')/E = %.6f ' char(177) ' %.6f ' char(8491) '\n'],resultsPS.sig_ef,resultsPS.sig_ef_err);
                        end    
                    end
                    str=[str; ' '];
                    str=[str; 'Aberations:'];
                    fprintf('Aberations:\n');
                    if GlobOpt.FitVerticalShift==1
                        str=[str; 'Sample vertical shift:'];
                        fprintf('Sample vertical shift:\n');
                        if handles.GroupsModels(g).Material.av==1
                            str=[str; sprintf(['sycos = %.4f ' char(177) ' %.4f'],resultsPS.sycos,resultsPS.sycos_err)];
                            fprintf(['sycos = %.4f ' char(177) ' %.4f\n'],resultsPS.sycos,resultsPS.sycos_err);
                        else
                            str=[str; sprintf(['sycos = %.4f ' char(177) ' %.4f'],resultsPS.sycos_ef,resultsPS.sycos_ef_err)];
                            fprintf(['sycos = %.4f ' char(177) ' %.4f\n'],resultsPS.sycos_ef,resultsPS.sycos_ef_err);
                        end
                    else
                        str=[str; 'Sample vertical shift: not fitted'];
                        fprintf('Sample vertical shift: not fitted\n');
                    end
                    if GlobOpt.FitDetectorShift==1
                        str=[str; 'Detector zero shift:'];
                        fprintf('Detector zero shift:\n');
                        if handles.GroupsModels(g).Material.av==1
                            str=[str; sprintf(['zero = %.4f ' char(177) ' %.4f deg'],resultsPS.dth*180/pi,resultsPS.dth_err*180/pi)];
                            fprintf(['zero = %.4f ' char(177) ' %.4f deg\n'],resultsPS.dth*180/pi,resultsPS.dth_err*180/pi);
                        else
                            str=[str; sprintf(['zero = %.4f ' char(177) ' %.4f'],resultsPS.dth_ef,resultsPS.dth_ef_err)];
                            fprintf(['zero = %.4f ' char(177) ' %.4f\n'],resultsPS.dth_ef,resultsPS.dth_ef_err);
                        end
                    else
                        str=[str; 'Detector zero shift: not fitted'];
                        fprintf('Detector zero shift: not fitted\n');
                    end
                case 'hexagonal'
                    if (strcmp(handles.GroupsModels(g).PeakShiftModel,'VookWitt') | strcmp(handles.GroupsModels(g).PeakShiftModel,'Reuss')) && handles.GroupsModels(g).Material.av==1 && isnan(handles.GroupsModels(g).Material.E_100)
                        str=[str; [char(957) char(8321) char(8320) char(8320) ' and E' char(8321) char(8320) char(8320) ' are not available, the tabled averages were used instead.']];
                    end
                    if handles.GroupsModels(g).Material.av==1 && ~strcmp(handles.GroupsModels(g).PeakShiftModel,'NoStress')
                        str=[str; sprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
                        str=[str; sprintf(['c' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.c0,resultsPS.c0_err)];
                        str=[str; sprintf([char(963) ' = %.4f ' char(177) ' %.4f GPa'],resultsPS.sig,resultsPS.sig_err)];
                        fprintf(['a' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);
                        fprintf(['c' char(8320) ' = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.c0,resultsPS.c0_err);
                        fprintf([char(963) ' = %.4f ' char(177) ' %.4f GPa\n'],resultsPS.sig,resultsPS.sig_err);
                    else
                        if strcmp(handles.GroupsModels(g).PeakShiftModel,'NoStress')
                            str=[str; sprintf(['a = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
                            str=[str; sprintf(['c = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.c0,resultsPS.c0_err)];
                        
                            fprintf(['a = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);
                            fprintf(['c = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.c0,resultsPS.c0_err);
                        else
                            str=[str; sprintf(['a = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.a0,resultsPS.a0_err)];
                            str=[str; sprintf(['c = %.4f ' char(177) ' %.4f ' char(8491)],resultsPS.c0,resultsPS.c0_err)];
                            str=[str; sprintf([char(963) '_ef = %.6f ' char(177) ' %.6f '],resultsPS.sig_ef,resultsPS.sig_ef_err)];
                        
                            fprintf(['a = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.a0,resultsPS.a0_err);
                            fprintf(['c = %.4f ' char(177) ' %.4f ' char(8491) '\n'],resultsPS.c0,resultsPS.c0_err);
                            fprintf([char(963) '_ef = %.6f ' char(177) ' %.6f \n'],resultsPS.sig_ef,resultsPS.sig_ef_err);
                        end    
                    end
            end
        end
        else
            str=[str; 'Not available.'];
            fprintf('Not available.\n');
        end
            
        str=[str; ' '];
        str=[str; 'Peak broadening from FWHM:'];
        if handles.GroupsModels(g).SizeOptions.Fitted==1
            str=[str; 'Grain shape model: ' handles.GroupsModels(g).SizeOptions.GrainsShape];
        else
            str=[str; 'Grain shape model: none'];
        end
        if ~strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Sphere')
            if strcmp(handles.GroupsModels(g).SizeOptions.GrownDirectionOpt,'Sample normal')
                str=[str; 'Growth direction: Sample normal'];
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    str=[str; 'Direction perpenducular to growth: Direction in sample surface with the azimuth ' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionPerpAzimuth) ' deg'];
                end
            else
                str=[str; 'Growth direction: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL) ']'];
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    str=[str; 'Direction perpenducular to growth: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL2) ']'];
                end
            end
        end
        if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            str=[str; 'Microstrain model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel];
        else
            str=[str; 'Microstrain model: none'];
        end
        
        fprintf(['Peak broadening from FWHM:\n']);
        if handles.GroupsModels(g).SizeOptions.Fitted==1
            fprintf(['Grain shape model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel '\n']);
        else
            fprintf(['Grain shape model: none\n']);
        end
        if ~strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Sphere')
            if strcmp(handles.GroupsModels(g).SizeOptions.GrownDirectionOpt,'Sample normal')
                fprintf('Growth direction: Sample normal\n');
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    fprintf(['Direction perpenducular to growth: Direction in sample surface with the azimuth ' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionPerpAzimuth) ' deg\n']);
                end
            else
                fprintf(['Growth direction: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL) ']\n']);
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    fprintf(['Direction perpenducular to growth: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL2) ']\n']);
                end
            end
        end
        if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            fprintf(['Microstrain model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel '\n']);
        else
            fprintf(['Microstrain model: none\n']);
        end

        if handles.PeakBroadeningResults_fwhm(g).av==1
            resultsPB_fwhm=handles.PeakBroadeningResults_fwhm(g).results;
            if handles.GroupsModels(g).SizeOptions.Fitted==1
                switch handles.GroupsModels(g).SizeOptions.GrainsShape
                    case 'Sphere'
                        str=[str; sprintf(['D = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D/10,resultsPB_fwhm.D_err/10)];
                        fprintf(['D = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D/10,resultsPB_fwhm.D_err/10);
                    case 'Cylinder'
                        str=[str; sprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D_norm/10,resultsPB_fwhm.D_norm_err/10)];
                        str=[str; sprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D_par/10,resultsPB_fwhm.D_par_err/10)];
                        fprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D_norm/10,resultsPB_fwhm.D_norm_err/10);
                        fprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D_par/10,resultsPB_fwhm.D_par_err/10);
                    case 'Cylinder with elliptical base'
                        str=[str; sprintf(['D1 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D1/10,resultsPB_fwhm.D1_err/10)];
                        str=[str; sprintf(['D2 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D2/10,resultsPB_fwhm.D2_err/10)];
                        str=[str; sprintf(['D3 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D3/10,resultsPB_fwhm.D3_err/10)];
                        fprintf(['D1 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D1/10,resultsPB_fwhm.D1_err/10);
                        fprintf(['D2 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D2/10,resultsPB_fwhm.D2_err/10);
                        fprintf(['D3 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D3/10,resultsPB_fwhm.D3_err/10);
                    case 'Rotational ellipsoid'
                        str=[str; sprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D_norm/10,resultsPB_fwhm.D_norm_err/10)];
                        str=[str; sprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D_par/10,resultsPB_fwhm.D_par_err/10)];
                        fprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D_norm/10,resultsPB_fwhm.D_norm_err/10);
                        fprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D_par/10,resultsPB_fwhm.D_par_err/10);
                    case 'General ellipsoid'
                        str=[str; sprintf(['D1 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D1/10,resultsPB_fwhm.D1_err/10)];
                        str=[str; sprintf(['D2 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D2/10,resultsPB_fwhm.D2_err/10)];
                        str=[str; sprintf(['D3 = %.4f ' char(177) ' %.4f nm'],resultsPB_fwhm.D3/10,resultsPB_fwhm.D3_err/10)];
                        fprintf(['D1 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D1/10,resultsPB_fwhm.D1_err/10);
                        fprintf(['D2 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D2/10,resultsPB_fwhm.D2_err/10);
                        fprintf(['D3 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_fwhm.D3/10,resultsPB_fwhm.D3_err/10);
                end
            end
            if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            switch handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel
                case 'Isotropic'
                    str=[str; sprintf(['e = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.e,resultsPB_fwhm.e_err)];
                    fprintf(['e = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.e,resultsPB_fwhm.e_err);
                case 'Dislocations (polycrystall)'
                    switch handles.GroupsModels(g).sys
                        case {'cubic_fcc','cubic_bcc'}
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err)];
                            str=[str; sprintf(['q = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.q_dis,resultsPB_fwhm.q_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err);
                            fprintf(['q = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.q_dis,resultsPB_fwhm.q_dis_err);
                        case 'hexagonal'
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err)];
                            str=[str; sprintf(['al = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.al_dis,resultsPB_fwhm.al_dis_err)];
                            str=[str; sprintf(['be = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.be_dis,resultsPB_fwhm.be_dis_err)];
                            str=[str; sprintf(['ga = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.ga_dis,resultsPB_fwhm.ga_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err);
                            fprintf(['al = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.al_dis,resultsPB_fwhm.al_dis_err);
                            fprintf(['be = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.be_dis,resultsPB_fwhm.be_dis_err);
                            fprintf(['ga = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.ga_dis,resultsPB_fwhm.ga_dis_err);
                    end
                case 'Dislocations (specific Burgers vector)'
                    str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err)];
                    str=[str; sprintf(['PsiB = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.PsiB,resultsPB_fwhm.PsiB_err)];
                    str=[str; sprintf(['PhiB = %.3f ' char(177) ' %.3f'],resultsPB_fwhm.PhiB,resultsPB_fwhm.PhiB_err)];
                    fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.e_dis_100,resultsPB_fwhm.e_dis_100_err);
                    fprintf(['PsiB = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.PsiB,resultsPB_fwhm.PsiB_err);
                    fprintf(['PsiB = %.3f ' char(177) ' %.3f\n'],resultsPB_fwhm.PhiB,resultsPB_fwhm.PhiB_err);
            end
            end
            if handles.GroupsModels(g).FitSF==1
                switch handles.GroupsModels(g).sys
                    case {'cubic_fcc','cubic_bcc'}
                        % the same for fcc and bcc
                        if isfield(resultsPS,'a0')
                            a=resultsPS.a0;
                        else
                            a=(resultsPS.a_par-resultsPS.a_norm)*0.4+resultsPS.a_norm; % sin2psi0=2/3/(1+1/3)=0.4
                        end
                        str=[str; sprintf(['al_ef = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.al_ef*a,resultsPB_fwhm.al_ef_err*a)];
                        fprintf(['al_ef = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.al_ef*a,resultsPB_fwhm.al_ef_err*a);
                    case 'hexagonal'
                        str=[str; sprintf(['al_ef = %.5f ' char(177) ' %.5f'],resultsPB_fwhm.al_ef,resultsPB_fwhm.al_ef_err)];
                        fprintf(['al_ef = %.5f ' char(177) ' %.5f\n'],resultsPB_fwhm.al_ef,resultsPB_fwhm.al_ef_err);
                end
            end
        else
            str=[str; 'Not available.'];
            fprintf('Not available.\n');
        end
        str=[str; ' '];
        str=[str; 'Peak broadening from ' char(946) ':'];
        if handles.GroupsModels(g).SizeOptions.Fitted==1
            str=[str; 'Grain shape model: ' handles.GroupsModels(g).SizeOptions.GrainsShape];
        else
            str=[str; 'Grain shape model: none'];
        end
        if ~strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Sphere')
            if strcmp(handles.GroupsModels(g).SizeOptions.GrownDirectionOpt,'Sample normal')
                str=[str; 'Growth direction: Sample normal'];
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    str=[str; 'Direction perpenducular to growth: Direction in sample surface with the azimuth ' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionPerpAzimuth) ' deg'];
                end
            else
                str=[str; 'Growth direction: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL) ']'];
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    str=[str; 'Direction perpenducular to growth: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL2) ']'];
                end
            end
        end
        if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            str=[str; 'Microstrain model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel];
        else
            str=[str; 'Microstrain model: none'];
        end
        
        fprintf(['Peak broadening from ' char(946) ':\n']);
        if handles.GroupsModels(g).SizeOptions.Fitted==1
            fprintf(['Grain shape model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel '\n']);
        else
            fprintf(['Grain shape model: none\n']);
        end
        if ~strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Sphere')
            if strcmp(handles.GroupsModels(g).SizeOptions.GrownDirectionOpt,'Sample normal')
                fprintf('Growth direction: Sample normal\n');
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    fprintf(['Direction perpenducular to growth: Direction in sample surface with the azimuth ' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionPerpAzimuth) ' deg\n']);
                end
            else
                fprintf(['Growth direction: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL) ']\n']);
                if strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'Cylinder with elliptical base') | strcmp(handles.GroupsModels(g).SizeOptions.GrainsShape,'General ellipsoid')
                    fprintf(['Direction perpenducular to growth: Crystallographic direction [' num2str(handles.GroupsModels(g).SizeOptions.GrownDirectionHKL2) ']\n']);
                end
            end
        end
        if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            fprintf(['Microstrain model: ' handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel '\n']);
        else
            fprintf(['Microstrain model: none\n']);
        end
        if handles.PeakBroadeningResults_beta(g).av==1
            resultsPB_beta=handles.PeakBroadeningResults_beta(g).results;
            if handles.GroupsModels(g).SizeOptions.Fitted==1
                switch handles.GroupsModels(g).SizeOptions.GrainsShape
                    case 'Sphere'
                        str=[str; sprintf(['D = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D/10,resultsPB_beta.D_err/10)];
                        fprintf(['D = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D/10,resultsPB_beta.D_err/10);
                    case 'Cylinder'
                        str=[str; sprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D_norm/10,resultsPB_beta.D_norm_err/10)];
                        str=[str; sprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D_par/10,resultsPB_beta.D_par_err/10)];
                        fprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D_norm/10,resultsPB_beta.D_norm_err/10);
                        fprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D_par/10,resultsPB_beta.D_par_err/10);
                    case 'Cylinder with elliptical base'
                        str=[str; sprintf(['D1 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D1/10,resultsPB_beta.D1_err/10)];
                        str=[str; sprintf(['D2 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D2/10,resultsPB_beta.D2_err/10)];
                        str=[str; sprintf(['D3 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D3/10,resultsPB_beta.D3_err/10)];
                        fprintf(['D1 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D1/10,resultsPB_beta.D1_err/10);
                        fprintf(['D2 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D2/10,resultsPB_beta.D2_err/10);
                        fprintf(['D3 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D3/10,resultsPB_beta.D3_err/10);
                    case 'Rotational ellipsoid'
                        str=[str; sprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D_norm/10,resultsPB_beta.D_norm_err/10)];
                        str=[str; sprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D_par/10,resultsPB_beta.D_par_err/10)];
                        fprintf(['D' char(8869) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D_norm/10,resultsPB_beta.D_norm_err/10);
                        fprintf(['D' char(8741) ' = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D_par/10,resultsPB_beta.D_par_err/10);
                    case 'General ellipsoid'
                        str=[str; sprintf(['D1 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D1/10,resultsPB_beta.D1_err/10)];
                        str=[str; sprintf(['D2 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D2/10,resultsPB_beta.D2_err/10)];
                        str=[str; sprintf(['D3 = %.4f ' char(177) ' %.4f nm'],resultsPB_beta.D3/10,resultsPB_beta.D3_err/10)];
                        fprintf(['D1 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D1/10,resultsPB_beta.D1_err/10);
                        fprintf(['D2 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D2/10,resultsPB_beta.D2_err/10);
                        fprintf(['D3 = %.4f ' char(177) ' %.4f nm\n'],resultsPB_beta.D3/10,resultsPB_beta.D3_err/10);
                end
            end
            if handles.GroupsModels(g).MicrostrainOptions.Fitted==1
            switch handles.GroupsModels(g).MicrostrainOptions.MicrostrainModel
                case 'Isotropic'
                    str=[str; sprintf(['e = %.5f ' char(177) ' %.5f'],resultsPB_beta.e,resultsPB_beta.e_err)];
                    fprintf(['e = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e,resultsPB_beta.e_err);
                case 'Dislocations (polycrystall)'
                    switch handles.GroupsModels(g).sys
                        case {'cubic_fcc';'cubic_bcc'}
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err)];
                            str=[str; sprintf(['q = %.3f ' char(177) ' %.3f'],resultsPB_beta.q_dis,resultsPB_beta.q_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err);
                            fprintf(['q = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.q_dis,resultsPB_beta.q_dis_err);
                        case 'hexagonal'
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err)];
                            str=[str; sprintf(['al = %.3f ' char(177) ' %.3f'],resultsPB_beta.al_dis,resultsPB_beta.al_dis_err)];
                            str=[str; sprintf(['be = %.3f ' char(177) ' %.3f'],resultsPB_beta.be_dis,resultsPB_beta.be_dis_err)];
                            str=[str; sprintf(['ga = %.3f ' char(177) ' %.3f'],resultsPB_beta.ga_dis,resultsPB_beta.ga_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err);
                            fprintf(['al = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.al_dis,resultsPB_beta.al_dis_err);
                            fprintf(['be = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.be_dis,resultsPB_beta.be_dis_err);
                            fprintf(['ga = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.ga_dis,resultsPB_beta.ga_dis_err);
                    end
                    switch handles.GroupsModels(g).sys
                        case {'cubic_fcc','cubic_bcc'}
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err)];
                            str=[str; sprintf(['q = %.3f ' char(177) ' %.3f'],resultsPB_beta.q_dis,resultsPB_beta.q_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err);
                            fprintf(['q = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.q_dis,resultsPB_beta.q_dis_err);
                        case 'hexagonal'
                            str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err)];
                            str=[str; sprintf(['al = %.3f ' char(177) ' %.3f'],resultsPB_beta.al_dis,resultsPB_beta.al_dis_err)];
                            str=[str; sprintf(['be = %.3f ' char(177) ' %.3f'],resultsPB_beta.be_dis,resultsPB_beta.be_dis_err)];
                            str=[str; sprintf(['ga = %.3f ' char(177) ' %.3f'],resultsPB_beta.ga_dis,resultsPB_beta.ga_dis_err)];
                            fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err);
                            fprintf(['al = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.al_dis,resultsPB_beta.al_dis_err);
                            fprintf(['be = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.be_dis,resultsPB_beta.be_dis_err);
                            fprintf(['ga = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.ga_dis,resultsPB_beta.ga_dis_err);
                    end
                case 'Dislocations (specific Burgers vector)'
                    str=[str; sprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err)];
                    str=[str; sprintf(['PsiB = %.3f ' char(177) ' %.3f'],resultsPB_beta.PsiB,resultsPB_beta.PsiB_err)];
                    str=[str; sprintf(['PhiB = %.3f ' char(177) ' %.3f'],resultsPB_beta.PhiB,resultsPB_beta.PhiB_err)];
                    fprintf(['e' char(8321) char(8320) char(8320) ' = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.e_dis_100,resultsPB_beta.e_dis_100_err);
                    fprintf(['PsiB = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.PsiB,resultsPB_beta.PsiB_err);
                    fprintf(['PsiB = %.3f ' char(177) ' %.3f\n'],resultsPB_beta.PhiB,resultsPB_beta.PhiB_err);
            
            end
            end
            if handles.GroupsModels(g).FitSF==1
                % the same for fcc and bcc
                switch handles.GroupsModels(g).sys
                    case {'cubic_fcc','cubic_bcc'}
                        if isfield(resultsPS,'a0')
                            a=resultsPS.a0;
                        else
                            a=(resultsPS.a_par-resultsPS.a_norm)*0.4+resultsPS.a_norm; % sin2psi0=2/3/(1+1/3)=0.4
                        end
                        str=[str; sprintf(['al_ef = %.5f ' char(177) ' %.5f'],resultsPB_beta.al_ef*a,resultsPB_beta.al_ef_err*a)];
                        fprintf(['al_ef = %.5f ' char(177) ' %.5fu\n'],resultsPB_beta.al_ef*a,resultsPB_beta.al_ef_err*a);
                    case 'hexagonal'
                        str=[str; sprintf(['al_ef = %.5f ' char(177) ' %.5f'],resultsPB_beta.al_ef,resultsPB_beta.al_ef_err)];
                        fprintf(['al_ef = %.5f ' char(177) ' %.5f\n'],resultsPB_beta.al_ef,resultsPB_beta.al_ef_err);
                end
            end
        else
            str=[str; 'Not available.'];
            fprintf('Not available.\n');
        end
        if GlobOpt.FitGroups==1
            str=[str; '----------------------------------'];
            fprintf('----------------------------------\n');
        end
    end
end


function Wavelength_ed_Callback(hObject, eventdata, handles)
% hObject    handle to Wavelength_ed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wavelength_ed as text
%        str2double(get(hObject,'String')) returns contents of Wavelength_ed as a double
iserr=0;
try
    hlp=str2num(handles.Wavelength_ed.String);
    if isempty(hlp) || hlp<0
        iserr=1;
    end
catch
    iserr=1;
end
if iserr==1
    errordlg('Wavelength should be a number greater than zero!','Wrong input')
end

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


% --- Executes on button press in CijMat_rad.
function CijMat_rad_Callback(hObject, eventdata, handles)
% hObject    handle to CijMat_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CijMat_rad
handles.MaterialConstants_tab.Enable='off';
handles.Stiffness_tab.Enable='on';
handles.Stiffness1_txt.Enable='on';
handles.Stiffness1_txt.String='C=(';
handles.Stiffness2_txt.Enable='on';
handles.MaterialDatabase_list.Enable='on';
handles.MaterialDatabase_txt.String='Material database (Sij/Cij)';
handles.MaterialDatabase_list.String=handles.StiffnessMatrices_database.Label;
handles.MaterialDatabase_list.Value=1;
switch handles.MaterialDatabase_list.String{1}
    case 'cubic_fcc'
        handles.System_pop.Value=1;
    case 'cubic_bcc'
        handles.System_pop.Value=2;
    case 'hexagonal'
        handles.System_pop.Value=3;
end
MaterialDatabase_list_Callback(handles.MaterialDatabase_list, eventdata, handles);

I=handles.GroupSelected;
handles.GroupsModels(I).Material.DatabaseRow=1;
handles.GroupsModels(I).Material.av=1;
handles.GroupsModels(I).Material.Option='Cij';
tab=handles.MaterialConstants_tab.Data;
handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.Cij=cell2mat(handles.Stiffness_tab.Data);
handles.GroupsModels(I).Material.Sij=inv(handles.GroupsModels(I).Material.Cij);


% --- Executes when entered data in editable cell(s) in Stiffness_tab.
function Stiffness_tab_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Stiffness_tab (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
val=hObject.Data{ind(1),ind(2)};
switch handles.System_pop.String{handles.System_pop.Value}
    case {'cubic_fcc','cubic_bcc'}
        if ind(1)==ind(2)
            if ind(1)<=3
                hObject.Data{1,1}=val;
                hObject.Data{2,2}=val;
                hObject.Data{3,3}=val;
            else
                hObject.Data{4,4}=val;
                hObject.Data{5,5}=val;
                hObject.Data{6,6}=val;
            end
        else
            if ind(1)<=3 && ind(2)<=3
                hObject.Data{1,2}=val;
                hObject.Data{1,3}=val;
                hObject.Data{2,1}=val;
                hObject.Data{2,3}=val;
                hObject.Data{3,1}=val;
                hObject.Data{3,2}=val;
            else
                hObject.Data{ind(1),ind(2)}=0;
            end
        end
    case 'hexagonal'
        if ind(1)==ind(2)
            if ind(1)<=2
                hObject.Data{1,1}=val;
                hObject.Data{2,2}=val;
            elseif ind(1)==3
                hObject.Data{3,3}=val;
            elseif ind(1)==4 || ind(1)==5
                hObject.Data{4,4}=val;
                hObject.Data{5,5}=val;
            else
                hObject.Data{6,6}=val;
                if handles.Sij_rad.Value==1
                    % s66=2*(s11-s12) -> s12=s11-0.5*s66
                    s12=hObject.Data{1,1}-0.5*val;
                    hObject.Data{1,2}=s12;
                    hObject.Data{2,1}=s12;
                else % Cij
                    % c66=0.5*(c11-c12) -> c12=c11-2*c66
                    c12=hObject.Data{1,1}-2*val;
                    hObject.Data{1,2}=c12;
                    hObject.Data{2,1}=c12;
                end
            end
        else
            if ind(1)+ind(2)==3
                hObject.Data{1,2}=val;
                hObject.Data{2,1}=val;
                if handles.Sij_rad.Value==1
                    % s66=2*(s11-s12)
                    s66=2*(hObject.Data{1,1}-val);
                    hObject.Data{6,6}=s66;
                else % Cij
                    % c66=0.5*(c11-c12)
                    c66=0.5*(hObject.Data{1,1}-val);
                    hObject.Data{6,6}=c66;
                end
            elseif (ind(1)==3 && (ind(2)==1 || ind(2)==2)) || (ind(2)==3 && (ind(1)==1 || ind(1)==2))
                hObject.Data{1,3}=val;
                hObject.Data{2,3}=val;
                hObject.Data{3,1}=val;
                hObject.Data{3,2}=val;
            else
                hObject.Data{ind(1),ind(2)}=0;
            end
        end
end
guidata(hObject,handles)
if handles.SijMat_rad.Value==1
    S=cell2mat(hObject.Data);
    C=inv(S);
    nu100=-hObject.Data{2,1}/hObject.Data{1,1};
    E100=1/hObject.Data{1,1};
    G100=1/hObject.Data{4,4};
elseif handles.CijMat_rad.Value==1
    C=cell2mat(hObject.Data);
    S=inv(C);

    nu100=-S(2,1)/(2*S(1,1));
    E100=1/S(1,1);
    G100=1/S(4,4);
end
[nu,E,G,~]=ElasticConstants_VoigtReussHill(C,handles.System_pop.String{handles.System_pop.Value});
handles.MaterialConstants_tab.Data={nu E G nu100 E100 G100};

I=handles.GroupSelected;
%     handles.GroupsModels(I).Material.av=1;
switch handles.GroupsModels(I).Material.Option
    case 'Sij'
        tab=handles.MaterialConstants_tab.Data;
        handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
        handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
        handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
        handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
        handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
        handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)

        handles.GroupsModels(I).Material.Sij=cell2mat(handles.Stiffness_tab.Data);
        handles.GroupsModels(I).Material.Cij=inv(handles.GroupsModels(I).Material.Sij);
    case 'Cij'
        tab=handles.MaterialConstants_tab.Data;
        handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
        handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
        handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
        handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
        handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
        handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)

        handles.GroupsModels(I).Material.Cij=cell2mat(handles.Stiffness_tab.Data);
        handles.GroupsModels(I).Material.Sij=inv(handles.GroupsModels(I).Material.Cij);
end
guidata(hObject,handles)


% --- Executes when entered data in editable cell(s) in DataTable_tab.
function DataTable_tab_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to DataTable_tab (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
G=handles.Groups_list.Value; % group number
I=eventdata.Indices; % peak number
if I(2)==10 % for safety, but other option should not occur, since other columns are not editable (instead of label)
    handles.DataGroup(I(1),G)=1*handles.DataTable_tab.Data{I(1),I(2)};
end
guidata(hObject,handles)


% --- Executes on button press in SaveFigures_btn.
function SaveFigures_btn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigures_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastFilePath')
    [FileName,PathName,FilterIndex] = uiputfile({'*.png','PNG image (*.png)';'*.tiff','TIFF image (*.tiff)';'*.pdf','PDF file (*.pdf)';'*.fig','Matlab figure (*.fig)'},'Save figures',handles.LastFilePath);
    handles.LastFilePath=PathName;
else
    [FileName,PathName,FilterIndex] = uiputfile({'*.png','PNG image (*.png)';'*.tiff','TIFF image (*.tiff)';'*.pdf','PDF file (*.pdf)';'*.fig','Matlab figure (*.fig)'},'Save figures');
    handles.LastFilePath=PathName;
end
if ~isempty(FileName)
    [~,filename,ext]=fileparts(FileName);
    if strcmp(ext,'.fig')
        fignew = figure('Visible','off'); % Invisible figure
        newAxes = copyobj(handles.PeakShift_ax,fignew); % Copy the appropriate axes
        set(newAxes,'Units','normalized','Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
        set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
        savefig(fignew,[PathName '\' filename '_PeakShift' ext]);
        delete(fignew);

        fignew = figure('Visible','off'); % Invisible figure
        newAxes = copyobj(handles.PeakBroadening_ax,fignew); % Copy the appropriate axes
        set(newAxes,'Units','normalized','Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
        set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
        savefig(fignew,[PathName '\' filename '_PeakBroadening' ext]);
        delete(fignew);
    else  
        exportgraphics(handles.PeakShift_ax,[PathName '\' filename '_PeakShift' ext],'Resolution',300);
        exportgraphics(handles.PeakBroadening_ax,[PathName '\' filename '_PeakBroadening' ext],'Resolution',300);
    end
end

% --- Executes on button press in SaveResults_btn.
function SaveResults_btn_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResults_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastFilePath')
    [FileName,PathName,FilterIndex] = uiputfile({'*.txt','Text file (*.txt)'},'Save figures',handles.LastFilePath);
    handles.LastFilePath=PathName;
else
    [FileName,PathName,FilterIndex] = uiputfile({'*.txt','Text file (*.txt)'},'Save figures');
    handles.LastFilePath=PathName;
end
if ~isempty(FileName)
    mylines={['Processed file: ' handles.Filename];
        ['Date and time: ' datestr(datetime)];
        ['Peak shift model: ' handles.GlobOpt.PeakShiftModel];
        ['Peak broadening model: ' handles.GlobOpt.MicrostrainModel];
        ['Crystalographic system: ' handles.GlobOpt.sys];
        ['Fit of vertical shift: ' num2str(handles.GlobOpt.FitVerticalShift)];
        ['Fit of vertical shift: ' num2str(handles.GlobOpt.FitDetectorShift)];
        ['Fit coherently diffracted domains size: ' num2str(handles.GlobOpt.FitSize)];
        ['Fit of stacking faults: ' num2str(handles.GlobOpt.FitSF)];
        '---------------------------------';
        ' ';
        'DATA:';
        '---------------------------------';
        sprintf('h\tk\tl\ttth (deg)\tpsi (deg)\tFWHM (deg)\tbeta (deg)\tlabel');
        '---------------------------------'};
    for r=1:handles.DataPointsNumber
        str=sprintf('%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s',handles.Data.h(r),handles.Data.k(r),handles.Data.l(r),handles.Data.tth(r),handles.Data.psi(r),handles.Data.fwhm(r),handles.Data.beta(r),handles.Data.label{r});
        mylines=[mylines; str];
    end
    mylines=[mylines;
        '---------------------------------';
        ' ';
        'RESULTS:';
        handles.Status_list.String];
    writelines(mylines,[PathName '\' FileName]);
end

% --- Executes on button press in SubstractInstrumental_chbx.
function SubstractInstrumental_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to SubstractInstrumental_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SubstractInstrumental_chbx


% --- Executes on button press in InstrumentalFunction_btn.
function InstrumentalFunction_btn_Callback(hObject, eventdata, handles)
% hObject    handle to InstrumentalFunction_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'IP_Parameters')
    handles.IP_Parameters=InstrumentalProfile(handles.IP_Parameters);
else
    handles.IP_Parameters=InstrumentalProfile;
end
handles.SubstractInstrumental_chbx.Enable='on';
guidata(hObject, handles);

function [FWHM,beta]=ComputeInstrumentalBroadening(IP,tth,chi)
% U(chi)
switch IP.UonChi
    case '1'
        U=IP.u0+IP.u1*1;
    case 'chi'
        U=IP.u0+IP.u1.*chi;
    case 'sin(chi)'
        U=IP.u0+IP.u1.*sind(chi);
    case 'cos(chi)'
        U=IP.u0+IP.u1.*cosd(chi);
    case 'tan(chi)'
        U=IP.u0+IP.u1.*tand(chi);
    case 'cotan(chi)'
        U=IP.u0+IP.u1./tand(chi);
end

% V(chi)
switch IP.VonChi
    case '1'
        V=IP.v0+IP.v1*1;
    case 'chi'
        V=IP.v0+IP.v1.*chi;
    case 'sin(chi)'
        V=IP.v0+IP.v1.*sind(chi);
    case 'cos(chi)'
        V=IP.v0+IP.v1.*cosd(chi);
    case 'tan(chi)'
        V=IP.v0+IP.v1.*tand(chi);
    case 'cotan(chi)'
        V=IP.v0+IP.v1./tand(chi);
end

% W
switch IP.WonChi
    case '1'
        W=IP.w0+IP.w1*1;
    case 'chi'
        W=IP.w0+IP.w1.*chi;
    case 'sin(chi)'
        W=IP.w0+IP.w1.*sind(chi);
    case 'cos(chi)'
        W=IP.w0+IP.w1.*cosd(chi);
    case 'tan(chi)'
        W=IP.w0+IP.w1.*tand(chi);
    case 'cotan(chi)'
        W=IP.w0+IP.w1./tand(chi);
end

switch IP.Eta0onChi
    case '1'
        eta0=IP.eta00+IP.eta01*1;
    case 'chi'
        eta0=IP.eta00+IP.eta01.*chi;
    case 'sin(chi)'
        eta0=IP.eta00+IP.eta01.*sind(chi);
    case 'cos(chi)'
        eta0=IP.eta00+IP.eta01.*cosd(chi);
    case 'tan(chi)'
        eta0=IP.eta00+IP.eta01.*tand(chi);
    case 'cotan(chi)'
        eta0=IP.eta00+IP.eta01./tand(chi);
end

switch IP.Eta1onChi
    case '1'
        eta1=IP.eta10+IP.eta11*1;
    case 'chi'
        eta1=IP.eta10+IP.eta11.*chi;
    case 'sin(chi)'
        eta1=IP.eta10+IP.eta11.*sind(chi);
    case 'cos(chi)'
        eta1=IP.eta10+IP.eta11.*cosd(chi);
    case 'tan(chi)'
        eta1=IP.eta10+IP.eta11.*tand(chi);
    case 'cotan(chi)'
        eta1=IP.eta10+IP.eta11./tand(chi);
end

% Asym parameter has no usage in this version

FWHM=sqrt(U.*(tand(tth/2)).^2+V.*tand(tth/2)+W); % Cagliotti formula
eta=eta0+eta1.*tth*pi/180; % tth should be in radians
beta_L=pi/2.*FWHM;
beta_G=0.5*sqrt(pi/log(2)).*FWHM;
beta=eta.*beta_L+(1-beta_G).*beta_G;


% --- Executes on button press in Microstrain_chbx.
function Microstrain_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to Microstrain_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Microstrain_chbx
I=handles.GroupSelected;
if hObject.Value==0
    handles.GroupsModels(I).MicrostrainOptions.Fitted=0;
else
    handles.GroupsModels(I).MicrostrainOptions.Fitted=1;
end
guidata(hObject, handles);


% --- Executes on button press in VerticalSampleShift_chbx.
function VerticalSampleShift_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to VerticalSampleShift_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VerticalSampleShift_chbx


% --- Executes on button press in DetectorShift_chbx.
function DetectorShift_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to DetectorShift_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DetectorShift_chbx


% --- Executes on button press in FitHcpLP_chbx.
function FitHcpLP_chbx_Callback(hObject, eventdata, handles)
% hObject    handle to FitHcpLP_chbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FitHcpLP_chbx


function N=EstimateStrainSizeEffect(tth,broad,lam)
broad_rs=1/lam*cosd(tth/2).*broad*pi/360;
% b=b_size+b_strain*sin(tth)
pT1=polyfit(sind(tth),broad_rs,1);
N1=tanh(log(pT1(1)/pT1(2)))/2+1.5; % empirical function going from value 1 to 2

% b.^2=b_size.^2+b_strain.^2*sin(tth).^2
pT2=polyfit(sind(tth).^2,broad_rs.^2,1);
N2=tanh(log(pT2(1)/pT2(2)))/2+1.5; % empirical function going from value 1 to 2

if pT1(2)<0 || pT2(2)<0
    % size is probably to big, which reduces the size term in the equation,
    % therefore the strain effect is dominant 
    N=2;
elseif pT1(1)<0 || pT2(1)<0
    % microstrain is probably to low, therefore the size effect is dominant 
    N=1;
elseif (imag(N1)+imag(N2)>0)
    % one the pT parameters is negative, both cannot be, because the
    % broadening would be negative
    N=NaN;
else
    N=(N1+N2)/2;
end


% --- Executes on button press in AdvancedShiftOpt_btn.
function AdvancedShiftOpt_btn_Callback(hObject, eventdata, handles)
% hObject    handle to AdvancedShiftOpt_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.GroupSelected;
InputData.MyOptions=handles.GroupsModels(I).AdvancedShiftOptions;
InputData.DataGroupN=size(handles.DataGroup,2)-1;
OutputData=AdvancedShiftOpt_GUI(InputData);
% handles.AdvancedShiftOptions=OutputData;
handles.GroupsModels(I).AdvancedShiftOptions=OutputData;
guidata(hObject, handles);


% --- Executes on button press in SizeOptions_btn.
function SizeOptions_btn_Callback(hObject, eventdata, handles)
% hObject    handle to SizeOptions_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.GroupSelected;
InputData=handles.GroupsModels(I).SizeOptions;
% InputData=handles.SizeOptions;
OutputData=SizeOptions_GUI(InputData);
% handles.SizeOptions=OutputData;
% disp('OutputData:')
% disp(OutputData)
handles.GroupsModels(I).SizeOptions.GrainsShape=OutputData.GrainsShape;
handles.GroupsModels(I).SizeOptions.ParGuess=OutputData.ParGuess;
handles.GroupsModels(I).SizeOptions.ParFix=OutputData.ParFix;
handles.GroupsModels(I).SizeOptions.ParLB=OutputData.ParLB;
handles.GroupsModels(I).SizeOptions.ParUB=OutputData.ParUB;
handles.GroupsModels(I).SizeOptions.GrownDirectionOpt=OutputData.GrownDirectionOpt;
handles.GroupsModels(I).SizeOptions.GrownDirectionHKL=OutputData.GrownDirectionHKL;
% handles.GroupsModels(I).SizeOptions.GrownDirectionOpt2=OutputData.GrownDirectionOpt2;
handles.GroupsModels(I).SizeOptions.GrownDirectionHKL2=OutputData.GrownDirectionHKL2;
handles.GroupsModels(I).SizeOptions.GrownDirectionPerpAzimuth=OutputData.GrownDirectionPerpAzimuth;
guidata(hObject, handles);
% disp('Actual size parameters:')
% numel(handles.GroupsModels)
% numel(handles.DataGroupMaxNumber)
% for g=1:numel(handles.GroupsModels)
%     disp(g)
%     disp(handles.GroupsModels(g).SizeOptions.GrainsShape);
%     disp(handles.GroupsModels(g).SizeOptions.ParGuess);
%     disp(handles.GroupsModels(g).SizeOptions.ParFix);
%     disp(handles.GroupsModels(g).SizeOptions.ParLB);
%     disp(handles.GroupsModels(g).SizeOptions.ParUB);
% end
% guidata(hObject, handles);


% --- Executes on selection change in GroupModelSelection_pop.
function GroupModelSelection_pop_Callback(hObject, eventdata, handles)
% hObject    handle to GroupModelSelection_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupModelSelection_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupModelSelection_pop


% --- Executes during object creation, after setting all properties.
function GroupModelSelection_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupModelSelection_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NoStress_rad.
function NoStress_rad_Callback(hObject, eventdata, handles)
% hObject    handle to NoStress_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NoStress_rad
I=handles.GroupSelected;
handles.GroupsModels(I).PeakShiftModel='NoStress';
guidata(hObject, handles);


% --- Executes on button press in IsotropicStress_rad.
function IsotropicStress_rad_Callback(hObject, eventdata, handles)
% hObject    handle to IsotropicStress_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IsotropicStress_rad
I=handles.GroupSelected;
handles.GroupsModels(I).PeakShiftModel='IsotropicStress';
guidata(hObject, handles);


% --- Executes on button press in VookWitt_rad.
function VookWitt_rad_Callback(hObject, eventdata, handles)
% hObject    handle to VookWitt_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VookWitt_rad
I=handles.GroupSelected;
handles.GroupsModels(I).PeakShiftModel='VookWitt';
guidata(hObject, handles);


% --- Executes on button press in Reuss_rad.
function Reuss_rad_Callback(hObject, eventdata, handles)
% hObject    handle to Reuss_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Reuss_rad
I=handles.GroupSelected;
handles.GroupsModels(I).PeakShiftModel='Reuss';
guidata(hObject, handles);


% --- Executes on button press in AdvancedShiftOpt_rad.
function AdvancedShiftOpt_rad_Callback(hObject, eventdata, handles)
% hObject    handle to AdvancedShiftOpt_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AdvancedShiftOpt_rad
I=handles.GroupSelected;
handles.GroupsModels(I).PeakShiftModel='AdvancedShiftOpt';
guidata(hObject, handles);


% % --- Executes on button press in IsotropicMicrostrain_rad.
% function IsotropicMicrostrain_rad_Callback(hObject, eventdata, handles)
% % hObject    handle to IsotropicMicrostrain_rad (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of IsotropicMicrostrain_rad
% I=handles.GroupSelected;
% handles.GroupsModels(I).MicrostrainOptions.MicrostrainModel='Isotropic';
% guidata(hObject, handles);


% % --- Executes on button press in Dislocations_rad.
% function Dislocations_rad_Callback(hObject, eventdata, handles)
% % hObject    handle to Dislocations_rad (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of Dislocations_rad
% I=handles.GroupSelected;
% handles.GroupsModels(I).MicrostrainOptions.MicrostrainModel='Dislocations';
% guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in MaterialConstants_tab.
function MaterialConstants_tab_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to MaterialConstants_tab (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
I=handles.GroupSelected;
tab=handles.MaterialConstants_tab.Data;
handles.GroupsModels(I).Material.nu=tab{1,1}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E=tab{1,2}; % Poisson ration
handles.GroupsModels(I).Material.G=tab{1,3}; % Shear modulus (GPa)
handles.GroupsModels(I).Material.nu_100=tab{1,4}; % Young modulus (GPa)
handles.GroupsModels(I).Material.E_100=tab{1,5}; % Poisson ration
handles.GroupsModels(I).Material.G_100=tab{1,6}; % Shear modulus (GPa)
guidata(hObject, handles);


function txt = myupdatefcn(~,event_obj,Data,DataGroup)
% Own text for the data cursor in the figure
% It is neccessary to distinguish between the substrate and layer peaks -
% this code is quite "dirty", but it works at least
I=get(event_obj,'DataIndex');
tar=get(event_obj,'Target');
ind1=strfind(tar.Tag,'_');
ind2=strfind(tar.Tag,'#');
DataOrFit=tar.Tag(1:ind1-1);
PSOrFWHMOrBeta=tar.Tag(ind1+1:ind2-2);
g=str2num(tar.Tag(ind2+1:end));

DataLabel=['Group #' num2str(g-1) ': ' DataOrFit];
% We need to take care of the fact, that the plotted data might not be
% sorted as in the order from the table
indG=find(DataGroup(:,g)==1);
switch PSOrFWHMOrBeta
    case 'PS'
        sin2psi=sind(Data.psi(indG)).^2;
        [~,indsort]=sort(sin2psi);
    otherwise
        sin2th=sind(Data.tth(indG)/2).^2;
        [~,indsort]=sort(sin2th);
end
ind=indG(indsort);
mylabel=['h k l: ',num2str(Data.h(ind(I))) ' ' num2str(Data.k(ind(I))) ' ' num2str(Data.l(ind(I)))];

txt = {DataLabel,...
       mylabel,...
       ['2Theta: ' num2str(round(Data.tth(ind(I)),4)) ' deg'],...
       ['psi: ', num2str(round(Data.psi(ind(I)),2)) ' deg'],...
       ['phi: ', num2str(round(Data.phi(ind(I)),2)) ' deg'],...
       ['FWHM: ', num2str(round(Data.fwhm(ind(I)),4)) ' deg'],...
       ['beta: ', num2str(round(Data.beta(ind(I)),4)) ' deg'],...
      };


function lam=ScanHeaderForWavelength(filename,nh)
% Function will try to scan file header, if there is an available
% information about wavelength
% lam - wavelength value, NaN if not found
hlp=readlines(filename);
MyHeader=hlp(1:nh);
KeywordFound=0;
WavelengthFound=0;
for r=1:numel(MyHeader)
    mystr=lower(MyHeader{r});
    if KeywordFound
       mystr(strfind(mystr, '=')) = [];
       KeyWord   = 'wavelength';
       Index = strfind(mystr,KeyWord);
       value = sscanf(mystr(Index(1) + length(KeyWord):end), '%g', 1);
       if ~isempty(value)
           WavelengthFound=1;
           break;
       end
    else
        if contains(mystr,'wavelength')
            KeywordFound=1;
            mystr(strfind(mystr, '=')) = [];
            KeyWord   = 'wavelength';
            Index = strfind(mystr,KeyWord);
            value = sscanf(mystr(Index(1) + length(KeyWord):end), '%g', 1);
            if ~isempty(value)
                WavelengthFound=1;
                break;
            end
        end
    end
end
if WavelengthFound
    lam=value;
else
    lam=NaN;
end


% --- Executes on button press in MicrostrainOptions_btn.
function MicrostrainOptions_btn_Callback(hObject, eventdata, handles)
% hObject    handle to MicrostrainOptions_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.GroupSelected;
InputData=handles.GroupsModels(I).MicrostrainOptions;
OutputData=MicrostrainOptions_GUI(InputData);
handles.GroupsModels(I).MicrostrainOptions.MicrostrainModel=OutputData.MicrostrainModel;
handles.GroupsModels(I).MicrostrainOptions.ParGuess=OutputData.ParGuess;
handles.GroupsModels(I).MicrostrainOptions.ParFix=OutputData.ParFix;
handles.GroupsModels(I).MicrostrainOptions.ParLB=OutputData.ParLB;
handles.GroupsModels(I).MicrostrainOptions.ParUB=OutputData.ParUB;
guidata(hObject, handles);


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('Documentation\Introduction.html','-new')

% --------------------------------------------------------------------
function OpenProject_menuit_Callback(hObject, eventdata, handles)
% hObject    handle to OpenProject_menuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('All unsaved changes will be lost. Do you want to save your project?','Save the project','Yes','No','Yes');
if strcmp(answer,'Yes')
    if isfield(handles,'LastFilePath')
        [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project',handles.LastFilePath);
        handles.LastFilePath=PathName;
    else
        [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project');
        handles.LastFilePath=PathName;
    end
    save([PathName FileName],'handles','-mat');
    handles.figure1.Name=['PSB_GUI: ' FileName];
    guidata(hObject, handles);
end
if strcmp(answer,'No') || FileName~=0
    if isfield(handles,'LastFilePath')
        [FileName,PathName,~] = uigetfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Load your previous project',handles.LastFilePath);
        handles.LastFilePath=PathName;
    else
        [FileName,PathName,~] = uigetfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Load your previous project');
        handles.LastFilePath=PathName;
    end
    if FileName~=0
        hf=gcf;
        load([PathName FileName],'-mat');
        handles.figure1.Name=['PSB_GUI: ' FileName];
        guidata(hObject, handles);
        pause(0.001); % wait little bit to close the current window after opening the new one

        delete(hf);
    end
end


% --------------------------------------------------------------------
function SaveProject_menuit_Callback(hObject, eventdata, handles)
% hObject    handle to SaveProject_menuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'LastFilePath')
    [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project',handles.LastFilePath);
    handles.LastFilePath=PathName;
else
    [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project');
    handles.LastFilePath=PathName;
end
if FileName~=0
    save([PathName FileName],'handles','-mat');
    handles.figure1.Name=['PSB_GUI: ' FileName];
    guidata(hObject, handles);
end


% --- Executes on button press in DataImportHelp_btn.
function DataImportHelp_btn_Callback(hObject, eventdata, handles)
% hObject    handle to DataImportHelp_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('Documentation\Data import.html','-new')


% --- Executes on button press in GroupSeparatorHelp_btn.
function GroupSeparatorHelp_btn_Callback(hObject, eventdata, handles)
% hObject    handle to GroupSeparatorHelp_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('Documentation\Groups separation.html','-new')


% --- Executes on button press in MaterialInformationHelp_btn.
function MaterialInformationHelp_btn_Callback(hObject, eventdata, handles)
% hObject    handle to MaterialInformationHelp_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('Documentation\Material information.html','-new')


% --- Executes on button press in FitAndResultsHelp_btn.
function FitAndResultsHelp_btn_Callback(hObject, eventdata, handles)
% hObject    handle to FitAndResultsHelp_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('Documentation\Fit and Results.html','-new')


% --- Executes on button press in QuickSelector_btn.
function QuickSelector_btn_Callback(hObject, eventdata, handles)
% hObject    handle to QuickSelector_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InputData.Labels=handles.Data.label;
G=handles.Groups_list.Value; % group number
InputData.CurrentSelection=handles.DataGroup(:,G);
OutputData=QuickSelector_GUI(InputData);

handles.DataGroup(:,G)=OutputData.NewSelection;
% display immediately in the Group table
for n=1:numel(handles.Data.label)
    handles.DataTable_tab.Data(n,end)={logical(handles.DataGroup(n,G))}; % there might be better option than doing a for cycle, but I haven't found it
end
guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
answer = questdlg('All unsaved changes will be lost. Do you want to save your project?','Save or close','Yes','No','Yes');
if strcmp(answer,'Yes')
    if isfield(handles,'LastFilePath')
        [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project',handles.LastFilePath);
        handles.LastFilePath=PathName;
    else
        [FileName,PathName,~] = uiputfile({'*.psbproj','PSB_GUI project files (*.psbproj)'},'Save your current project');
        handles.LastFilePath=PathName;
    end
    if FileName~=0
        save([PathName FileName],'handles','-mat');
        handles.figure1.Name=['PSB_GUI: ' FileName];
        guidata(hObject, handles);
        % the GUI will be closed only if the saving was not canceled
        delete(hObject);
    end
else
    delete(hObject);
end


% --------------------------------------------------------------------
function About_menu_Callback(hObject, eventdata, handles)
% hObject    handle to About_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg({'PSB_GUI is Matlab GUI application developed for the refinement of residual stress, microstrain and crystallite sizes in cubic materials.';' ';'Version: 1.01';' ';['Created by Petr Cejpek ' char(169) ' 2025']},'About PSB_GUI');
