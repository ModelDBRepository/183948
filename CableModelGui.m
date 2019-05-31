function varargout = CableModelGui(varargin)
% CABLEMODELGUI MATLAB code for CableModelGui.fig
%      CABLEMODELGUI, by itself, creates a new CABLEMODELGUI or raises the existing
%      singleton*.
%
%      H = CABLEMODELGUI returns the handle to a new CABLEMODELGUI or the handle to
%      the existing singleton*.
%
%      CABLEMODELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CABLEMODELGUI.M with the given input arguments.
%
%      CABLEMODELGUI('Property','Value',...) creates a new CABLEMODELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applifed to the GUI before CableModelGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CableModelGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CableModelGui

% Last Modified by GUIDE v2.5 06-Jul-2015 17:22:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CableModelGui_OpeningFcn, ...
                   'gui_OutputFcn',  @CableModelGui_OutputFcn, ...
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


% --- Executes just before CableModelGui is made visible.
function CableModelGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CableModelGui (see VARARGIN)

% Choose default command line output for CableModelGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes10); hold all
axis off
box off

% UIWAIT makes CableModelGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CableModelGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

    d    = str2double(get(handles.edit51,'string'));
    Cm   = str2double(get(handles.edit1,'string'));
    Rm   = str2double(get(handles.edit2,'string'));
    Ri   = str2double(get(handles.edit5,'string'));
    tau  = Rm * Cm * 1e-3; % ms
    lami = sqrt( (Rm/Ri)*((d*1e-4)/4)) * 1e4;
    set(handles.edit19,'String',num2str(.01*round(100*tau))); 
    set(handles.edit18,'String',num2str(round(lami))); 

  
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
    d    = str2double(get(handles.edit51,'string'));
    Cm   = str2double(get(handles.edit1,'string'));
    Rm   = str2double(get(handles.edit2,'string'));
    Ri   = str2double(get(handles.edit5,'string'));
    tau  = Rm * Cm * 1e-3; % ms
    lami = sqrt( (Rm/Ri)*((d*1e-4)/4)) * 1e4;
    set(handles.edit19,'String',num2str(.01*round(100*tau))); 
    set(handles.edit18,'String',num2str(round(lami)));  


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1. %%% SOLVE %%%%
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.togglebutton6,'Value',0);

    FS = 12;
    
    % POPULATION % 
    d    = str2double(get(handles.edit51,'string'));
    Cm   = str2double(get(handles.edit1,'string'));
    Rm   = str2double(get(handles.edit2,'string'));
    Ri   = str2double(get(handles.edit5,'string'));
    tau  = Rm * Cm * 1e-3; % ms
    lami = sqrt( (Rm/Ri)*((d*1e-4)/4)) * 1e4;
    set(handles.edit19,'String',num2str(.01*round(100*tau))); 
    set(handles.edit18,'String',num2str(round(lami))); 

    % TEST NEURON % 
    CmTEST   = str2double(get(handles.edit38,'string'));
    RmTEST   = str2double(get(handles.edit39,'string'));
    RiTEST   = str2double(get(handles.edit41,'string'));
    tauTEST  = RmTEST * CmTEST * 1e-3; %ms
    lamiTEST = sqrt( (RmTEST/RiTEST) * ((d*1e-4)/4) ) * 1e4; % um
    set(handles.edit43,'String',num2str(.01*round(100*tauTEST))); 
    set(handles.edit42,'String',num2str(round(lamiTEST))); 

    % SET PARAMETER STRUCTURE
    P.lami   = lami*1e-4; % cm
    P.tau    = tau * 1e-3;  % [s]
    P.rm     = Rm / (pi*d*1e-4);  % Ohm cm 
    P.L      = str2num(get(handles.edit50,'String'))*1e-4; % Length of cables [cm]   

    P.lamiTEST   = lamiTEST * 1e-4;  % cm
    P.tauTEST    = tauTEST * 1e-3;  % [s]
    P.rmTEST     = RmTEST/ (pi*d*1e-4); % Ohm cm
    P.LTEST      = P.L; % Length of test cable same as population  

    P.DG     = str2num(get(handles.edit48,'String')) * 1e-4; % Distance to ground [cm]   
    P.kappa  = str2num(get(handles.edit49,'String')); % coupling constant

    P.w    = str2num(get(handles.edit44,'String'));  % Hz 
    P.x0   = str2num(get(handles.edit4,'String')) * 1e-4;   % Input position [cm]
    P.Iapp = str2num(get(handles.edit45,'String')) * 1e-3;   %[mA/cm]

    % SOLVE
    [x0,y] = sinusoidBVP(P);
    x =  x0 / P.lami;
    xTEST = x0 / P.lamiTEST;
    
    % PLOT AMPLITUDE
    axes(handles.axes1); hold all
    plot(x,abs ( y(1,:) - y(2,:) ), 'LineWidth',2)
    xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
    ylabel('Vm','fontsize',FS)

    axes(handles.axes2); hold all
    plot(xTEST,abs ( y(3,:) - y(2,:) ), 'LineWidth',2)
    xlabel('$x/\hat{\lambda}$','Interpreter','latex','fontsize',FS); 
    ylabel('Vm TEST','fontsize',FS)
    
    axes(handles.axes3); hold all
    plot(x,abs ( y(2,:) ), 'LineWidth',2)
    xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
    ylabel('Ve','fontsize',FS)

    % PLOT PHASE 
    axes(handles.axes4); hold all
    plot(x,phase ( y(1,:) - y(2,:) ), 'LineWidth',2)
    xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
    ylabel('Vm','fontsize',FS)

    axes(handles.axes5); hold all
    plot(xTEST,phase ( y(3,:) - y(2,:) ), 'LineWidth',2)
    xlabel('$x/\hat{\lambda}$','Interpreter','latex','fontsize',FS); 
    ylabel('Vm TEST','fontsize',FS)
    
    axes(handles.axes6); hold all
    plot(x,phase ( y(2,:) ), 'LineWidth',2)
    xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
    ylabel('Ve','fontsize',FS)

    
    axes(handles.axes10); hold all
    set(handles.text84,'String', num2str(str2num(get(handles.text84,'String'))+4) );
    axis off
    box off
    ylim([0 100])
    xlim([0 1])
    plot([.02 .1],100-[str2num(get(handles.text84,'String')) str2num(get(handles.text84,'String'))],'linewidth',4)
    text(.12,99,['  \lambda           \lambda_{TEST}        \tau        \omega         \kappa']);
    text(.12,100-str2num(get(handles.text84,'String')),...
               [num2str(round(P.lami*1e4)),'\mum  ',...
                num2str(round(P.lamiTEST*1e4)),'\mum    ',...
                num2str(.01*round(100*P.tau*1e3)),'ms   ',...
                num2str(P.w),'Hz    ',...
                num2str(P.kappa)]);
    
    
    % ANIMATIONS
    if P.w>0;
        t = 0 ; 
        dt = 2e-5;
            
        while ~get(handles.togglebutton6,'value')

           t = t + dt;
           w = P.w;

           axes(handles.axes7);
           Vm = real( (y(1,:) - y(2,:)) .*exp (2 * pi * sqrt(-1) * w * t) );
           plot(x,Vm,'k','linewidth',2)
           ylim(1.1*max(abs(y(1,:)-y(2,:)))*[-1,1])
           xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
           ylabel('Vm','fontsize',FS)

           axes(handles.axes8);
           VmTEST = real( (y(3,:) - y(2,:)) .*exp (2 * pi * sqrt(-1) * w * t) );
           plot(xTEST,VmTEST,'k','linewidth',2)
           ylim(1.1*max(abs(y(3,:)-y(2,:)))*[-1,1])
           xlabel('$x/\hat{\lambda}$','Interpreter','latex','fontsize',FS); 
           ylabel('Vm TEST','fontsize',FS)

           axes(handles.axes9);
           Ve = real( y(2,:) .*exp (2 * pi * sqrt(-1) * w * t) );
           plot(x,Ve,'k','linewidth',2)
           ylim(1.1*max(abs(y(2,:)))*[-1,1])
           xlabel('$x/\lambda$','Interpreter','latex','fontsize',FS); 
           ylabel('Ve','fontsize',FS)

           set(handles.text80,'String',num2str(round(t*1e3)));

    %        pause(0.001)


        end
    end
    

    
% --- Executes on button press in pushbutton3.  %%% CLEAR %%%%
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.togglebutton6,'Value',1);
    axes(handles.axes1); cla reset;
    axes(handles.axes2); cla reset;
    axes(handles.axes3); cla reset;
    axes(handles.axes4); cla reset;
    axes(handles.axes5); cla reset;
    axes(handles.axes6); cla reset;
    axes(handles.axes7); cla reset;
    axes(handles.axes8); cla reset;
    axes(handles.axes9); cla reset;
    axes(handles.axes10); cla reset; axis off; box off;     
    set(handles.text84,'String', num2str(0));
    set(handles.text80,'String',num2str(0));


% CABLE MODEL
function [X, Y] = sinusoidBVP(P)

    % spatial domain
%     nx = 1e3;
%     xint = linspace(0,P.L,nx);
    nx = 1e3;
    dx = P.L / (nx-1);
    xint = ([1:nx]-1) * (P.L)/(nx-1);

    % INITIAL GUESS
    solinit = bvpinit(xint,@initfun);

    % SOLVE PROBLEM
    options = bvpset('reltol',1e-6,'abstol',1e-6);
    sol = bvp4c(@odefun,@bcfun,solinit,options,P);
    Sxint = deval(sol,xint);

    % Amplitude profiles
    AmpVi = abs ( Sxint(1,:) );
    AmpVe = abs ( Sxint(2,:) );
    AmpVm = abs ( Sxint(1,:) - Sxint(2,:) ); 
    AmpViTEST = abs ( Sxint(3,:) ); 
    AmpVmTEST = abs ( Sxint(3,:) - Sxint(2,:) ); 

    % Phase Profiles
    PhVi =  phase( Sxint(1,:) );  
    PhVe =  phase( Sxint(2,:) ); 
    PhVm =  phase( Sxint(1,:) -  Sxint(2,:));
    PhViTEST =  phase( Sxint(3,:) );  
    PhVmTEST =  phase( Sxint(3,:) -  Sxint(2,:)); 

    % Output Variables
    X = xint;
    Y = Sxint;


function dY = odefun(x,Y,P)
    dY = [ Y(4), ... % Ai
           Y(5), ... % Ae
           Y(6), ... % Ai TEST
           (1/P.lami^2)      *(  (2*pi*P.tau*P.w*1i + 1)    *(Y(1) - Y(2)) - P.rm*stim(x,P.x0, P.Iapp) ), ... % dAi
           (P.kappa/P.lami^2)*( -(2*pi*P.tau*P.w*1i + 1)    *(Y(1) - Y(2)) + P.rm*stim(x,P.x0, P.Iapp) ), ...  % dAe
           (1/P.lamiTEST^2)  *(  (2*pi*P.tauTEST*P.w*1i + 1)*(Y(3) - Y(2))  )]; % , ... % dAi TEST

      
function res = bcfun(YL,YR,P)
    res = [YL(4) , ... % dAi=0 (L)
           YR(4) , ....% dAi=0 (R)
           P.DG*YL(5) - YL(2) , ... % dAe mixed (L)
           P.DG*YR(5) + YR(2) , ... % dAe mixed (R)
           YL(6) , ... % dAiTEST=0 (L)
           YR(6) ]; % ....% dAiTEST=0 (R)
           

function Y0 = initfun(x)
    Y0 = zeros([6,length(x)]);
      
function Istim = stim(x,x0, I) 
  Istim = 0;
  if abs(x-x0)<.001
      Istim = I;
  end
  
    



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
    d    = str2double(get(handles.edit51,'string'));
    Cm   = str2double(get(handles.edit1,'string'));
    Rm   = str2double(get(handles.edit2,'string'));
    Ri   = str2double(get(handles.edit5,'string'));
    tau  = Rm * Cm * 1e-3; % ms
    lami = sqrt( (Rm/Ri)*((d*1e-4)/4)) * 1e4;
    set(handles.edit19,'String',num2str(.01*round(100*tau))); 
    set(handles.edit18,'String',num2str(round(lami))); 
 


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double

    % TEST NEURON % 
    d    = str2double(get(handles.edit51,'string'));
    CmTEST   = str2double(get(handles.edit38,'string'));
    RmTEST   = str2double(get(handles.edit39,'string'));
    RiTEST   = str2double(get(handles.edit41,'string'));
    tauTEST  = RmTEST * CmTEST * 1e-3; %ms
    lamiTEST = sqrt( (RmTEST/RiTEST) * ((d*1e-4)/4) ) * 1e4; % um
    set(handles.edit43,'String',num2str(.01*round(100*tauTEST))); 
    set(handles.edit42,'String',num2str(round(lamiTEST)));

% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double

    % TEST NEURON % 
    d    = str2double(get(handles.edit51,'string'));
    CmTEST   = str2double(get(handles.edit38,'string'));
    RmTEST   = str2double(get(handles.edit39,'string'));
    RiTEST   = str2double(get(handles.edit41,'string'));
    tauTEST  = RmTEST * CmTEST * 1e-3; %ms
    lamiTEST = sqrt( (RmTEST/RiTEST) * ((d*1e-4)/4) ) * 1e4; % um
    set(handles.edit43,'String',num2str(.01*round(100*tauTEST))); 
    set(handles.edit42,'String',num2str(round(lamiTEST)));
    
    
% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double

    % TEST NEURON % 
    d    = str2double(get(handles.edit51,'string'));
    CmTEST   = str2double(get(handles.edit38,'string'));
    RmTEST   = str2double(get(handles.edit39,'string'));
    RiTEST   = str2double(get(handles.edit41,'string'));
    tauTEST  = RmTEST * CmTEST * 1e-3; %ms
    lamiTEST = sqrt( (RmTEST/RiTEST) * ((d*1e-4)/4) ) * 1e4; % um
    set(handles.edit43,'String',num2str(.01*round(100*tauTEST))); 
    set(handles.edit42,'String',num2str(round(lamiTEST)));

    
% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in togglebutton2.
% function togglebutton2_Callback(hObject, eventdata, handles)
% % hObject    handle to togglebutton2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on button press in togglebutton6.
function togglebutton6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton6



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double

    
    % POPULATION % 
    d    = str2double(get(handles.edit51,'string'));
    Cm   = str2double(get(handles.edit1,'string'));
    Rm   = str2double(get(handles.edit2,'string'));
    Ri   = str2double(get(handles.edit5,'string'));
    tau  = Rm * Cm * 1e-3; % ms
    lami = sqrt( (Rm/Ri)*((d*1e-4)/4)) * 1e4;
    set(handles.edit19,'String',num2str(.01*round(100*tau))); 
    set(handles.edit18,'String',num2str(round(lami))); 

    % TEST NEURON % 
    CmTEST   = str2double(get(handles.edit38,'string'));
    RmTEST   = str2double(get(handles.edit39,'string'));
    RiTEST   = str2double(get(handles.edit41,'string'));
    tauTEST  = RmTEST * CmTEST * 1e-3; %ms
    lamiTEST = sqrt( (RmTEST/RiTEST) * ((d*1e-4)/4) ) * 1e4; % um
    set(handles.edit43,'String',num2str(.01*round(100*tauTEST))); 
    set(handles.edit42,'String',num2str(round(lamiTEST))); 



% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
