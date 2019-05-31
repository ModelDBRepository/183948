% MSO_dae.m
% By: Joshua H. Goldwyn [jhgoldwyn@gmail.com]
% Date: July 1, 2015

% This code is called by runEphapticMSO.m
% Computes VmPOP of MSO neuron in population
%          VmTEST of MSO test neuron
%          Ve extracellular voltage

% Ve model is adapted from:
% "A model of the medial superior olive explains spatiotemporal features of local field potentials"
% JH Goldwyn, M Mc Laughlin, E Verschooten, PX Joris, J Rinzel
%
% MSO model neuron is adapted from:
% "Control of submillisecond synaptic timing in binaural coincidence detectors by Kv1 channels"
% Paul J Mathews, Pablo E Jercog,	John Rinzel, Luisa L Scott, Nace L Golding
% Nature Neuroscience 13 601-609 (2010)

function OutStruct = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa)
    
    % Simulation time
    ParamStruct.t0       = 0;     % ms
    ParamStruct.tEnd     = tEnd;  % ms
    ParamStruct.dt       = 1e-3;  % ms
    
    % Neuron Parameter Values
    ParamStruct.Cm       = 0.9;  % micro F / cm2
    ParamStruct.Ri       = 200.; % axial resistivity [Ohm cm]
    ParamStruct.Glk      = 0.3;  % mS / cm2
    ParamStruct.Vlk      = -60.; % mV used for dendrites and soma
    ParamStruct.VK       = -106.;% mV
    ParamStruct.Vh       = -43.; % mV
    ParamStruct.Vrest    = -59.72;% mV
    
    % Soma Parameters
    ParamStruct.nS       = 3;     % # soma compartments
    ParamStruct.GKLTS    = 17.;   % mS / cm2
    ParamStruct.GhS      = 0.86;  % mS / cm2 
    ParamStruct.diamS    = 20E-4; % soma diameter (cm)
    ParamStruct.lS       = 20.E-4;% soma length (cm)
    ParamStruct.dxS      = ParamStruct.lS/ParamStruct.nS; % soma compartment length (cm)
    ParamStruct.SurfaceS = pi * ParamStruct.diamS * ParamStruct.dxS; % soma compartment surface area (cm2)
    ParamStruct.CrossS   = pi * (ParamStruct.diamS/2.).^2. ; % soma compartment cross section area (cm2)
    
    % Dendrite Parameters
    ParamStruct.nD       = 10;       % # dendrite compartments
    ParamStruct.GKLTD    = 3.57;     % mS / cm2
    ParamStruct.GhD      = .18;      % mS / cm2
    ParamStruct.diamD    = 3.5E-4;   % dendrite diameter (cm)
    ParamStruct.lD       = 150.E-4;  % dendrite length (cm)
    ParamStruct.dxD      = ParamStruct.lD/ParamStruct.nD; % dendrite compartment length (cm)
    ParamStruct.SurfaceD = pi * ParamStruct.diamD * ParamStruct.dxD ;% dendrite compartment surface area (cm2)
    ParamStruct.CrossD   = pi * (ParamStruct.diamD/2.)^2.; % dendrite compartment cross section area (cm2)

    % Spike Initiation Zone Parameters
    ParamStruct.siz      = siz;   % attach SIZ to test neuron soma? (0=NO, 1-23=Compartment Location)
    ParamStruct.iCenter  = 12;    % compartment number of connections between SIZ and soma
    ParamStruct.GlkSIZ   = 200;   % lk conducatance (mS/cm2)
    ParamStruct.VlkSIZ   = -59.81;% lk reversal potential (mV)
    ParamStruct.GNaSIZ   = 75000; % Na max conductance (mS/cm2)
    ParamStruct.VNaSIZ   = 55;    % Na reversal potential (mV)
    ParamStruct.Gax      = 0.00006; % axial conductance between SIZ and soma (mS)
	ParamStruct.SurfaceSIZ = 3.14e-8; % surface area of SIZ (cm2)
    
    % Synapse parameters
    ParamStruct.tausyn   = 0.2;  % synaptic time constant (ms) [for alpha]
    ParamStruct.Vsyn     = 0;  % excitation reversal potential (mV)
    ParamStruct.csyn     = [2 22]; % compartments receiving synaptic inputs
    ParamStruct.stimPOP  = stimPOP; % stimulus type for MSO population ('alpha' or 'rectSine')
    ParamStruct.gPOP     = gPOP;  % peak input conductance density for population neuron (mS/cm2)
    ParamStruct.freqPOP  = freqPOP; % input frequency for population neuron (Hz)
    ParamStruct.tsynPOP  = tsynPOP; % Start time of input to population (ms)
    ParamStruct.stimTEST = stimTEST; % stimulus type for test neuron ('alpha' or 'rectSine')
    ParamStruct.gTEST    = gTEST; % peak input conductance density for test neuron (mS/cm2)
    ParamStruct.freqTEST = freqTEST;% input frequency for test neuron (Hz)
    ParamStruct.tsynTEST = tsynTEST; % Start time of input to test neuron (ms)
    
    % Extracellular parameters
    ParamStruct.dxG      = 0.1; % Distance to ground (cm)
    ParamStruct.kappaS   = kappa(1); % Ephaptic coupling constant in soma
    ParamStruct.kappaD   = kappa(2); % Ephaptic coupling constant in dendrite
    ParamStruct.CrossG   = 450e-8; %(cm2)
    
    % Compartment measurements
    ParamStruct.nC       = 2*ParamStruct.nD + ParamStruct.nS;
    ParamStruct.GKLT     = [ParamStruct.GKLTD+zeros(1,ParamStruct.nD)     , ParamStruct.GKLTS+zeros(1,ParamStruct.nS)     , ParamStruct.GKLTD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Gh       = [ParamStruct.GhD+zeros(1,ParamStruct.nD)       , ParamStruct.GhS+zeros(1,ParamStruct.nS)       , ParamStruct.GhD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.dx       = [ParamStruct.dxD+zeros(1,ParamStruct.nD)       , ParamStruct.dxS+zeros(1,ParamStruct.nS)       , ParamStruct.dxD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.diam     = [ParamStruct.diamD+zeros(1,ParamStruct.nD)     , ParamStruct.diamS+zeros(1,ParamStruct.nS)     , ParamStruct.diamD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Cross    = [ParamStruct.CrossD+zeros(1,ParamStruct.nD)    , ParamStruct.CrossS+zeros(1,ParamStruct.nS)    , ParamStruct.CrossD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.Surface  = [ParamStruct.SurfaceD+zeros(1,ParamStruct.nD)  , ParamStruct.SurfaceS+zeros(1,ParamStruct.nS)  , ParamStruct.SurfaceD+zeros(1,ParamStruct.nD) ]';
    ParamStruct.kappa    = [ParamStruct.kappaD+zeros(1,ParamStruct.nD)    , ParamStruct.kappaS+zeros(1,ParamStruct.nS)    , ParamStruct.kappaD+zeros(1,ParamStruct.nD) ]';
    
    %%%%%%%%%% RUN SOLVER %%%%%%%%%%
    [t,y] = DAESolveIt(ParamStruct);
    
    % Get results
    ViPOP   = y(1:ParamStruct.nC,:);       
    wPOP    = y(ParamStruct.nC+1:2*ParamStruct.nC,:);
    zPOP    = y(2*ParamStruct.nC+1:3*ParamStruct.nC,:);   
    ViTEST  = y(3*ParamStruct.nC+1:4*ParamStruct.nC,:);
    wTEST   = y(4*ParamStruct.nC+1:5*ParamStruct.nC,:);
    zTEST   = y(5*ParamStruct.nC+1:6*ParamStruct.nC,:);
    Ve      = y(6*ParamStruct.nC+1:7*ParamStruct.nC,:);       
    VmSIZ   = y(7*ParamStruct.nC+1,:);
    mSIZ    = y(7*ParamStruct.nC+2,:);
    hSIZ    = y(7*ParamStruct.nC+3,:);

    VmPOP   = ViPOP  - Ve;
    VmTEST  = ViTEST - Ve;
    
    % center locations of compartments (micro m)
    x = 1e4 * [-(ParamStruct.lD+ParamStruct.lS/2)+ParamStruct.dxD/2+[0:ParamStruct.nD-1]*ParamStruct.dxD ...
                ParamStruct.dxS*[-1 0 1] ...
                ParamStruct.lS/2+ParamStruct.dxD/2+[0:ParamStruct.nD-1]*ParamStruct.dxD];
    
    % OUTPUT
    OutStruct.ParamStruct = ParamStruct;
    OutStruct.t           = t;
    OutStruct.x           = x;
    OutStruct.VmPOP       = VmPOP;
    OutStruct.VmTEST      = VmTEST;
    OutStruct.VmSIZ       = VmSIZ;
    OutStruct.Ve          = Ve;


% GATING FUNCTIONS
function out   = gate_func()

    % KLT gating variables
	out.winf   = @(V) 1 ./(1.+exp((V-(-57.34))/(-11.7)));
    out.tauw   = @(V) 21.5 ./ (6.*exp((V+60.)/7.) + 24.*exp(-(V+60.)/50.6)) + 0.35;
    out.zinf   = @(V) (1.-.27) ./ (1 + exp((V-(-67))/6.16)) + .27;
    out.tauz   = @(V) 170 ./ (5*exp((V+60.)/10.) + exp(-(V+70.)/8.)) + 10.7;
    
    % Na gating variables
	out.minf   = @(V) 1 ./ (1.+exp(-(V+38.)/7.));
    out.taum   = @(V) 0.24 .* (10 ./ (5*exp((V+60.)/18.)+36.*exp(-(V+60.)/25.))+0.04);
    out.hinf   = @(V) 1 ./ (1.+exp((V+71)/6.)); % left shifted 6mV
    out.tauh   = @(V) 0.24 .* (100 ./ (7.*exp((V+66.)/11.)+10.*exp(-(V+66.)/25.))+0.6); % left shifted 6mV
    
   
function [t,y] = DAESolveIt(ParamStruct)

    % Gating functions
    gate = gate_func;

    % number of variables
    nVar = 7*ParamStruct.nC+3;
    id = ones(nVar,1);

    % IDA solver options
    options = IDASetOptions('RelTol',1.e-6,...
                            'AbsTol',1.e-8,...
                            'VariableTypes',id,...
                            ...%  'MaxStep', [],...
                            ...%'MaxNumSteps', [],...
                            'UserData', ParamStruct); % Structure of parameter values

    % Initialize variables
    Vrest   = ParamStruct.Vrest*ones(ParamStruct.nC,1);%-58; % Resting potential (mV)
    ViPOP0  = Vrest;
    ViTEST0 = Vrest;
    Ve0     = zeros(size(Vrest));
    wPOP0   = gate.winf(Vrest);
    wTEST0  = gate.winf(Vrest);
    zPOP0   = gate.zinf(Vrest);
    zTEST0  = gate.zinf(Vrest);
    VmSIZ0  = ParamStruct.Vrest;
    m0      = gate.minf(ParamStruct.Vrest);
    h0      = gate.hinf(ParamStruct.Vrest);
    y0 = [ViPOP0  ; wPOP0  ; zPOP0  ; ... % POPULATION NEURON
          ViTEST0 ; wTEST0 ; zTEST0 ; ... % TEST NEURON
          Ve0 ;                       ... % EXTRACELLULAR VOLTAGE
          VmSIZ0 ; m0 ; h0];              % SIZ
    yp0 = zeros(size(y0)); % derivatives

    % Initialize solver
    IDAInit(@res_dae,ParamStruct.t0,y0,yp0,options);
    [status, y0_mod, yp0_mod] = IDACalcIC(ParamStruct.dt, 'FindAlgebraic');

    % Run IDA solver
    [status,t, y] = IDASolve([ParamStruct.dt:ParamStruct.dt:ParamStruct.tEnd],'Normal');
 
    % output vectors, including initial values
    t = [ParamStruct.t0 t];
    y = [y0 y];
 
    % Deallocate IDAS memory
    IDAFree


function [res, flag, new_data] = res_dae(t,y,yp,ParamStruct) 
    % Solve equations:
    % Vmp: Cm*Vmp - (Ilk + Iion + Isyn)) - IlongI =0
    % Vmp: Cm*Vmp - (Ilk + Iion + Isyn)) + IlongE =0
    % VmTEST: Cm*VmpTEST - (Iion + Isyn)) + IlongE =0
    % Gating variables: xp - (xinf-x)/taux = 0 [x=m,h,w,z]

    % Gating functions
    gate = gate_func;

    % Variables [Dendrite Soma Dendrite]
    ViPOP   = y(1:ParamStruct.nC);       
    wPOP    = y(ParamStruct.nC+1:2*ParamStruct.nC);
    zPOP    = y(2*ParamStruct.nC+1:3*ParamStruct.nC);   
    ViTEST  = y(3*ParamStruct.nC+1:4*ParamStruct.nC);
    wTEST   = y(4*ParamStruct.nC+1:5*ParamStruct.nC);
    zTEST   = y(5*ParamStruct.nC+1:6*ParamStruct.nC);
    Ve      = y(6*ParamStruct.nC+1:7*ParamStruct.nC);       
    VmSIZ   = y(7*ParamStruct.nC+1);
    mSIZ    = y(7*ParamStruct.nC+2);
    hSIZ    = y(7*ParamStruct.nC+3);
    
    % Derivatives
    ViPOPp  = yp(1:ParamStruct.nC);       
    wPOPp   = yp(ParamStruct.nC+1:2*ParamStruct.nC);
    zPOPp   = yp(2*ParamStruct.nC+1:3*ParamStruct.nC);   
    ViTESTp = yp(3*ParamStruct.nC+1:4*ParamStruct.nC);
    wTESTp  = yp(4*ParamStruct.nC+1:5*ParamStruct.nC);
    zTESTp  = yp(5*ParamStruct.nC+1:6*ParamStruct.nC);
    Vep     = yp(6*ParamStruct.nC+1:7*ParamStruct.nC);         
    VmSIZp  = yp(7*ParamStruct.nC+1);
    mSIZp   = yp(7*ParamStruct.nC+2);
    hSIZp   = yp(7*ParamStruct.nC+3);     

    % Membrane potential (mV)
    VmPOP   = ViPOP  - Ve;
    VmPOPp  = ViPOPp - Vep;
    VmTEST  = ViTEST  - Ve;
    VmTESTp = ViTESTp - Vep;

    % Population neuron membrane current per unit length (mA / cm)
    IlkPOP  = (ParamStruct.Glk/1000)  .* pi .* ParamStruct.diam .* (VmPOP - ParamStruct.Vlk);
    IhPOP   = (ParamStruct.Gh/1000)   .* pi .* ParamStruct.diam .* (VmPOP - ParamStruct.Vh);
    IkltPOP = (ParamStruct.GKLT/1000) .* pi .* ParamStruct.diam .* wPOP.^4. .* zPOP.* (VmPOP - ParamStruct.VK); 

    % Test neuron membrane current per unit length (mA / cm)
    IlkTEST  = (ParamStruct.Glk/1000)  .* pi .* ParamStruct.diam .* (VmTEST - ParamStruct.Vlk);
    IhTEST   = (ParamStruct.Gh/1000)   .* pi .* ParamStruct.diam .* (VmTEST - ParamStruct.Vh);
    IkltTEST = (ParamStruct.GKLT/1000) .* pi .* ParamStruct.diam .* wTEST.^4 .* zTEST.* (VmTEST - ParamStruct.VK); 

    % SIZ membrane current density per unit length (mA / cm)
    IlkSIZ   = (ParamStruct.GlkSIZ /1000.) * (VmSIZ-ParamStruct.VlkSIZ);
    INaSIZ   = (ParamStruct.GNaSIZ /1000.) * mSIZ.^3 * hSIZ *(VmSIZ-ParamStruct.VNaSIZ);
    
    % Intracellular (axial) current in population neuron (mA)
    nC = ParamStruct.nC; dx = ParamStruct.dx; Cross = ParamStruct.Cross;
    IinPOP = zeros(size(IlkPOP));
    IinPOP(1)      = (1/ParamStruct.Ri) * ( ViPOP(2)-ViPOP(1)) / (dx(1)/Cross(1));    % left dendrite boundary
    IinPOP(2:nC-1) = (1/ParamStruct.Ri) * ( (ViPOP(3:nC)  -ViPOP(2:nC-1))./((dx(3:nC)./Cross(3:nC)+dx(2:nC-1)./Cross(2:nC-1))/2)  + ...
                                            (ViPOP(1:nC-2)-ViPOP(2:nC-1))./((dx(2:nC-1)./Cross(2:nC-1)+dx(1:nC-2)./Cross(1:nC-2))/2) ); 
    IinPOP(nC)     = (1/ParamStruct.Ri) * ( ViPOP(nC-1)-ViPOP(nC) ) / ( dx(nC)/Cross(nC) ) ; % Right dendrite boundary 
    
    % Intracellular (axial) current in test neuron (mA)
    IinTEST = zeros(size(IlkTEST));
    IinTEST(1)      = (1/ParamStruct.Ri) * ( ViTEST(2)-ViTEST(1)) / (dx(1)/Cross(1));    % left dendrite boundary
    IinTEST(2:nC-1) = (1/ParamStruct.Ri) * ( (ViTEST(3:nC)-ViTEST(2:nC-1))./((dx(3:nC)./Cross(3:nC)+dx(2:nC-1)./Cross(2:nC-1))/2)  + ...
                                             (ViTEST(1:nC-2)-ViTEST(2:nC-1))./((dx(2:nC-1)./Cross(2:nC-1)+dx(1:nC-2)./Cross(1:nC-2))/2) ); 
    IinTEST(nC)     = (1/ParamStruct.Ri) * ( ViTEST(nC-1)-ViTEST(nC) ) / ( dx(nC)/Cross(nC) ) ; % Right dendrite boundary 

    % Axial current between SIZ and soma in test neuron (mA)
    if ParamStruct.siz>0
        ViSIZ           = VmSIZ + Ve(ParamStruct.siz); % Vi in SIZ. Ve depends on SIZ location 
    else
        ViSIZ = VmSIZ;
    end
    
    IaxSIZ          = zeros(size(IlkTEST));
    if ParamStruct.siz>0 % attach SIZ to test neuron soma
        IaxSIZ(ParamStruct.iCenter) = (ParamStruct.Gax/1000.) * (ViTEST(ParamStruct.iCenter)-ViSIZ) ; 
    end
    
    % Exctracellular current (mA)
    Iex = zeros(size(IlkPOP));
    kappa    = ParamStruct.kappa;
    kappa(kappa==0) = 1e-12; % make kappa small but non-zero to avoid divide by 0 error
    Iex(1)          = (1/ParamStruct.Ri) * ( (Ve(2)-Ve(1))/(kappa(1)*dx(1)/Cross(1))    +   (0-Ve(1))/(kappa(1)*ParamStruct.dxG/Cross(1) + (kappa(1)*dx(1)/Cross(1))/2)   ); % Left dendrite boundary
    Iex(2:nC-1)     = (1/ParamStruct.Ri) * ( ( Ve(3:nC)  -Ve(2:nC-1) ) ./ ( (kappa(3:nC).*dx(3:nC)./Cross(3:nC)+kappa(2:nC-1).*dx(2:nC-1)./Cross(2:nC-1))/2 ) + ...
                                             ( Ve(1:nC-2)-Ve(2:nC-1) ) ./ ( (kappa(1:nC-2).*dx(1:nC-2)./Cross(1:nC-2)+kappa(2:nC-1).*dx(2:nC-1)./Cross(2:nC-1))/2 )  );
    Iex(nC)         = (1/ParamStruct.Ri) * ( ( 0-Ve(nC) ) / ( kappa(nC)*ParamStruct.dxG/Cross(nC) + (kappa(nC)*dx(nC)/Cross(nC))/2 ) + ( Ve(nC-1)-Ve(nC) ) / ( kappa(nC)*dx(nC)/Cross(nC))  ); % Right dendrite boundary

    % Input current to population (mA/cm2)
    IstimPOP  = 0*IlkPOP;
    tt1 = t - ParamStruct.tsynPOP(1);
    tt2 = t - ParamStruct.tsynPOP(2);
    switch(ParamStruct.stimPOP)
        case('alpha')
            IstimPOP(ParamStruct.csyn(1)) =  (tt1>0)*(ParamStruct.gPOP(1)/1000)*pi*ParamStruct.diamD* (tt1/ParamStruct.tausyn)*exp(1-tt1/ParamStruct.tausyn) * (VmPOP(ParamStruct.csyn(1)) - ParamStruct.Vsyn);
            IstimPOP(ParamStruct.csyn(2)) =  (tt2>0)*(ParamStruct.gPOP(2)/1000)*pi*ParamStruct.diamD* (tt2/ParamStruct.tausyn)*exp(1-tt2/ParamStruct.tausyn) * (VmPOP(ParamStruct.csyn(2)) - ParamStruct.Vsyn);
        case('rectSine')
            IstimPOP(ParamStruct.csyn(1)) =  (tt1>0)*(ParamStruct.gPOP(1)/1000)*pi*ParamStruct.diamD*max(0, sin(2*pi*(tt1/1000)*ParamStruct.freqPOP(1)))* (VmPOP(ParamStruct.csyn(1)) - ParamStruct.Vsyn);
            IstimPOP(ParamStruct.csyn(2)) =  (tt2>0)*(ParamStruct.gPOP(2)/1000)*pi*ParamStruct.diamD*max(0, sin(2*pi*(tt2/1000)*ParamStruct.freqPOP(2)))* (VmPOP(ParamStruct.csyn(2)) - ParamStruct.Vsyn);
    end
    
    % Input current to test neuron (mA/cm2)
    IstimTEST = 0*IlkTEST;
    tt1 = t - ParamStruct.tsynTEST(1);
    tt2 = t - ParamStruct.tsynTEST(2);
    switch(ParamStruct.stimTEST)
        case('alpha')
            IstimTEST(ParamStruct.csyn(1)) =  (tt1>0)*(ParamStruct.gTEST(1)/1000)*pi*ParamStruct.diamD* (tt1/ParamStruct.tausyn)*exp(1-tt1/ParamStruct.tausyn) * (VmTEST(ParamStruct.csyn(1)) - ParamStruct.Vsyn);
            IstimTEST(ParamStruct.csyn(2)) =  (tt2>0)*(ParamStruct.gTEST(2)/1000)*pi*ParamStruct.diamD* (tt2/ParamStruct.tausyn)*exp(1-tt2/ParamStruct.tausyn) * (VmTEST(ParamStruct.csyn(2)) - ParamStruct.Vsyn);
        case('rectSine')
            IstimTEST(ParamStruct.csyn(1)) =  (tt1>0)*(ParamStruct.gTEST(1)/1000)*pi*ParamStruct.diamD*max(0, sin(2*pi*(tt1/1000)*ParamStruct.freqTEST(1)))* (VmTEST(ParamStruct.csyn(1)) - ParamStruct.Vsyn);
            IstimTEST(ParamStruct.csyn(2)) =  (tt2>0)*(ParamStruct.gTEST(2)/1000)*pi*ParamStruct.diamD*max(0, sin(2*pi*(tt2/1000)*ParamStruct.freqTEST(2)))* (VmTEST(ParamStruct.csyn(2)) - ParamStruct.Vsyn);
    end

    % residuals for DAE solver
    dx = ParamStruct.dx;
    res =   [dx.* ( (ParamStruct.Cm/1000.).*pi.*ParamStruct.diam.*VmPOPp + IkltPOP + IhPOP + IlkPOP + IstimPOP ) -  IinPOP
             gate.tauw(VmPOP).*wPOPp   - (gate.winf(VmPOP) - wPOP )
             gate.tauz(VmPOP).*zPOPp   - (gate.zinf(VmPOP) - zPOP )
             dx.* ( (ParamStruct.Cm/1000.).*pi.*ParamStruct.diam.*VmTESTp + IkltTEST + IhTEST + IlkTEST + IstimTEST ) -  IinTEST + IaxSIZ
             gate.tauw(VmTEST).*wTESTp - (gate.winf(VmTEST) - wTEST )
             gate.tauz(VmTEST).*zTESTp - (gate.zinf(VmTEST) - zTEST )
             dx.*( (ParamStruct.Cm/1000.).*pi.*ParamStruct.diam.*VmPOPp + IkltPOP + IhPOP + IlkPOP + IstimPOP ) + Iex
             ParamStruct.SurfaceSIZ*(ParamStruct.Cm/1000. * VmSIZp +  INaSIZ + IlkSIZ) - IaxSIZ(ParamStruct.iCenter)
             gate.taum(VmSIZ)*mSIZp    - (gate.minf(VmSIZ) - mSIZ)
             gate.tauh(VmSIZ)*hSIZp    - (gate.hinf(VmSIZ) - hSIZ) ];
    flag = 0;
    new_data = [];

   
    
    