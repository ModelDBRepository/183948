% runEphapticMSO.m
% Joshua H. Goldwyn [jhgoldwyn@gmail.com]
% July 1, 2015

% Accompanies "Neuronal coupling by endogenous electric fields: Cable theory and applications to coincidence detector neurons in the auditory brainstem" 
% by Goldwyn and Rinzel

%%% THIS CODE CALLS MSO_dae.m %%%
%%% THIS CODE REQUIRES uses sundialsTB. Available at http://computation.llnl.gov/casc/sundials/main.html %%%

% Computes Vm of MSO neuron model 
%          VmTEST of MSO test neuron (with optional spike initiation zone SIZ)
%          Ve extracellular volrage
%
% Ve model is adapted from
% "A model of the medial superior olive explains spatiotemporal features of local field potentials"
% JH Goldwyn, M Mc Laughlin, E Verschooten, PX Joris, J Rinzel
%
% MSO model neuron is adpapted from:
% "Control of submillisecond synaptic timing in binaural coincidence detectors by Kv1 channels"
% Paul J Mathews, Pablo E Jercog,	John Rinzel, Luisa L Scott, Nace L Golding
% Nature Neuroscience 13 601-609 (2010)

close all
clear all

%% MONOLATERAL EXCITATION. NO SIZ INCLUDED. FIGURE 8 IN MANUSCRIPT.
    
%%% Set Parameters %%%
tEnd       = 10;       % simulation duration [ms]
stimPOP    = 'alpha';  % conductance input to population neurons ['alpha' or 'rectsine']
gPOP       = [27 0];   % peak excitatory conductance on left and right dendrite of population neurons [mS / cm2]
freqPOP    = [0 0];    % input frequency (Hz) on left and right dendrite
tsynPOP    = [0 0];    % start time of input to population [ms]
stimTEST   = 'none';   % conductance input to test neuron ['alpha' or 'rectsine']
gTEST      = [0 0];    % peak excitatory conductance on left and right dendrite of test neuron [mS / cm2]
freqTEST   = [0 0];    % input frequency (Hz) on left and right dendrite of test neuron
tsynTEST   = [0 0];    % start time of input to test neuron [ms]
siz        = 0;        % include spike initiation zone in test neuron? [0=No, 1=Yes]
kappa      = []; % SEE BELOW: ephaptic coupling strength in soma and dendrite 

%%% Run model without ephaptic %%%
kappa = [0 0]; 
A = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Run model with ephaptic %%%
kappa = [7 .12]; 
% kappa = [7 .066]; 
B = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);
 
%%% Plot results %%%
FS = 10; LW = 2;
ix = [2 12 22];
figure(8), clf

subplot(2,3,1), hold all, 
    plot(B.t,B.VmPOP(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-60:10:-30,'fontsize',FS)
    set(gca,'position',[.1 .6 .2 .25]); 
    xlim([-.2,3])
    ylim([-65,-27])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m (mV)','fontsize',FS)
    ti=title('A','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.6,0,0])
    text(-1.5,-12,'Example time courses','fontsize',12);
    
subplot(2,3,2), hold all, plot(B.t,B.VmTEST(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-60:1:-58)
    set(gca,'position',[.425 .6 .2 .25]);
    xlim([-.2,3]) 
    ylim([-60.7,-57.8])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m Test (mV)','fontsize',FS)
    ti=title('B','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.9,0,0])
    
subplot(2,3,3), hold all, plot(B.t,B.Ve(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-2:2:2,'fontsize',FS)
    set(gca,'position',[.75 .6 .2 .25]);
    xlim([-.2,3])
    ylim([-2,2.7])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_e (mV)','fontsize',FS)
    ti=title('C','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.5,0,0])

subplot(2,3,4), hold all, 
    plot(A.x,max(A.VmPOP - repmat(A.VmPOP(:,end),1,length(A.t)),[],2),'k','linewidth',LW)
    plot(B.x,max(B.VmPOP - repmat(B.VmPOP(:,end),1,length(B.t)),[],2),'b','linewidth',LW)
    set(gca,'xtick',[-150 0 150],'ytick',0:10:30,'fontsize',FS)
    set(gca,'position',[.1 .1 .2 .25]); 
    xlim([-160,160])
    ylim([0,32])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_m (mV)','fontsize',FS)
    ti=title('D','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-272,0,0])
    text(-295,45,'Maximum deviation from rest','fontsize',12);

ZmTEST = zeros(1,length(B.x));
maxVmTEST = max(B.VmTEST,[],2);  minVmTEST = min(B.VmTEST,[],2); VmTEST0 = B.VmTEST(:,end);
maxVmTEST0 = maxVmTEST - VmTEST0; minVmTEST0 = minVmTEST - VmTEST0;
ZmTEST( maxVmTEST-VmTEST0 >= abs(minVmTEST-VmTEST0) ) = maxVmTEST0( maxVmTEST-VmTEST0 >= abs(minVmTEST-VmTEST0) );
ZmTEST( maxVmTEST-VmTEST0 < abs(minVmTEST-VmTEST0) ) = minVmTEST0( maxVmTEST-VmTEST0 < abs(minVmTEST-VmTEST0) );
subplot(2,3,5), hold all, 
    plot(B.x,zeros(size(B.x)),'k','linewidth',LW)
    plot(B.x,ZmTEST,'linewidth',LW)
    set(gca,'xtick',-150:150:150,'ytick',-1:2)
    set(gca,'position',[.425 .1 .2 .25]);
    xlim([-160,160]) 
    ylim([-1,2])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_m Test (mV)','fontsize',FS)
    ti=title('E','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-295,0,0])

Ze = zeros(1,length(B.x));
maxVe = max(B.Ve,[],2); minVe = min(B.Ve,[],2); Ve0 = B.Ve(:,end);
Ze( maxVe-Ve0 >= abs(minVe-Ve0) ) = maxVe( maxVe-Ve0 >= abs(minVe-Ve0) );
Ze( maxVe-Ve0 < abs(minVe-Ve0) ) = minVe( maxVe-Ve0 < abs(minVe-Ve0) );
subplot(2,3,6), hold all, 
    plot(A.x,A.Ve(:,end),'k','linewidth',LW)
    plot(B.x,Ze,'linewidth',LW)
    set(gca,'xtick',-150:150:150,'ytick',-2:2:2,'fontsize',FS)
    set(gca,'position',[.75 .1 .2 .25]);
    xlim([-160,160])
    ylim([-2,2.7])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_e (mV)','fontsize',FS)
    ti=title('F','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-258,0,0])







%% BILATERAL EXCITATION. NO SIZ INCLUDED. FIGURE 9 IN MANUSCRIPT.

%%% Set Parameters %%%
tEnd       = 5;       % simulation duration [ms]
stimPOP    = 'alpha';  % conductance input to population neurons ['alpha' or 'rectsine']
gPOP       = [27 27];  % peak excitatory conductance on left and right dendrite of population neurons [mS / cm2]
freqPOP    = [0 0];    % input frequency (Hz) on left and right dendrite
tsynPOP    = [0 0];    % start time of input to population [ms]
stimTEST   = 'none';   % conductance input to test neuron ['alpha' or 'rectsine']
gTEST      = [0 0];    % peak excitatory conductance on left and right dendrite of test neuron [mS / cm2]
freqTEST   = [0 0];    % input frequency (Hz) on left and right dendrite of test neuron
tsynTEST   = [0 0];    % start time of input to population [ms]
siz        = 0;        % include spike initiation zone in test neuron? [0=No, 1=Yes]
kappa      = []; % SEE BELOW: ephaptic coupling strength in soma and dendrite 

%%% Run model without ephaptic %%%
kappa = [0 0]; 
A = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Run model with ephaptic %%%
kappa = [7 .12]; 
B = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);
 
%%% Plot results %%%
FS = 10; LW = 2;
ix = [2 12];
figure(9), clf

subplot(2,3,1), hold all, 
    plot(B.t,B.VmPOP(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-60:10:-30,'fontsize',FS)
    set(gca,'position',[.1 .6 .2 .25]); 
    xlim([-.2,3])
    ylim([-65,-27])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m (mV)','fontsize',FS)
    ti=title('A','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.6,0,0])
    text(-1.5,-12,'Example time courses','fontsize',12);
    
subplot(2,3,2), hold all, plot(B.t,B.VmTEST(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-60:1:-58)
    set(gca,'position',[.425 .6 .2 .25]);
    xlim([-.2,3]) 
    ylim([-60.7,-57.8])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m Test (mV)','fontsize',FS)
    ti=title('B','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.9,0,0])
    
subplot(2,3,3), hold all, plot(B.t,B.Ve(ix,:),'linewidth',LW)
    set(gca,'xtick',0:3,'ytick',-2:2:2,'fontsize',FS)
    set(gca,'position',[.75 .6 .2 .25]);
    xlim([-.2,3])
    ylim([-2,2.7])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_e (mV)','fontsize',FS)
    ti=title('C','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-2.5,0,0])

subplot(2,3,4), hold all, 
    plot(A.x,max(A.VmPOP - repmat(A.VmPOP(:,end),1,length(A.t)),[],2),'k','linewidth',LW)
    plot(B.x,max(B.VmPOP - repmat(B.VmPOP(:,end),1,length(B.t)),[],2),'b','linewidth',LW)
    set(gca,'xtick',[-150 0 150],'ytick',0:10:30,'fontsize',FS)
    set(gca,'position',[.1 .1 .2 .25]); 
    xlim([-160,160])
    ylim([0,32])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_m (mV)','fontsize',FS)
    ti=title('D','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-272,0,0])
    text(-295,45,'Maximum deviation from rest','fontsize',12);

ZmTEST = zeros(1,length(B.x));
maxVmTEST = max(B.VmTEST,[],2);  minVmTEST = min(B.VmTEST,[],2); VmTEST0 = B.VmTEST(:,end);
maxVmTEST0 = maxVmTEST - VmTEST0; minVmTEST0 = minVmTEST - VmTEST0;
ZmTEST( maxVmTEST-VmTEST0 >= abs(minVmTEST-VmTEST0) ) = maxVmTEST0( maxVmTEST-VmTEST0 >= abs(minVmTEST-VmTEST0) );
ZmTEST( maxVmTEST-VmTEST0 < abs(minVmTEST-VmTEST0) ) = minVmTEST0( maxVmTEST-VmTEST0 < abs(minVmTEST-VmTEST0) );
subplot(2,3,5), hold all, 
    plot(B.x,zeros(size(B.x)),'k','linewidth',LW)
    plot(B.x,ZmTEST,'linewidth',LW)
    set(gca,'xtick',-150:150:150,'ytick',-1:2)
    set(gca,'position',[.425 .1 .2 .25]);
    xlim([-160,160]) 
    ylim([-1,2])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_m Test (mV)','fontsize',FS)
    ti=title('E','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-295,0,0])

Ze = zeros(1,length(B.x));
maxVe = max(B.Ve,[],2); minVe = min(B.Ve,[],2); Ve0 = B.Ve(:,end);
Ze( maxVe-Ve0 >= abs(minVe-Ve0) ) = maxVe( maxVe-Ve0 >= abs(minVe-Ve0) );
Ze( maxVe-Ve0 < abs(minVe-Ve0) ) = minVe( maxVe-Ve0 < abs(minVe-Ve0) );
subplot(2,3,6), hold all, 
    plot(A.x,A.Ve(:,end),'k','linewidth',LW)
    plot(B.x,Ze,'linewidth',LW)
    set(gca,'xtick',-150:150:150,'ytick',-2:2:2,'fontsize',FS)
    set(gca,'position',[.75 .1 .2 .25]);
    xlim([-160,160])
    ylim([-2,2.7])
    xlabel('x (\mum)','fontsize',FS)
    ylabel('V_e (mV)','fontsize',FS)
    ti=title('F','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-258,0,0])







%% BILATERAL EXCITATION. SIZ INCLUDED. FIGURE 10 IN MANUSCRIPT.

%%% Set Parameters %%%
tEnd       = 30;        % simulation duration [ms]
stimPOP    = 'rectSine';% conductance input to population neurons ['alpha' or 'rectsine']
gPOP       = [20 20];   % peak excitatory conductance on left and right dendrite of population neurons [mS / cm2]
freqPOP    = [200 200]; % input frequency (Hz) on left and right dendrite
tsynPOP    = [0 0];     % start time of input to population [ms]
stimTEST   = 'none';    % conductance input to test neuron ['alpha' or 'rectsine']
gTEST      = [0 0];     % peak excitatory conductance on left and right dendrite of test neuron [mS / cm2]
freqTEST   = [0 0];     % input frequency (Hz) on left and right dendrite of test neuron
tsynTEST   = [0 0];    % start time of input to test neuron [ms]
siz        = [];        % include spike initiation zone in test neuron? [0=No, 1-23=Compartment number]
kappa      = []; % SEE BELOW: ephaptic coupling strength in soma and dendrite 

%%% Run model without ephaptic %%%
kappa = [0 0]; 
siz = 12; % center
A = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Run model with ephaptic -- CENTER SIZ %%%
kappa = [7 .12]; 
siz = 12; % center
B = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);
 
%%% Run model with ephaptic -- OFF-CENTER SIZ %%%
kappa = [7 .12]; 
siz = 20; % off-center
C = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Plot results %%%
FS = 10; LW = 2;
ix = [2 12];
figure(10), clf

subplot(2,3,1), hold all, 
    plot(B.t,B.VmPOP(ix,:),'linewidth',LW)
    set(gca,'xtick',15:5:30,'ytick',-60:10:-30,'fontsize',FS)
    set(gca,'position',[.1 .6 .2 .27]); 
    xlim([14.8,30])
    ylim([-65,-27])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m (mV)','fontsize',FS)
    ti=title('A','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-12.5,0,0])

subplot(2,3,2), hold all
    plot(B.t,B.VmTEST(ix,:),'linewidth',LW)
    set(gca,'xtick',15:5:30,'ytick',-60:1:-58)
    set(gca,'position',[.425 .6 .2 .27]);
    xlim([14.8,30])
    ylim([-60.7,-57.8])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m Test (mV)','fontsize',FS)
    ti=title('B','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-12.5,0,0])

subplot(2,3,3), hold all
    plot(B.t,B.Ve(ix,:),'linewidth',LW)
    set(gca,'xtick',15:5:30,'ytick',[-2 0 2],'fontsize',FS)
    set(gca,'position',[.75 .6 .2 .27]);
    xlim([14.8,30])
    ylim([-2,2.7])
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_e (mV)','fontsize',FS)
    ti=title('C','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-11.5,0,0])

subplot(2,3,4)
    F = imread('../figure/cartoon.jpeg');
    image(F)
    box off
    axis off
    set(gca,'position',[.1 .12 .2 .27]); 
    ti=title('D','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-240,0,0])

subplot(2,3,5), hold all, 
    plot(A.t,A.VmSIZ,'k','linewidth',LW)
    plot(B.t,B.VmSIZ,'color',[0 0.74 0.74],'linewidth',LW)
    plot(C.t,C.VmSIZ,'r','linewidth',LW)
    xlim([14.8,30])
    ylim([-60.2,-58.8])
    set(gca,'xtick',15:5:30,'ytick',-60:.5:-59,'fontsize',FS)
    set(gca,'position',[.425 .12 .2 .27]);
    xlabel('time (ms)','fontsize',FS)
    ylabel('V_m SIZ (mV)','fontsize',FS)
    ti=title('E','fontweight','bold','fontsize',12); 
    set(ti,'position',get(ti,'position')+[-14,0,0])
    
    
    
    

%% EXAMPLE: EPHAPTIC COUPLING PREVENTS SPIKE FOR "CENTER SIZ".  
% Refer to Fig 10F for thresholds, and see note below: re converting gTEST to nS

%%% Set Parameters %%%
tEnd       = 20;        % simulation duration [ms]
stimPOP    = 'rectSine';% conductance input to population neurons ['alpha' or 'rectsine']
gPOP       = [20 20];   % peak excitatory conductance on left and right dendrite of population neurons [mS / cm2]
freqPOP    = [200 200]; % input frequency (Hz) on left and right dendrite
tsynPOP    = [0 0];     % start time of input to population [ms]
stimTEST   = 'alpha';   % conductance input to test neuron ['alpha' or 'rectsine']
gTEST      = [16.5 16.5]; % peak excitatory conductance on left and right dendrite of test neuron [mS / cm2].
freqTEST   = [0 0];     % input frequency (Hz) on left and right dendrite of test neuron
tsynTEST   = [10.1 9.9];     % start time of input to test neuron [ms]
siz        = [];        % include spike initiation zone in test neuron? [0=No, 1-23=Compartment number]
kappa      = []; % SEE BELOW: ephaptic coupling strength in soma and dendrite 

%%% Note: Fig 11 thresholds are in units of nS
%%% ton covert gTEST to nS multiply by compartment surface area in dendrite and scale by 1e6 for nS:   gTEST * (pi * 15e-4cm  * 3.5e-4cm) * 1e6
%%% Run model without ephaptic %%%
kappa = [0 0]; 
siz = 12; % center
A = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Run model with ephaptic -- CENTER SIZ %%%
kappa = [7 .12]; 
siz = 12; % center
B = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);
 
%%% Run model with ephaptic -- OFF-CENTER SIZ %%%
kappa = [7 .12]; 
siz = 20; % off-center
C = MSO_dae(tEnd, stimPOP, gPOP, freqPOP, tsynPOP, stimTEST, gTEST, freqTEST, tsynTEST, siz, kappa);

%%% Plot results %%%
LW = 2;
iSoma = 12;
figure(1), clf, hold all

plot(A.t, A.VmTEST(iSoma,:),'k','linewidth',LW)
plot(B.t, B.VmTEST(iSoma,:),'color',[0 0.74 0.74],'linewidth',LW)
plot(C.t, C.VmTEST(iSoma,:),'r','linewidth',LW)
leg = legend({'No Ve' ,'Center SIZ', 'Off-Center SIZ'},'fontsize',20);
xlabel('Time (ms)','fontsize',20)
ylabel('Vm Test (soma) (mV)','fontsize',20)
set(gca,'fontsize',20)
xlim([9 15])
    
