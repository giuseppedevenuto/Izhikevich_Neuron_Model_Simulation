%% Computational Project 
% © De Venuto Giuseppe, Impollonia Roberta 2023

close all
clear
clc

%% Simulation parameters

Tdur = 700;                                             % duration of the simulation [ms]
dt = 0.1;                                               % integration time stamp [ms] (time resolution)
t = 0:dt:Tdur;                                          % vector of time-points
ton = 50;                                               % time to begin applied current (onset)
step_length=500;                                        % length of the current step applied [ms]
toff = ton+step_length;                                 % time to end applied current (offset)
stim_time=toff-ton;                                     % duration of the stimulation [ms]
non = round(ton/dt);                                    % time-point index of current onset
noff = round(toff/dt);                                  % time-point index of current offset

E_rest = -70;                                           % resting potential [mV]
Vspike = 30;                                            % spike peak amplitude [mV]

Ntrials = 401;                                          % number of trials per model = number of applied current steps 
Istim_max_ini=20;                                       % maximal current step amplitude applied [pA]
Istim_min_ini=0;                                        % minimal current step amplitude applied [pA]
dI_ini=(Istim_max_ini-Istim_min_ini)/(Ntrials-1);       % amplitude interval between two consecutive current step [pA] (current resolution)
Istim = Istim_min_ini:dI_ini:Istim_max_ini;             % vector of the amplitudes of the applied current steps [pA]
Istim_in=Istim'*ones(1,noff+1-non);                     % matrix with the vectors of current of the different trials between onset and offset of the stimulation (amplitude of the step)

f=0;                                                    % figure index

%% (a) 
% Simulate the model and plot the shape/trend of membrane potential Vmem as a function of the amplitude 
% of the current clamp (dc stimulation). Consider the two working conditions, of under-threshold and 
% supra-threshold.

% Regular Spiking Neuron model
[~,Vout2PlotRS, SpikesCountRS, IstimOutRS] = Izhikevich_RS(dt, E_rest, Vspike, Istim_in,t,non,noff,Ntrials);

% Plot of the trend of the membrane potential in under-threshold  and
% supra-threshold condition for the Regular Spiking neuron model
f=f+1;
figure(f)
set(gcf,'position',[0 100 2000 350])
axes('position',[0.03 0.13 0.93 0.80])
yyaxis left
plot(t,IstimOutRS(201,:))           % supra-threshold applied current plot 
hold on
plot(t,IstimOutRS(41,:))            % under-threshold applied current plot 
ylabel('Applied Current [pA]')
ylim([-1.5 80])
hold off
text(ton-40, Istim(201),sprintf('I_a_p_p: %.2f pA',Istim(201)),FontSize=8)
text(ton-40, Istim(61),sprintf('I_a_p_p: %.2f pA',Istim(41)),FontSize=8)
yyaxis right
plot(t,Vout2PlotRS(201,:))          % supra-threshold voltage plot
hold on
plot(t,Vout2PlotRS(41,:))           % under-threshold voltage plot
ylabel('Membrane Potential [mV]')
ylim([-100 32.5])
hold off
title('Under-threshold & Supra-threshold conditions for RS neuron')
xlabel('time [ms]')
box off
legend('Supra-threshold stimulus','Under-threshold stimulus','Supra-threshold mem. pot. behaviour','Under-threshold mem. pot. behaviour')
legend box off

% Chattering Neuron model
[~,Vout2PlotCH, SpikesCountCH, IstimOutCH] = Izhikevich_CH(dt, E_rest, Vspike, Istim_in,t,non,noff,Ntrials);

% Plot of the trend of the membrane potential in under-threshold  and
% supra-threshold condition for the Chattering neuron model
f=f+1;
figure(f)
set(gcf,'position',[0 100 2000 350])
axes('position',[0.03 0.13 0.93 0.80])
yyaxis left
plot(t,IstimOutCH(201,:))        % supra-threshold applied current plot
hold on
plot(t,IstimOutCH(41,:))        % under-threshold applied current plot
ylabel('Applied Current [pA]') 
ylim([-1.5 80])
hold off
text(ton-40, Istim(201),sprintf('I_a_p_p: %.2f pA',Istim(201)),FontSize=8)
text(ton-40, Istim(61),sprintf('I_a_p_p: %.2f pA',Istim(41)),FontSize=8)
yyaxis right
plot(t,Vout2PlotCH(201,:))       % supra-threshold voltage plot
hold on
plot(t,Vout2PlotCH(41,:))       % under-threshold voltage plot
ylabel('Membrane Potential [mV]')
ylim([-100 32.5])
hold off
title('Under-threshold & Supra-threshold conditions for CH neuron')
xlabel('time [ms]')
box off
legend('Supra-threshold stimulus','Under-threshold stimulus','Supra-threshold mem. pot. behaviour','Under-threshold mem. pot. behaviour')
legend box off

%% (b)
% Determine for the two models the correspondent voltages threshold VTH.
 
n_cycles=6;                                         % number of iterations for the "voltage threshold-finding method"
                                    
% Regular Spiking neuron case
n_of_spikesRS=sum(SpikesCountRS,2);                 % number of the spikes per trial
ThRS_idx=find(n_of_spikesRS, 1 );                   % find the index of the first current(= minimum current) that make the neuron fire
IThRS=Istim(ThRS_idx);                              % current threshold for spiking
fprintf('0) %.2f \n',IThRS)
VThRS=max(Vout2PlotRS(ThRS_idx-1,:));               % voltage threshold for spiking
fprintf('0) %.2f \n',VThRS)

% To look to a narrow current window with more resolution
dI_cycle=dI_ini;                                    % initialization of the amplitude interval between two consecutive current step in the iterative threshold-finding method [pA]
IThRS_cycle=IThRS;                                  % initialization of the current threshold from which the iterative process starts
for i=1:n_cycles
    Istim_max=IThRS_cycle;                         % update the maximal current step amplitude applied to use in the iterative process [pA]
    Istim_min=IThRS_cycle-dI_cycle;                         % update the minimal current step amplitude applied to use in the iterative process [pA]
    dI_cycle=(Istim_max-Istim_min)/Ntrials;                 % update the amplitude interval between two consecutive current step in the iterative process [pA] (current resolution)
    Istim_cycle = Istim_min:dI_cycle:Istim_max;             % vector of the amplitudes of the applied current steps [pA]
    Istim_in_cycle=Istim_cycle'*ones(1,noff+1-non);         % matrix with the vectors of current of the different trials between onset and offset of the stimulation (amplitude of the step)

    [~,Vout2PlotRS_cycle, SpikesCountRS_cycle, IstimOutRS_cycle] = Izhikevich_RS(dt, E_rest, Vspike, Istim_in_cycle,t,non,noff,Ntrials);
    n_of_spikesRS_cycle=sum(SpikesCountRS_cycle,2);         % update the number of the spikes per trial
    ThRS_idx_cycle=find(n_of_spikesRS_cycle, 1 );           % update the index of the first current(= minimum current) that make the neuron fire
    IThRS_cycle=Istim_cycle(ThRS_idx_cycle);                % update current threshold for spiking
    fprintf('%d) %.2f \n',i,IThRS_cycle)
    VThRS=max(Vout2PlotRS_cycle(ThRS_idx_cycle-1,:));       % update voltage threshold for spiking
    fprintf('%d) %.2f \n',i,VThRS)
end

% Plot of the trend of the membrane potential in under-threshold and
% supra-threshold condition for the Regular Spiking neuron model with the
% reference of voltage threshold and current threshold
f=f+1;
figure(f)
set(gcf,'position',[0 100 1400 350])
axes('position',[0.03 0.13 0.93 0.80])
yyaxis left
plot(t,IstimOutRS(201,:))                        % supra-threshold applied current plot
hold on
plot(t,IstimOutRS(41,:))                        % under-threshold applied current plot
plot(t,IThRS_cycle*ones(size(t)),LineWidth=1)   % current threshold plot
ylabel('Applied Current [pA]') 
ylim([-1.5 80])
hold off
text(ton-40, Istim(201),sprintf('I_a_p_p: %.2f pA',Istim(201)),FontSize=8)
text(ton-40, Istim(61),sprintf('I_a_p_p: %.2f pA',Istim(41)),FontSize=8)
text(toff+100, IThRS_cycle+2,sprintf('I_t_h: %.2f pA',IThRS_cycle),FontSize=8)
yyaxis right
plot(t,Vout2PlotRS(201,:))                       % supra-threshold voltage plot
hold on
plot(t,Vout2PlotRS(41,:))                       % under-threshold voltage plot
plot(t,VThRS*ones(size(t)),LineWidth=1)         % voltage threshold plot
ylabel('Membrane Potential [mV]')
ylim([-100 32.5])
hold off
text(toff+100, VThRS+4,sprintf('V_t_h: %.2f mV',VThRS),FontSize=8)
title('Under-threshold & Supra-threshold conditions for RS neuron')
xlabel('time [ms]')
box off
legend(sprintf('Supra-threshold current: %.2f pA',Istim(201)),sprintf('Under-threshold current: %.2f pA',Istim(41)), ...
    'Current threshold','Supra-threshold mem. pot. behaviour','Under-threshold mem. pot. behaviour','Voltage threshold')
legend box off

% Chattering Spiking neuron case 
% (for a commented code see the Regular Firing neuron case appling the necessary 
% consideration(the code is now referred to a CH neuron) )
n_of_spikesCH=sum(SpikesCountCH,2);             
ThCH_idx=find(n_of_spikesCH, 1 );
IThCH=Istim(ThCH_idx);
fprintf('0) %.2f \n',IThCH)
VThCH=max(Vout2PlotCH(ThCH_idx-1,:));
fprintf('0) %.2f \n',VThCH)

dI_cycle=dI_ini;
IThCH_cycle=IThCH;
for i=1:n_cycles
    Istim_max=IThCH_cycle;
    Istim_min=IThCH_cycle-dI_cycle;
    dI_cycle=(Istim_max-Istim_min)/Ntrials;
    Istim_cycle = Istim_min:dI_cycle:Istim_max;              
    Istim_in_cycle=Istim_cycle'*ones(1,noff+1-non);

    [~,Vout2PlotCH_cycle, SpikesCountCH_cycle, IstimOutCH_cycle] = Izhikevich_CH(dt, E_rest, Vspike, Istim_in_cycle,t,non,noff,Ntrials);
    n_of_spikesCH_cycle=sum(SpikesCountCH_cycle,2);
    ThCH_idx_cycle=find(n_of_spikesCH_cycle, 1 );
    IThCH_cycle=Istim_cycle(ThCH_idx_cycle);
    fprintf('%d) %.2f \n',i,IThCH_cycle)
    VThCH=max(Vout2PlotCH_cycle(ThCH_idx_cycle-1,:));
    fprintf('%d) %.2f \n',i,VThCH)
end

f=f+1;
figure(f)
set(gcf,'position',[0 100 1400 350])
axes('position',[0.03 0.13 0.93 0.80])
yyaxis left
plot(t,IstimOutCH(201,:))
hold on
plot(t,IstimOutCH(41,:))
plot(t,IThCH_cycle*ones(size(t)),LineWidth=1)
ylabel('Applied Current [pA]') 
ylim([-1.5 80])
hold off
text(ton-40, Istim(201),sprintf('I_a_p_p: %.2f pA',Istim(201)),FontSize=8)
text(ton-40, Istim(61),sprintf('I_a_p_p: %.2f pA',Istim(41)),FontSize=8)
text(toff+100, IThCH_cycle+2,sprintf('I_t_h: %.2f pA',IThCH_cycle),FontSize=8)
yyaxis right
plot(t,Vout2PlotCH(201,:))
hold on
plot(t,Vout2PlotCH(41,:))
plot(t,VThCH*ones(size(t)),LineWidth=1)
ylabel('Membrane Potential [mV]')
ylim([-100 32.5])
hold off
text(toff+100, VThCH+4,sprintf('V_t_h: %.2f mV',VThCH),FontSize=8)
title('Under-threshold & Supra-threshold conditions for CH neuron')
xlabel('time [ms]')
box off
legend(sprintf('Supra-threshold current: %.2f pA',Istim(201)),sprintf('Under-threshold current: %.2f pA',Istim(41)), ...
    'Current threshold','Supra-threshold mem. pot. behaviour','Under-threshold mem. pot. behaviour','Voltage threshold')
legend box off

%% (c)
% Derive the gain function of the two models, by applying dc current pulses with different amplitudes. 
% In which class the two neurons can be identified? Are the developed models able to code low-frequency 
% neuronal activity? Are significant differences in the behavior of the two neurons?

% RS
firing_rateRS=n_of_spikesRS*10^3/stim_time;     % firing rate [Hz]
fminRS=firing_rateRS(ThRS_idx);                 % minimum frequency of spiking

% CH
firing_rateCH=n_of_spikesCH*10^3/stim_time;     % firing rate [Hz]
fminCH=firing_rateCH(ThCH_idx);                 % minimum frequency of spiking

% Comparison of the gain funtion of Regular Spiking neuron and Chattering neuron
f=f+1;
figure(f)
set(gcf,'position',[100 100 1200 500])
axes('position',[0.04 0.56 0.45 0.37])
plot(Istim,firing_rateRS)
hold on
plot(Istim,firing_rateCH)
plot(IThRS,fminRS,'+',LineWidth=1,MarkerSize=6,Color=[0, 0.4470, 0.7410])
plot(IThCH,fminCH,'x',LineWidth=1,MarkerSize=6,Color=[0.8500, 0.3250, 0.0980])
hold off
xlim([min(Istim) max(Istim)])
ylim([min([firing_rateRS;firing_rateCH]) max([firing_rateRS;firing_rateCH])])
ylabel('Firing Rate [Hz]')
xlabel('Applied Current [pA]') 
box off
legend('RS','CH',sprintf('I_t_h: %.2f pA \nf_m_i_n: %.2f Hz',IThRS,fminRS), ...
    sprintf('I_t_h: %.2f pA \nf_m_i_n: %.2f Hz',IThCH,fminCH),Location='northwest')
legend box off

%% (d)
% Evaluate again the gain function by considering a noisy stimulation (hint: add a Gaussian noise to the 
% current pulses) and discuss the achieved results compared to the ones obtained in (c)

rng(1)                                                                      % for repeatability
noise=randn(1,noff+1-non);                                                  % Gaussian noise to add to the stimulation
Istim_in_plus_noise=Istim_in+5*noise;                                       % matrix with the vectors of current of the different trials between onset and offset of the stimulation with noise

% Regular Spiking Neuron model with noisy current  
[~,Vout2PlotRS_plus_noise, SpikesCountRS_plus_noise, IstimOutRS_plus_noise] ...
    = Izhikevich_RS(dt, E_rest, Vspike, Istim_in_plus_noise,t,non,noff,Ntrials);
n_of_spikesRS_plus_noise=sum(SpikesCountRS_plus_noise,2);                   % count of the spike per trial
firing_rateRS_plus_noise=n_of_spikesRS_plus_noise*10^3/stim_time;           % firing rate [Hz]

ThRS_plus_noise_idx=find(n_of_spikesRS_plus_noise, 1 );                     % find the index of the first current(= minimum current) that make the neuron fire
IThRS_plus_noise=Istim(ThRS_plus_noise_idx);                                % current threshold for spiking
fprintf('0) %.2f \n',IThRS_plus_noise)
fminRS_plus_noise=firing_rateRS_plus_noise(ThRS_plus_noise_idx);            % minimum frequency of spiking

% Chattering Neuron model with noisy current  
[~,Vout2PlotCH_plus_noise, SpikesCountCH_plus_noise, IstimOutCH_plus_noise] ...
    = Izhikevich_CH(dt, E_rest, Vspike, Istim_in_plus_noise,t,non,noff,Ntrials);
n_of_spikesCH_plus_noise=sum(SpikesCountCH_plus_noise,2);                   % count of the spike per trial
firing_rateCH_plus_noise=n_of_spikesCH_plus_noise*10^3/stim_time;           % firing rate [Hz]

ThCH_plus_noise_idx=find(n_of_spikesCH_plus_noise, 1 );                     % find the index of the first current(= minimum current) that make the neuron fire
IThCH_plus_noise=Istim(ThCH_plus_noise_idx);                                % current threshold for spiking
fprintf('0) %.2f \n',IThCH_plus_noise)
fminCH_plus_noise=firing_rateCH_plus_noise(ThCH_plus_noise_idx);            % minimum frequency of spiking


% Comparison RS-CH gain functions both with noise
axes('position',[0.54 0.56 0.45 0.37])
plot(Istim,firing_rateRS_plus_noise)
hold on
plot(Istim,firing_rateCH_plus_noise)
plot(IThRS_plus_noise,fminRS_plus_noise,'+',LineWidth=1,MarkerSize=6,Color=[0, 0.4470, 0.7410])
plot(IThCH_plus_noise,fminCH_plus_noise,'x',LineWidth=1,MarkerSize=6,Color=[0.8500, 0.3250, 0.0980])
hold off
xlim([min(Istim) max(Istim)])
ylim([min([firing_rateRS_plus_noise;firing_rateCH_plus_noise]) max([firing_rateRS_plus_noise;firing_rateCH_plus_noise])])
ylabel('Firing Rate [Hz]')
xlabel('Applied Current [pA]') 
box off
legend('RS with noisy stimulation','CH with noisy stimulation',sprintf('I_t_h: %.2f pA \nf_m_i_n: %.2f Hz',IThRS_plus_noise,fminRS_plus_noise), ...
    sprintf('I_t_h: %.2f pA \nf_m_i_n: %.2f Hz',IThCH_plus_noise,fminCH_plus_noise),Location='northwest')
legend box off

% Comparison RS-RS with noise gain functions
axes('position',[0.04 0.09 0.45 0.37])
plot(Istim,firing_rateRS)
hold on
plot(Istim,firing_rateRS_plus_noise)
hold off
xlim([min(Istim) max(Istim)])
ylim([min([firing_rateRS;firing_rateRS_plus_noise]) max([firing_rateRS;firing_rateRS_plus_noise])])
ylabel('Firing Rate [Hz]')
xlabel('Applied Current [pA]') 
box off
legend('RS','RS with noisy stimulation',Location='northwest')
legend box off

% Comparison CH-CH with noise gain functions
axes('position',[0.54 0.09 0.45 0.37])
plot(Istim,firing_rateCH)
hold on
plot(Istim,firing_rateCH_plus_noise)
hold off
xlim([min(Istim) max(Istim)])
ylim([min([firing_rateCH;firing_rateCH_plus_noise]) max([firing_rateCH;firing_rateCH_plus_noise])])
ylabel('Firing Rate [Hz]')
xlabel('Applied Current [pA]') 
box off
legend('CH','CH with noisy stimulation',Location='northwest')
legend box off

sgtitle('Gain Functions','FontWeight','bold')