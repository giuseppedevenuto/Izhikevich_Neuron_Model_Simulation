function [Vout, Vout2Plot, SpikesCount, IstimOut] = Izhikevich_CH(dt, E_rest, Vspike, Istim_in,t,non,noff,Ntrials)
% Function for implementing the CH neuron with Izhikevich equations
% © De Venuto Giuseppe, Impollonia Roberta 2023

% Izhikevich parameters for Chattering
a=0.02;
b=0.2;
c=-50;                                                              % after-spike reset value of the mem. pot. V [mV]
d=2;                                                                % after-spike reset value of the recovery variable U

Vout=[];                                                            % Matrix with the raw mem. pot. of all the trials
Vout2Plot = [];                                                     % Matrix with the mem. pot. to be plotted of all the trials
SpikesCount = [];                                                   % Matrix with the spike trains of all the trials
IstimOut = [];                                                      % Matrix with the current applied in all the trials

% Loop through trials with different applied currents
for trial = 1:Ntrials

    I = zeros(size(t));                                             % Initialize current to zero
    I(non:noff) = Istim_in(trial,:);                                % Then generate the step current
    V = E_rest*ones(size(t));                                       % Initialize membrane potential V
    U = b*V;                                                        % Inizialize membrane recovery variable U
    V2plot=V;                                                       % Inizialize membrane potential to be plotted (spikes with the same peaks amplitude)
    spikes = zeros(size(t));                                        % Initialize vector of recorded spikes(spike train)

    % Loop through the simulation time
    for i = 2:length(t)

        V(i)= dt*(0.04*V(i-1)^2+5*V(i-1)+140-U(i-1)+I(i))+V(i-1);   % Forward Euler method for membrane potential V
        U(i)= dt*(a*(b*V(i-1)-U(i-1)))+U(i-1);                      % Forward Euler method for membrane recovery variable U

        if V(i) > Vspike                                            % If membrane potential reaches its apex
            spikes(i) = 1;                                          % Record a spike
            V2plot(i)=Vspike;                                       % Impose last mem. pot. before the reset equal to 30mV
            V(i) = c;                                               % Reset membrane potential value after the spike
            U(i) = U(i) + d;                                        % Reset recovery variable value after the spike
        else
            V2plot(i)=V(i);                                         % Mem. pot. to plot
        end
    end

    Vout = [Vout; V];                                               % Matrix with the raw mem. pot. of all the trials
    Vout2Plot = [Vout2Plot; V2plot];                                % Matrix with the mem. pot. to be plotted of all the trials
    SpikesCount = [SpikesCount; spikes];                            % Matrix with the spike trains of all the trials
    IstimOut = [IstimOut; I];                                       % Matrix with the current applied in all the trials
end