function [ sim_stats ] = simulate_network()
% Simulate neural network defined in create_network and set_parameters. 
fprintf('simulating network ... \n'); tic;
%% DECLARATIONS & INITIALISATIONS
% Declare global parameters
global par V_m I ADP osc MUA delay ref_period weight_matrix weight_matrix_STDP I_REC wc_pos wc_neg;
sympref('HeavisideAtOrigin', 1); % set heaviside (dirac) function to 1, for when t = t_spike

%% PRE-SIMULATION CALCULATIONS
% stores time since last spike for each neuron (used for ADP function)
last_spikes = zeros(par.network_size,1); 

% defines Theta STDP multiplier, 1 for trough 0 for peak
theta = 1-(osc.Hip/max(osc.Hip) + 1)/2;  

% find unique synaptic-time-constants in network and pre-calculate each EPSP
psp_u = unique(par.tau_syn); psp = cell(length(psp_u),1);
for i=1:length(psp_u)
    psp{i} = transpose(compute_psp(psp_u(i)));
end

% RECORD INPUT TO GROUPS OF NEURONS
I_REC.HIP_THETA = theta;
I_REC.NC_NC = zeros(1, par.sim_length); I_REC.NC_HIP = zeros(1, par.sim_length);
I_REC.HIP_NC = zeros(1, par.sim_length); I_REC.HIP_HIP_WIT = zeros(1, par.sim_length); 
I_REC.HIP_HIP_BET = zeros(1, par.sim_length); I_REC.ADP_HIP = zeros(1, par.sim_length);
I_REC.BG_NC = mean(I.SYN(1:par.n_NC,:)) + mean(I.OSC(1:par.n_NC,:)); 
I_REC.BG_HIP = mean(I.SYN(par.n_NC+1:end,:)) + mean(I.OSC(par.n_NC+1:end,:));
I_REC.LTP = zeros(par.network_size, par.network_size, par.sim_length);
I_REC.LTD = zeros(par.network_size, par.network_size, par.sim_length);
I_REC.HET = zeros(par.network_size, par.network_size, par.sim_length);
I_REC.STDP_decay = zeros(par.network_size, par.network_size, par.sim_length);

%% MAIN SIMULATION
for t= 2:par.sim_length
    %% AFTER-DEPOLARISAION (ADP) CURRENT
    if(par.Hip_ADP)
        I.ADP(par.n_NC+1 : par.n_NC+par.n_Hip, t) = ADP.Hip(min((t) - ...
            last_spikes(par.n_NC+1 : par.n_NC+par.n_Hip),length(ADP.Hip)));
        I_REC.ADP_HIP(:,t) = mean(ADP.Hip(min((t) - last_spikes(par.n_NC+1 : par.n_NC+par.n_Hip),length(ADP.Hip))));
    end
    
    %% VOLTAGE CHANGE
    % equation for integrate & fire voltage updates
    I.L(:,t) = par.G_m .* (par.E - V_m(:,t-1));                             % membrane decay
    I_all = I.L(:,t) + I.SYN(:,t) + I.OSC(:,t) + I.DC(:,t) + I.ADP(:,t);    % inputs
    V_m(:,t) = V_m(:,t-1) + I_all ./ par.C_m .* ref_period(:,t);            % update membrane potentials
    
    %% ADD SPIKES
    x = find(V_m(:,t) > par.V_peak); % find neurons over threshold
    MUA(x,t) = MUA(x,t) + 1; last_spikes(x) = t;  % add spike to MUA and update last spike time
    ref_period(MUA(:,t)==1, t+1:min(par.sim_length,t+par.V_ref)) = 0; % update refractory period (clamp)
    V_m(MUA(:,t)==1, t) = par.E; % reset voltage to E
    
%     % update weight matrix
    weight_matrix(:,:,t) = weight_matrix(:,:,t-1);
    x = find(MUA(:,t)>0); % find spike events
    wm = weight_matrix(x,:,t); [c, y]=find(wm>0); clear('wm'); % find active synapses at spiking neurons
    for i = 1:length(c) % add spikes to input of connected neurons
        % find relevant EPSP based on synaptic time constant of spiking neurons
        epsp = psp{find(psp_u==par.tau_syn(x(c(i)),y(i)))}*weight_matrix(x(c(i)),y(i),t);
        del = t + delay(x(c(i)),y(i)); % find delay between relevant neurons
        % add EPSP to future input of downstream neurons
        if(par.EC_modulation)
           if( x(c(i)) <= par.n_NC && y(i) > par.n_NC ) % NC -> HIP
               % no EC_mod after learning
               if(par.recall)
                   EC_mod = ones(1, length(epsp));
               else
%                    EC_mod = (theta(del:min((del-1)+length(epsp),par.sim_length)) + (1-par.EC_W)) / ( 1 + (1-par.EC_W) );
                   EC_mod = osc.EC(del:min((del-1)+length(epsp),par.sim_length));
               end
%                if(par.recall); EC_mod = 1-EC_mod; end
           else; EC_mod = ones(1, length(epsp));
           end
        else; EC_mod = ones(1, length(epsp));
        end
        I.SYN(y(i), del:min((del-1)+length(epsp),par.sim_length)) = ...
            I.SYN(y(i), del:min((del-1)+length(epsp),par.sim_length)) + ...
            epsp(1:min(length(epsp),par.sim_length-del+1)) ...
            .* EC_mod(1:min(length(epsp),par.sim_length-del+1)); 
        
        % RECORD INPUT STREAMS FOR GROUPS OF NEURONS (HIP & NC) - only necessary for Figure 3Cii
        if(x(c(i)) <= par.n_NC && y(i) <= par.n_NC) % NC -> NC
            I_REC.NC_NC(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.NC_NC(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) <= par.n_NC && y(i) > par.n_NC) % NC -> HIP
            I_REC.NC_HIP(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.NC_HIP(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) > par.n_NC && y(i) <= par.n_NC) % HIP -> NC
            I_REC.HIP_NC(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                I_REC.HIP_NC(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                epsp(1:min(length(epsp),par.sim_length-del+1));
        elseif(x(c(i)) > par.n_NC && y(i) > par.n_NC) % HIP -> HIP
            if(x(c(i)) <= par.n_NC + par.n_Hip/2 && y(i) <= par.n_NC + par.n_Hip/2 || ...
                    x(c(i)) > par.n_NC + par.n_Hip/2 && y(i) > par.n_NC + par.n_Hip/2) % WITHIN GROUPS
                I_REC.HIP_HIP_WIT(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                    I_REC.HIP_HIP_WIT(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                    epsp(1:min(length(epsp),par.sim_length-del+1));
            elseif(x(c(i)) <= par.n_NC + par.n_Hip/2 && y(i) > par.n_NC + par.n_Hip/2 || ...
                    x(c(i)) > par.n_NC + par.n_Hip/2 && y(i) <= par.n_NC + par.n_Hip/2) % BETWEEN GROUPS
                I_REC.HIP_HIP_BET(:, del:min((del-1)+length(epsp),par.sim_length)) = ...
                    I_REC.HIP_HIP_BET(:, del:min((del-1)+length(epsp),par.sim_length)) + ...
                    epsp(1:min(length(epsp),par.sim_length-del+1));
            end
        end
        clear('epsp');
    end
    
    %% STDP LEARNING
    spikes = find(MUA(:,t)==1);
    STDP( spikes, t, theta(t) );
    fprintf('\b\b\b\b\b%3.0f%%\n',t/par.sim_length*100);
end


%% SAVE SIMULATION DATA
% SPIKE DETECTOR
[s, d] = find(MUA>0); 
sim_stats.spike_detector = [reshape(s,[length(s),1]) reshape(d,[length(d),1])];
ID = transpose(1:par.network_size); 
% VOLTMETER
sim_stats.voltmeter = zeros(par.sim_length*par.network_size, 3);
for t=1:par.sim_length
   sim_stats.voltmeter((t-1)*length(ID)+1:t*length(ID),:) = [ID t*ones(par.network_size,1) V_m(:,t)];
end
% INPUT RECORDING STREAMS
I_REC.NC_NC = I_REC.NC_NC / par.n_NC; I_REC.HIP_NC = I_REC.HIP_NC / par.n_NC; 
I_REC.HIP_HIP_WIT = I_REC.HIP_HIP_WIT / par.n_Hip; I_REC.HIP_HIP_BET = I_REC.HIP_HIP_BET / par.n_Hip; 
I_REC.NC_HIP = I_REC.NC_HIP / par.n_Hip; 
I_REC.wc_pos = wc_pos; I_REC.wc_neg = wc_neg; 
I_REC.I = I;
sim_stats.I_REC = I_REC; clear('I_REC');
% REFRACTORY PERIOD
sim_stats.RP = ref_period; clear('ref_period');
% WEIGHT MATRIX AND ACTIVE SYNAPSES
sim_stats.weight_matrix = weight_matrix; sim_stats.weight_matrix_STDP = weight_matrix_STDP;

fprintf('\b\b\b\b\b finished in %.0f seconds\n', toc);
end

