function [ ] = set_parameters( )
% Create set of parameters used for simulating network.

%% SET PARAMETER SPACE FUNCTION
    global par;
     if(isempty(par)==1)
        par = struct;
     end

    %% INDEPENDENT PARAMETERS
    parameter_settings = ...
        ... % NEURONS PROPERTIES
         {'V_peak', -55, 'E', -70, 'V_ref', 2,'V_m_initial', -65, ... % intracellular
         'C_m', .9, 'G_m', .03, 'tau_m', 20, ... % intracellular
        'Hip_SR', 1500, 'Hip_SW', .015, 'Hip_ST', 1.5,... % Hip spike rate/weight/tau
        'NC_SR', 4000, 'NC_SW', .023, 'NC_ST', 1.5,... % NC spike rate/weight/tau
        'NC_amp', .1, 'NC_freq', 10, 'NC_reset', false, 'NC_ADP', false, ... % Neo-Cortex Alpha
        'Hip_amp', .25, 'Hip_freq', 4, 'Hip_reset', true, 'reset_ph', 90,...% Hippocampal Theta
        'Hip_ADP', true, 'Hip_ADP_amp', .2, 'Hip_ADP_freq', 4,... % Hippocampal ADP
        'rand_phase', true, ... % random phase for oscillation cosine waves ('on' or 'off')
        ... % NEURON GROUPS PROPERTIES
        'n_Items', 2, 'NC_per_item', 10, 'Hip_per_item', 5,... % group properties
        'EC_modulation', true, 'EC_W', 0.3, 'rand_EC_phase', false,...
        'NC_Hip_intra_syn_str', .35, 'NC_Hip_intra_R_conn', 100, 'NC_Hip_intra_delay', 2, 'NC_Hip_intra_tau_syn', 1.5,... % NC -> Hip connectivity
        'NC_Hip_inter_syn_str', 0, 'NC_Hip_inter_R_conn', 0, 'NC_Hip_inter_delay', 2, 'NC_Hip_inter_tau_syn', 1.5,... % NC -> Hip connectivity
        'Hip_NC_intra_syn_str', .08, 'Hip_NC_intra_R_conn', 100, 'Hip_NC_intra_delay', 2, 'Hip_NC_intra_tau_syn', 1.5,... % Hip -> NC connectivity
        'Hip_NC_inter_syn_str', 0, 'Hip_NC_inter_R_conn', 0, 'Hip_NC_inter_delay', 2, 'Hip_NC_inter_tau_syn', 1.5,... % Hip -> NC connectivity
        'NC_intra_syn_str', .3, 'NC_intra_R_conn', 25, 'NC_intra_delay', 2, 'NC_intra_tau_syn', 1.5,... % NC connectivity
        'NC_inter_syn_str', 0, 'NC_inter_R_conn', 0, 'NC_inter_delay', 2, 'NC_inter_tau_syn', 1.5,... % NC connectivity
        'Hip_intra_syn_str', .65, 'Hip_intra_R_conn', 50, 'Hip_intra_delay', 2, 'Hip_intra_tau_syn', 1.5,... % Hip connectivity
        'Hip_inter_syn_str', 0, 'Hip_inter_R_conn', 50, 'Hip_inter_delay', 2, 'Hip_inter_tau_syn', 1.5,... % Hip connectivity
        'W_m', .5, 'W_sd', .05, 'EC_SD', 0.5, ... 
        ... % STDP PROPERTIES
        'NC_STDP', false, 'Hip_STDP', true, ... 
        'theta_LTD', false, 'theta_mod', true, 'STDP_th', true, 'STDP_hetero', false ...
        'weight_decay', false, 'T_weight_decay', 0.1667, ...
        'weight_max', .65, 'STDP_tau', 5, ...
        'T_Ca', 20, 'T_WC', 10, 'a_pos', .65, 'a_neg', .65, ...
        'T_p', 1, 'T_d', 1, 'T_h', .8, 'G_p', 1.5, 'G_d', .75, 'G_h', .55, 'P_h', 1, ...
        ... % STIMULUS PRESENTATION PROPERTIES
        'idling_length', 1000, 'pre_stim_length', 2000, 'post_stim_length', 0, 'stim_length', 3000,... % length of simulation periods (ms)
        'stim_type', 'DC', ... % set input type: DC or alpha
        'stim_TS', 400, 'stimulus_SR', 6000, 'stimulus_SW', .1, 'stimulus_ST', 1.5, ... % for alpha modulated spikes
        'stim_mod', true, 'stim_range', false, 'stim_F', 4, 'stim_G1_ph', 270, 'stim_G2_ph', 180, ... % for oscillatory modulation
        'stim_n', 0, 'burst_n', [], 'stimulus_strength', 1.75}; % for strength & timing DC input (or burst)
   
    % SET INDEPENDENT par IF NOT ALREADY SUPPLIED
    for i=1:length(parameter_settings)/2
        %if(any(strcmp(parameter_settings{(i-1)*2+1},fieldnames(par)))==0)
            par.([parameter_settings{(i-1)*2+1}]) = parameter_settings{(i-1)*2+2};
        %end
    end
    
    %% DEPENDENT PARAMETERS   
    % GROUP SIZES & IDs
    par.n_Hip = par.Hip_per_item * par.n_Items;
    par.n_NC = par.NC_per_item * par.n_Items;
    par.network_size = par.n_NC + par.n_Hip;
    par.nG = {'NC', 'Hip'};
    
    % ORDER OF THE SIMULATION
    par.sim_order = {'idling', 'learning', 'NP presentation', 'P presentation', 'idling'}; 
    par.sim_order_n = {'idle-BL', 'DL', 'NP-AL', 'P-AL', 'idle-AL'}; % used for filenames
%     par.sim_order = {'learning'};
%     par.sim_order_n = {'DL'};

end



