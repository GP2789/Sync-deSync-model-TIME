function [ fstat_clouter, fstat_wang ] = make_figures( filename, p, freqs, phase_con )
% theta only model 
% in set_parameters.m: theta_LTD = false 
% adding par.theta_mod: true: theta phase modulation between -1 (peak) and
% 1 (trough) false: theta = 1 (no theta phase modulation)
% Here in the Theta-phase-learning only model, theta_mod = true
% adding par.stim_range: true: stimulus input at 4 Hz ranging between -1 and 1 
% false: ranging between 0 1
% Changes of inside functions:
% create_network.m: start with 0 synaptic connectivity in the hippocampus:
% line 140 str = zeros(size(syn));
% STDP.m: weights only depend on theta phases: constant = 0.02
% lines 85 86 and 93 94: dW = weight_matrix_STDP(:, spikes) * theta * 0.02; p_w(:, spikes) = p_w(:, spikes) + dW; for both directions
% To enable Theta-phase-learning dependent weight changes only, STDP_th = false
% filename: string
% p: a modifier for variable sampling & number of trials (0<p<=1).
% freqs: stimulus modulation frequency: 4
% phase_con: phase for 'stim_G2_ph' (has to be this order to reproduce the same figures): 270,0,90,180 (stim_G1_ph = 270)

    %% for figure 5
    % pure input frequency
    trials = ceil(p * 200);
    % Pure hippocampal frequency: 4 Hz
    hip_rf = ones(1,trials)*4;  
    
    analyse_filenames = cell(2,length(phase_con));
    % loop through phase offset conditions
    rf = ones(1,trials)*freqs;
    for ph = 1:length(phase_con)
        % load data from full model
        load(['../full_model/' filename '_freq4_phase' num2str(phase_con(ph)) '_' num2str(p) '/DEFAULT_' num2str(trials) 'T/variables.mat']);
        analyse_filenames{1,ph} = ['../full_model/' vars.test_f];
        filenametemp = [filename '_Theta_only_freq' num2str(freqs) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
        mkdir(filenametemp);
        clear vars
        vars = evaluate_network({'STDP_th', {false}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 5Ai and 5Bi 
        analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 5Ai and 5Bi 
        analyse_filenames{2,ph} = vars.test_f; 
        clear filenametemp
        clearvars -global;
    end
    [ fstat_clouter, fstat_wang ] = analyse_results(analyse_filenames, filename, {'WC'}, trials, phase_con ); % generate figure 5Ai and 5Bi
    