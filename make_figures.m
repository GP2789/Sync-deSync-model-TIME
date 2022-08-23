function [ fstat_clouter, fstat_wang ] = make_figures( filename, p, freqs, phase_con )
% STDP only model 
% in set_parameters.m: EC_modulation = false theta_LTD = false 
% adding par.theta_mod: true: theta phase modulation between -1 (peak) and
% 1 (trough) false: theta = 1 (no theta phase modulation)
% Here in the STDP only model, theta_mod = false
% adding par.stim_range: true: stimulus input at 4 Hz ranging between -1 and 1 
% false: ranging between 0 1
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
        filenametemp = [filename '_STDP_only_freq' num2str(freqs) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
        mkdir(filenametemp);
        % for figure 5Aii and 5Bii define 'stim_range' to be true. For
        % figure S2, define 'stim_range' to be false
        clear vars
        vars = evaluate_network({'Hip_amp', {0},'stim_range', {false}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 5Aii and 5Bii 
        analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 5Aii and 5Bii 
        analyse_filenames{2,ph} = vars.test_f; 
        clear filenametemp
        clearvars -global;
    end
    [ fstat_clouter, fstat_wang ] = analyse_results(analyse_filenames, [filename '_figS2_stim_range0_1'], {'WC'}, trials, phase_con ); % generate figure 5Aii and 5Bii
    