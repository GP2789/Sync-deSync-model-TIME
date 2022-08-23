function [  ] = make_figures( filename, p, freqs, phase_con, no_flicker )
% filename: string
% p: a modifier for variable sampling & number of trials (0<p<=1).
% freqs: stimulus modulation frequency (has to be this order to reproduce the same figures): 1.652, 4, 10.472 to plot figure 4
% phase_con: phase for 'stim_G2_ph' (has to be this order to reproduce the same figures): 270,0,90,180 (stim_G1_ph = 270)
% if no_flicker is true, stimulus is not frequency modulated: for figure 4B

    %% Full model: Weight change time series, Firing activity time series, Mean WC as a function of phase offset condition 
    % pure input frequency
    trials = ceil(p * 200);
    % Pure hippocampal frequency: 4 Hz
    hip_rf = ones(1,trials)*4;  
    
    analyse_filenames = cell(length(freqs),length(phase_con));
    % loop through each frequency and phase offset condition
    for nf = 1:length(freqs)
        rf = ones(1,trials)*freqs(nf);
        stim_str = 1.75*exp((freqs(nf)/20).^3); 
        for ph = 1:length(phase_con)
            filenametemp = [filename '_freq' num2str(freqs(nf)) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            vars = evaluate_network({'stimulus_strength', {stim_str}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
            analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 3A-Ci 
            analyse_filenames{nf,ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    end
    analyse_results(analyse_filenames, filename, {'FR','WC'}, trials, phase_con, freqs, false); % generate figure 3A 3B and Ci
    
    % for plotting figure 4 and figure s3
    single_trial_raster(analyse_filenames, filename, {'ST'}, trials, phase_con, freqs);
    
    %% Full model: Mean WC as a function of phase offset condition 
    % Adding noise to input frequency
    % Pure hippocampal frequency: 4 Hz
    hip_rf = ones(1,trials)*4;  
    
    % for figure 3Cii: freqs = 4 Hz std = 0.015*4
    sp = 0.015; %define the STD of the input frequencies distribution
    rf = randsample(randn(1,1000)*sp*freqs(2) + freqs(2), trials);
    analyse_filenames = cell(1,length(phase_con));
    % loop through each phase offset condition
    for ph = 1:length(phase_con)
        filenametemp = [filename '_rand_input_freq' num2str(freqs(2)) '_std' num2str(sp) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
        mkdir(filenametemp);

        vars = evaluate_network({'Default'}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
        analyse_network(vars.test_f, {'FR', 'WC'}); 
        analyse_filenames{ph} = vars.test_f; 
        clear filenametemp
        clearvars -global;
    end
    
    analyse_results(analyse_filenames, [filename '_rand_input_freq' num2str(freqs(2)) '_std' num2str(sp)], {'WC'}, trials, phase_con, freqs(2), false); 

    %% Full model: Weight change time series, Firing activity time series, Mean WC as a function of phase offset condition 
    % for figure 3Ciii & 5Aiii: Adding noise to Hippocampal frequency: mean = 4 std = 0.02
    hip_std = 0.02;
    hip_rf = randsample(randn(1,1000)*hip_std + 4, trials); 
    
    analyse_filenames = cell(length(freqs),length(phase_con));    
    % loop through each phase offset condition and phase offset condition
    for nf = 1:length(freqs)
        % pure input frequency: freqs = 1.652, 4 or 10.472 Hz
        rf = ones(1,trials)*freqs(nf);
        stim_str = 1.75*exp((freqs(nf)/20).^3); 
        for ph = 1:length(phase_con)
            filenametemp = [filename '_freq' num2str(freqs(nf)) '_rand_hipp_freq4_std' num2str(hip_std) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            % adding noise to EC phase offset distribution: mean = 180 degrees
            % std = 0.5/3
            vars = evaluate_network({'rand_EC_phase', {true},'stimulus_strength', {stim_str}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); 
            analyse_network(vars.test_f, {'FR', 'WC'});  
            analyse_filenames{nf,ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    end
    
    analyse_results(analyse_filenames, [filename '_rand_hipp_freq4_std' num2str(hip_std)], {'WC'}, trials, phase_con, freqs, false); 

    %% Full model: comparing with a non-flickered stimulus 
    % pure input frequency
    % Pure hippocampal frequency: 4 Hz
    hip_rf = ones(1,trials)*4;  
    % only compare 4 Hz 0 vs. 180 (actual phase offset for stim_G2_ph: 270
    % vs. 90) vs. non-flickered stimulus
    if no_flicker
        freqs = 4;
        phase_con = [270 90];
        analyse_filenames = cell(1,3);
        % loop through each phase offset condition
        rf = ones(1,trials)*freqs;
        for ph = 1:length(phase_con)
            filenametemp = [filename '_freq' num2str(freqs) '_vs_noflicker_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            vars = evaluate_network({'Default'}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
            analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 4Bii 
            analyse_filenames{ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    
        % run analyses for non-flickered stimulus
        filenametemp = [filename '_noflicker_' num2str(p)];
        mkdir(filenametemp);
        vars = evaluate_network({'stim_length', {1500},'stim_mod', {false}},trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
        analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 4Bii
        analyse_filenames{3} = vars.test_f; 
        clear filenametemp
        clearvars -global;
        
        analyse_results(analyse_filenames, filename, {'WC'}, trials, phase_con, freqs, true); % generate figure 4Bii
    end

    %% Full model: comparing with a non-flickered
    % pure input frequency
    % for figure 5Biii: Adding noise to Hippocampal frequency: mean = 4 std = 0.02
    hip_std = 0.02;
    hip_rf = randsample(randn(1,1000)*hip_std + 4, trials); 
    % only compare 4 Hz 0 vs. 180 (actual phase offset for stim_G2_ph: 270
    % vs. 90) vs. non-flickered stimulus
    if no_flicker
        freqs = 4;
        phase_con = [270 90];
        analyse_filenames = cell(1,3);
        % loop through each phase offset condition
        rf = ones(1,trials)*freqs;
        for ph = 1:length(phase_con)
            filenametemp = [filename '_rand_hipp_freq4_std' num2str(hip_std) '_freq' num2str(freqs) '_vs_noflicker_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            % adding noise to EC phase offset distribution: mean = 180 degrees
            % std = 0.5/3
            vars = evaluate_network({'rand_EC_phase', {true}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Biii 
            analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 4Biii 
            analyse_filenames{ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    
        % run analyses for non-flickered stimulus
        filenametemp = [filename '_rand_hipp_freq4_std' num2str(hip_std) '_noflicker_' num2str(p)];
        mkdir(filenametemp);
        % adding noise to EC phase offset distribution: mean = 180 degrees
        % std = 0.5/3
        vars = evaluate_network({'rand_EC_phase', {true},'stim_length', {1500},'stim_mod', {false}},trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Biii 
        analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 4Biii
        analyse_filenames{3} = vars.test_f; 
        clear filenametemp
        clearvars -global;
        
        analyse_results(analyse_filenames, [filename '_rand_hipp_freq4_std' num2str(hip_std)], {'WC'}, trials, phase_con, freqs, true); % generate figure 4Biii
    end

    %% Full model: Mean WC as a function of phase offset condition for wider frequency and phase offset range
    % for figure 6: freqs order: [1.652 4 10.472 18.335 41.236 71.771] 
    % phase conditions: [45 90 135 180 225 270 315 0] 
    % pure input frequency
    % Pure hippocampal frequency: 4 Hz
    trials = ceil(p * 200);
    hip_rf = ones(1,trials)*4;  
    
    analyse_filenames = cell(length(freqs),length(phase_con));
    % loop through each frequency and phase offset condition: delta theta
    % alpha
    for nf = 1:length(freqs)-3
        rf = ones(1,trials)*freqs(nf);
        stim_str = 1.75*exp((freqs(nf)/20).^3); 
        for ph = 1:length(phase_con)
            filenametemp = [filename '_freq' num2str(freqs(nf)) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            vars = evaluate_network({'stimulus_strength', {stim_str}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
            analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 3A-Ci 
            analyse_filenames{nf,ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    end
    
    for nf = length(freqs)-2:length(freqs)
        rf = ones(1,trials)*freqs(nf);
        stim_str = 2.2*log10(freqs(nf)); 
        for ph = 1:length(phase_con)
            filenametemp = [filename '_freq' num2str(freqs(nf)) '_phase' num2str(phase_con(ph)) '_' num2str(p)];
            mkdir(filenametemp);

            vars = evaluate_network({'stimulus_strength', {stim_str}}, trials, filenametemp, hip_rf, rf,  phase_con(ph)); % simulate Figures 4Bii 
            analyse_network(vars.test_f, {'FR', 'WC'}); % analyse data for Figures 3A-Ci 
            analyse_filenames{nf,ph} = vars.test_f; 
            clear filenametemp
            clearvars -global;
        end
    end
    
    analyse_results(analyse_filenames, filename, {'WC'}, trials, phase_con, freqs, false); % generate figure 3A 3B and Ci
