function [ ] = single_trial_raster( filename, savefilename, type, trials, phase_con, freqs)
freq.label = {'delta','theta','alpha'};
freq.value = [1.652 4 10.472];
con_label = {['0',char(176), 'phase offset'],  ['90',char(176), ' phase offset'],['180',char(176), ' phase offset'], ['270',char(176), ' phase offset']};
for f = 1:length(freqs)
    for con = 1:length(phase_con)
        %% INITIALISE DATA STRUCTURES
        load([filename{f,con} '/variables.mat']); par = vars.par;
        v1_L = max(1,length(vars.v1_vars)); 
        v2_L = max(1,length(vars.v2_vars));
        v3_L = max(1,length(vars.v3_vars));
    
        %% AVERAGE FIRING RATE PLOTS OVER TRIALS
        if(sum(cellfun(@sum,strfind(type,'ST')))>0) %Tyep = Single trial
            % loop through filenames
            for k=1:v3_L % loop variable 1
                if(isempty(vars.v3_vars)~=1)
                    if(ischar(vars.v3_vars{k})~=1); v3 = num2str(vars.v3_vars{k}); 
                    else; v3 = vars.v3_vars{k}; end
                else; v3 = []; 
                end; filename_1 = [vars.test_f '/' vars.v3_name v3];
                for i=1:v1_L % loop variable 2
                    if(isempty(vars.v1_vars)~=1)
                        if(ischar(vars.v1_vars{i})~=1); v1 = num2str(vars.v1_vars{i}); 
                        else; v1 = vars.v1_vars{i}; end
                    else; v1 = []; 
                    end; filename_2 = [filename_1 '/' vars.v1_name v1];
                    for j=1:v2_L % loop variable 3
                        if(isempty(vars.v2_vars)~=1)
                            if(ischar(vars.v2_vars{j})~=1); v2 = num2str(vars.v2_vars{j}); 
                            else; v2 = vars.v2_vars{j}; end
                            filename_3 = [filename_2 '_' vars.v2_name v2];
                        else; filename_3 = filename_2;
                        end
                        
                                            
                        %% EXTRACT DATA
                        t=trials;
                        x1 = 0: 1/1000 : par{k}{i,j}{t}.stim_length/1000-1/1000;

                        r_trial = ceil(rand()*trials); % choose random trial
                        % load single trial data
                        load([vars.test_f '/Simulations/' vars.v1_name num2str(vars.par{k}{i,j}{t}.stimulus_strength)...
                            '/T' int2str(r_trial) '/Data/DL.mat']);
                        
                        r_P_randseq = [];
                        while isempty(r_P_randseq)
                            hip_NP = ceil(rand() * par{k}{i,j}{t}.n_Hip/2);
                            r_NP = par{k}{i,j}{t}.n_NC + hip_NP + par{k}{i,j}{t}.n_Hip/2; % choose random NP neuron
                            % choose a random P that connects with NP
                            % load DL data for trial
                            r_P_seq = find(data.sim_stats.weight_matrix_STDP(r_NP,par{k}{i,j}{t}.n_NC+1:par{k}{i,j}{t}.n_NC+5)'...
                                &data.sim_stats.weight_matrix_STDP(par{k}{i,j}{t}.n_NC+1:par{k}{i,j}{t}.n_NC+5,r_NP))+par{k}{i,j}{t}.n_NC;
                            r_P_randseq = r_P_seq(randperm(length(r_P_seq)));
                        end
                        r_P = r_P_randseq(1);

                        clear data;
                        
                        % load DL data for trial
                        load([vars.test_f '/Simulations/' vars.v1_name num2str(vars.par{k}{i,j}{t}.stimulus_strength) '/T' int2str(r_trial) '/Data/DL.mat']); 
                        spikes = data.sim_stats.spike_detector; spikes(:,2) = spikes(:,2) - par{k}{i,j}{t}.pre_stim_length;
                        P_DL = spikes(spikes(:,1)==r_P,:); NP_DL = spikes(spikes(:,1)==r_NP,:);
                        P_DL(:,1) = P_DL(:,1) - r_P + 2; NP_DL(:,1) = NP_DL(:,1) - r_NP + 1;

                        % weights between groups
                        Hip_all = par{k}{i,j}{t}.n_NC + 1 : par{k}{i,j}{t}.network_size;
                        Hip_G1 = Hip_all(1 : par{k}{i,j}{t}.Hip_per_item);
                        Hip_G2 = Hip_all(par{k}{i,j}{t}.Hip_per_item + 1 : end);

                        n_inter = [sum(data.sim_stats.weight_matrix_STDP(Hip_G1, Hip_G2),'all') ...
                        sum(data.sim_stats.weight_matrix_STDP(Hip_G2, Hip_G1),'all')];
                        WM_inter_NP = (squeeze(data.sim_stats.weight_matrix( r_NP, r_P, : )))'/ par{k}{i,j}{t}.weight_max;
                        WM_inter_P = (squeeze(data.sim_stats.weight_matrix( r_P, r_NP, : )))'/ par{k}{i,j}{t}.weight_max;
                        STDP{con}.wm_STDP = data.sim_stats.weight_matrix_STDP;
                        clear('data');
                        % plot hipp theta 
                        theta_vis = cos(2*pi*par{k}{i,j}{t}.Hip_freq*(x1+1/par{k}{i,j}{t}.Hip_freq*par{k}{i,j}{t}.reset_ph/360))*1/2+1/2+7;
                        theta_aud = cos(2*pi*par{k}{i,j}{t}.Hip_freq*(x1+1/par{k}{i,j}{t}.Hip_freq*par{k}{i,j}{t}.reset_ph/360))*1/2+1/2+5;
                        mod_P = cos(2*pi*freqs(f)*(x1+1/freqs(f)*par{k}{i,j}{t}.stim_G1_ph/360))*1/2+1/2+7;
                        mod_NP = cos(2*pi*freqs(f)*(x1+1/freqs(f)*par{k}{i,j}{t}.stim_G2_ph/360))*1/2+1/2+5;
                        
                        mean_WM(con,1) = mean(WM_inter_NP(2751:end));
                        if con == 1
                           fst(f) = figure;
                           set(fst(f),'position',[0        0        2000       1000]);
                        end
                        T = -2:1/1000:3-1/1000; %(1:1:sum(par{k}{i,j}{t}.sim_length_total(2)))/1000; %DL sim length
                        figure(fst(f));
                        subplot(4,2,(con-1)*2+1); hold on;  
                        plot(T, WM_inter_NP,'k-','linewidth',1); 
                        
                        plot(T, WM_inter_P,'k--','linewidth',2); 
                        ylabel('norm. weights'); ylim([0 1]);  
                        xlim([-1 3]); 
                        xticklabels({'-1','0','1','2','3'});
                        if con == 1
                           title([{'Stimulus Frequency: ' freq.label{freq.value...
                            == freqs(f)}},{'Hipp. synaptic change between groups ', con_label{con}}]);
                           legend('Auditory --> Visual neuron','Visual --> Auditory neuron','EdgeColor','none',...
                               'Color','none','Location',[0.42 0.93 0.15 0.07]);
                        else
                            title(['Hipp. synaptic change between groups ', con_label{con}]);
                        end
                        xlabel('Time (s)');
                        ax = gca; ax.FontSize = 12;
                        ax.FontName = 'Arial';

                        subplot(4,2,(con-1)*2+2); hold on; ylim([0.5 8.5]);
                        xticks(-1000:1000:3000); xticklabels({'-1','0','1','2','3'});
                        for s = 1:length(P_DL); plot([P_DL(s,2) P_DL(s,2)], [3 4],'color',[0 0.4 0.6]); end
                        for s = 1:length(NP_DL); plot([NP_DL(s,2) NP_DL(s,2)], [1 2],'color',[0.8 0.2 0]); end
                        r(1) = plot(x1*1000,mod_NP,'linewidth',2,'color',[0.8 0.4 0]);
                        r(2) = plot(x1*1000,mod_P,'linewidth',2,'color',[0 0.45 0.7]);
                        r(3) = plot(x1*1000,theta_aud,'color',[0.5 0.5 0.5]);
                        plot(x1*1000,theta_vis,'color',[0.5 0.5 0.5]);
                        if con == 1
                           title([{'Stimulus Frequency: ' freq.label{freq.value...
                            == freqs(f)}},{'Hipp. neuronal activity during learning ', con_label{con}}]); 
                           legend(r,{'Auditory neuron','Visual neuron','Hipp Theta'},'EdgeColor','none','Color',...
                               'none','Location',[0.87 0.90 0.11 0.09]);
                        else
                           title(['Hipp. neuronal activity during learning ', con_label{con}]);                     
                        end
                        yticks([1.5 3.5 5.5 7.5]); yticklabels({'Aud','Vis', 'Aud', 'Vis'}); xlim([-1000 3000]); 
                        xlabel('Time (s)'); 
                        ax = gca; ax.FontSize = 12; 
                        ax.FontName = 'Arial';
                        
                        if con == length(phase_con)
                           saveas(fst(f),['analyse_figures/' savefilename '/single_trial_raster_', freq.label{freq.value == freqs(f)},'.tiff']); close(fst(f));
                        end
                    end
                end
            end
        end
    end
end


                        

