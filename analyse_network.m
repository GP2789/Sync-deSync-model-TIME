function [ ] = analyse_network( filename, type )
% Analyse simulation data in the directory specified by filename. Analyse
% by type ('FR', 'WC' for Activity, Weight Change, respectively).

%% INITIALISE DATA STRUCTURES
load([filename '/variables.mat']); par = vars.par;
v1_L = max(1,length(vars.v1_vars)); 
v2_L = max(1,length(vars.v2_vars));
v3_L = max(1,length(vars.v3_vars));
N = length(par{1}{1,1}{1}.sim_order_n);
n_sim = cell(1,N); n_group = {'NC', 'Hip'};
n_syn_g = {'P','NP','NP_P','P_NP'}; n_stage = {'DL','AL'};
n_sub_g = {'NC_P','NC_NP','Hip_P','Hip_NP','NC_All','Hip_All'}; 
n_BS = 1000;

for n=1:N
    n_sim{n} = strrep(par{1}{1,1}{1}.sim_order_n{n},'-','_');
end
n_sim{N+1} = 'All';

%% LOAD, CALCULATE & STORE SPIKE & WEIGHT DATA
% check to see if function has already been performed
% if(exist([vars.test_f '/Analysis/firing-rate-data.mat'],'file')~=2 || exist([vars.test_f '/Analysis/weight-change-data.mat'],'file')~=2) 
    h1 = waitbar(0, 'Parameter Progression', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);
    for t=1:vars.trials % loop trials
       % LOOP THROUGH FILENAMES
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
                    filename = [filename_3 '/T' int2str(t)];
                    waitbar(((t-1)*v3_L*v1_L*v2_L + (k-1)*v1_L*v2_L + ...
                        (i-1)*v2_L + j) / (vars.trials * v3_L * v1_L * v2_L),h1); % increment progress bar

                    NCPI = par{k}{i,j}{t}.NC_per_item; HPI = par{k}{i,j}{t}.Hip_per_item; % extract network properties
                    for n=1:N % loop through stages of simulation
                         % for the simulatiton that turns off STDP after
                         % learning load STDP matrix from DL
                         load([filename '/Data/DL.mat']) % load data
                         syn = data.sim_stats.weight_matrix_STDP; 
                         clear data
                         load([filename '/Data/' par{k}{i,j}{t}.sim_order_n{n} '.mat']) % load data
                         % extract spikes/weight/synapse data
                         spikes = data.sim_stats.spike_detector; weights_t = data.sim_stats.weight_matrix; 
%                          syn = data.sim_stats.weight_matrix_STDP; 
                         sim_length = unique(data.sim_stats.voltmeter(:,2));
                         
                         %% EXTRACT INPUT RECORDING DATA
                         rec = fieldnames(data.sim_stats.I_REC);
                         rec = rec(1:9);
                         for r = 1:length(rec)
                             if(t==1); FR.I.([n_sim{n}]).([rec{r}]){k}{i,j} = zeros(1, length(sim_length)); end % initialise data
                             FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  = FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  + data.sim_stats.I_REC.([rec{r}]); % increment data
                             if(t==vars.trials); FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  = (FR.I.([n_sim{n}]).([rec{r}]){k}{i,j}  / vars.trials) / par{k}{i,j}{t}.C_m; end% average data
                         end
                         clear('rec');
                         
                         %% EXTRACT ALL HIPPOCAMPAL RELATED DATA
                         for p = 1:2 
                             %% EXTRACT WEIGHT CHANGE OF HIP GROUPS
                             for w=1:2 % loop through weight groups
                                 if(w==p); v=0; else; v=1; end % extract all weight changes in group
                                 S = syn(NCPI*2 + HPI*v+1:NCPI*2 + HPI*(v+1) , NCPI*2 + HPI*(p-1)+1:NCPI*2 + HPI*p);
                                 wc = weights_t(NCPI*2 + HPI*v+1:NCPI*2 + HPI*(v+1) , NCPI*2 + HPI*(p-1)+1:NCPI*2 + HPI*p,:);
                                 [wy, wx] = find(S==1); % find non-zeros/active synapses in simulation
                                 W = zeros(length(sim_length),length(wx));
                                 for l=1:length(wx)
                                     W(:,l) = squeeze(wc(wy(l),wx(l),:));
                                 end % add to temporary matrix 
                                 if(t == 1)
                                     WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = [];
                                 end
                                 WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = ...
                                     [WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} mean(W,2)]; % concatenate matrix
                                 if(t == vars.trials)
                                     % save WC of all trials
                                     WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2},'_raw']){k}{i,j} = WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j};
                                     WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j} = mean(WC.([n_sim{n}]).([n_syn_g{p+(w-1)*2}]){k}{i,j},2);
                                 end
                             end
                             
                             %% EXTRACT Y-DATA FOR HIPPOCAMPAL ACTIVITY PLOTS
                             % extract P/NP spikes
                             FR.ACT.bin_width = 20;
                             hip_spikes = spikes(spikes(:,1) > NCPI*2 + HPI *(p-1),:);
                             hip_spikes = hip_spikes(hip_spikes(:,1) <= HPI*p + NCPI*2,:);
                             if(t==1);FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S = zeros(vars.trials,round(length(sim_length)/FR.ACT.bin_width));end
                             FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S(t, :) = histcounts(hip_spikes(:,2),length(sim_length)/FR.ACT.bin_width,'binlimits',[0 length(sim_length)]) ...
                                 / (par{k}{i,j}{t}.n_Hip/2) * (1000 / FR.ACT.bin_width);
                             if(t==vars.trials)
                                 if(sum(cellfun(@sum,strfind(type,'FR')))>0)
                                 [FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.CI_L, FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.CI_U, ...
                                     FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.M] = ...
                                     bootstrap(FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}.S, n_BS, []);
                                 end
%                                  FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j} = ...
%                                      rmfield(FR.ACT.([n_sim{n}]).([n_syn_g{p}]){k}{i,j},'S'); % remove field for data size
                             end
                             
                             %% EXTRACT AND CALCULATE DEGREES/RADIANS PHASE DATA FOR HIPPOCAMPAL NEURONS
                             hip_spikes(:,2) = hip_spikes(:,2) + par{k}{i,j}{t}.Hip_r_phase_n(n)*1000; % adjust for random Theta phase
                             phase_deg = [hip_spikes(:,1) (hip_spikes(:,2)/250 - floor(hip_spikes(:,2)/250))*360]; % degrees
                             phase_rad = [hip_spikes(:,1) phase_deg(:,2)*(pi/180)]; % radians
                             FR.PHASE.([n_sim{n}]).([n_syn_g{p}]){k}{i,j}{t} = [phase_deg phase_rad(:,2)]; % store data
                         end
                    end
                    %% EXTRACT AFTER-LEARNING (AL) WEIGHTS FOR HIPPOCAMPAL GROUPS
                    for p=1:4
                        WC.AL_W.([n_syn_g{p}]){k}{i,j}(t,1) = WC.DL.([n_syn_g{p}]){k}{i,j}(end);
                        if(t == vars.trials)
                            WC.AL_W.([n_syn_g{p}]){k}{i,j}(t,1) = WC.DL.([n_syn_g{p},'_raw']){k}{i,j}(end);
                        end                                                             
                        WC.REC_W.([n_syn_g{p}]){k}{i,j}(t,1) = mean([WC.P_AL.([n_syn_g{p}]){k}{i,j}(end) WC.NP_AL.([n_syn_g{p}]){k}{i,j}(end)]);
                    end
                end
            end
       end
    end
    close(h1)
    % SAVE ANALYSIS DATA
    mkdir([vars.test_f '/Analysis'])
    save([vars.test_f '/Analysis/firing-rate-data.mat'],'FR','-v7.3') 
    save([vars.test_f '/Analysis/weight-change-data.mat'],'WC','-v7.3')
% end
