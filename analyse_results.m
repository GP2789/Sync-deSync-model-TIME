function [ ] = analyse_results( filename, savefilename, type, trials, phase_con, freqs, no_flicker )
% Analyse simulation data in the directory specified by filename. Analyse
% by type ('FR', 'WC' for Activity, Weight Change, respectively).
line_color = {[0.8 0.6 0.7; 0.8 0.6 0.7] [0 0 0; 0 0 0] [0.8 0.4 0; 0.8 0.4 0] [0.95 0.9 0.25; 0.95 0.9 0.25] [0.2 1 1; 0 1 1] [1 0.2 1; 1 0 1]...
    [1 1 0.2; 1 1 0] [0.85 0.2 0; 0.8500 0.3250 0.0980]};
con_label = {['0',char(176), 'phase offset'], ['90',char(176), ' phase offset'], ['180',char(176), ' phase offset'], ['270',char(176), ' phase offset']};
        
mkdir(['analyse_figures/' savefilename]);
for cond = 1:length(phase_con)
    %% INITIALISE DATA STRUCTURES
    load([filename{freqs == 4,cond} '/variables.mat']); par = vars.par;
    v1_L = max(1,length(vars.v1_vars)); 
    v2_L = max(1,length(vars.v2_vars));
    v3_L = max(1,length(vars.v3_vars));
    N = length(par{1}{1,1}{1}.sim_order_n);
    n_sim = cell(1,N); 
    n_syn_g = {'P','NP','NP_P','P_NP'};
    n_syn_bg = {'A --> V','V --> A'};
    n_g = {'visual','auditory' };
    
    for n=1:N
        n_sim{n} = strrep(par{1}{1,1}{1}.sim_order_n{n},'-','_');
    end
    n_sim{N+1} = 'All';
    %% AVERAGE FIRING RATE PLOTS OVER TRIALS
    if(sum(cellfun(@sum,strfind(type,'FR')))>0)
        load([vars.test_f '/Analysis/firing-rate-data.mat']); % load activity data
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
                    
                    %% EXTRACT & SMOOTH DATA
                    t=vars.trials;
                    x1 = -par{k}{i,j}{t}.pre_stim_length + 1: 1 : par{k}{i,j}{t}.stim_length;
                    sim_length = par{k}{i,j}{t}.pre_stim_length + par{k}{i,j}{t}.stim_length; bx = 10;
                    I = FR.I.DL.BG_HIP{k}{i,j} + FR.I.DL.NC_HIP{k}{i,j} ...
                        + FR.I.DL.HIP_HIP_WIT{k}{i,j} + FR.I.DL.HIP_HIP_BET{k}{i,j} ...
                        + FR.I.DL.ADP_HIP{k}{i,j};
                    I_F = zeros(1, sim_length);
                    I_NC_F = zeros(1, sim_length); I_BG_F = zeros(1, sim_length); I_ADP_F = zeros(1, sim_length);
                    I_Hip_W_F = zeros(1, sim_length); I_Hip_B_F = zeros(1, sim_length);
                    for l = 1:sim_length % smooth data with box 
                        I_F(l) = mean(I(max(1,l-bx):min(sim_length,l+bx)));
                        I_BG_F(l) = mean(FR.I.DL.BG_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                        I_NC_F(l) = mean(FR.I.DL.NC_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                        I_ADP_F(l) = mean(FR.I.DL.ADP_HIP{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                        I_Hip_W_F(l) = mean(FR.I.DL.HIP_HIP_WIT{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                        I_Hip_B_F(l) = mean(FR.I.DL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
    %                     I_Hip_BL(l) = mean(FR.I.NP_BL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx))) + ...
    %                         mean(FR.I.P_BL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                        I_Hip_AL(l) = mean(FR.I.NP_AL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx))) + ...
                            mean(FR.I.P_AL.HIP_HIP_BET{k}{i,j}(max(1,l-bx):min(sim_length,l+bx)));
                    end

                    x = -par{k}{i,j}{t}.pre_stim_length : FR.ACT.bin_width : par{k}{i,j}{t}.stim_length-FR.ACT.bin_width;
                    x = x / 1000; x1 = x1/1000;
                    FR.ACT.DL.P{k}{i,j}.M = mean(FR.ACT.DL.P{k}{i,j}.S);
                    FR.ACT.DL.P{k}{i,j}.SE = std(FR.ACT.DL.P{k}{i,j}.S)/sqrt(trials);
                    FR.ACT.DL.NP{k}{i,j}.M = mean(FR.ACT.DL.NP{k}{i,j}.S);
                    FR.ACT.DL.NP{k}{i,j}.SE = std(FR.ACT.DL.NP{k}{i,j}.S)/sqrt(trials);
                    DL_P = [FR.ACT.DL.P{k}{i,j}.M - FR.ACT.DL.P{k}{i,j}.SE; ...
                        FR.ACT.DL.P{k}{i,j}.M; ...
                        FR.ACT.DL.P{k}{i,j}.M + FR.ACT.DL.P{k}{i,j}.SE];
                    DL_NP = [FR.ACT.DL.NP{k}{i,j}.M - FR.ACT.DL.NP{k}{i,j}.SE; ...
                        FR.ACT.DL.NP{k}{i,j}.M; ...
                        FR.ACT.DL.NP{k}{i,j}.M + FR.ACT.DL.NP{k}{i,j}.SE];

                    %% PLOT DL ACTIVITY
                    for p = 1:2
                        if cond == 1
                           f_fr(p) = figure;
                        end
                        figure(f_fr(p));
                        hold on
                        temp_plot = eval(['DL_',n_syn_g{p}]);
                        r{p}(cond) = fill([x fliplr(x)],[temp_plot(1,:) fliplr(temp_plot(3,:))], line_color{cond}(1,:),'edgecolor',line_color{cond}(1,:)); alpha(0.2); %alpha(a); % plot mean distribution
                        plot(x,temp_plot(2,:),'color', line_color{cond}(2,:), 'linestyle', '-','linewidth',2); % plot average
                        if cond == 1
                            xlim([-2 3]);  title(['Hippocampal Activity During Learning: ',{[n_g{p} ' neurons']}])
                            ylabel('Activity (Hz)'); 
                            xlabel('Time (s)');     
                        end
                        ax = gca; ax.FontSize = 16;
                        ax.FontName = 'Arial';

                        if cond == length(phase_con)
                           legend(r{p},con_label,'EdgeColor','none','Location','NorthWest');
                           saveas(f_fr(p),['analyse_figures/' savefilename '/hipp_activity_DL_', n_g{p},'.tiff']); close(f_fr(p));
                        end
                        clear temp_plot
                    end
                end
            end
        end
    end 
end

x_freqspec1 = [1 1.04]; 
x_freqspec2 = [1.2 1.24]; 
x_freqspec3 = [1.4 1.44];
for num_freqs = 1:length(freqs)
    data_post = cell(1,2); %zeros(trials, length(phase_con));
    data_pre = data_post;

    for cond = 1:length(phase_con)
        %% AVERAGE WEIGHT CHANGE PLOTS
        if(sum(cellfun(@sum,strfind(type,'WC')))>0)
            load([filename{num_freqs,cond} '/variables.mat']); par = vars.par;
            load([vars.test_f '/Analysis/weight-change-data.mat']);
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

                        %% EXTRACT AND PLOT HIPPOCAMPAL WEIGHT CHANGE
                        % plot weight change between groups only
                        p = 2;
                        for w=1:2 % aud -> visual & vis -> aud
                            n=2; % DL
                            % weight change data over simulations
                            average_y = transpose(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2}]){k}{i,j});
                            ste_y = transpose(std(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j},0,2)/sqrt(trials));
                            recall_AL_cue.([n_syn_g{w}])(1,:) = average_y - ste_y;
                            recall_AL_cue.([n_syn_g{w}])(2,:) = average_y;
                            recall_AL_cue.([n_syn_g{w}])(3,:) = average_y + ste_y;
                            x = -par{k}{i,j}{1}.pre_stim_length+1 : par{k}{i,j}{1}.stim_length;
                            x = x/1000;
                            if length(phase_con) == 4
                                if cond == 1
                                    f_wc(w) = figure;
                                end
                                figure(f_wc(w));
                                hold on;
                                rwc{w}(cond) = fill([x fliplr(x)],[recall_AL_cue.([n_syn_g{w}])(1,:) fliplr(recall_AL_cue.([n_syn_g{w}])(3,:))], line_color{cond}(1,:),'edgecolor',line_color{cond}(1,:)); alpha(0.2); %alpha(a); % plot mean distribution
                                plot(x,recall_AL_cue.([n_syn_g{w}])(2,:),'color', line_color{cond}(2,:), 'linestyle', '-','linewidth',2); % plot average
                                if cond == 1
                                    xlim([-2 3]); ylim([0 0.7]); title(['Weight Change Between Groups: ', n_syn_bg{w}])
                                    ylabel('Synaptic Efficacy'); 
                                    xlabel('Time (s)');     
                                end
                                ax = gca; ax.FontSize = 16;
                                ax.FontName = 'Arial';

                                if cond == length(phase_con)
                                   if w == 1
                                      fill([x(stim_end - 250+1) x(stim_end - 250+1) x(stim_end) x(stim_end)],[0 0.7 0.7 0],[0.8 0.8 0.8],'FaceAlpha',.4,'EdgeAlpha',0)
                                   end
                                   legend(rwc{w},con_label,'EdgeColor','none','Location','NorthWest');
                                   saveas(f_wc(w),['analyse_figures/' savefilename '/weight_change_' n_g{w} '_freq' num2str(freqs(num_freqs)) '.tiff']); close(f_wc(w));
                                end
                            end
                        end
                        stim_onset = find(x==0);
                        stim_end = find(x==x(end));
                        for w = 1:2 %A--> V and V-->A
                            % post learning WC: average across post stimulus end - 250 ms
                            % baseline WC: average across prestimulus 250 ms to onset 
                            data_post{w}(:,cond) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_end - 250+1 : stim_end,:),1));       
                            data_pre{w}(:,cond) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_onset-250 : stim_onset -1,:),1));  
                        end

                    end
                end
            end
        end
        clear recall_AL_cueNP x
    end

    %% average across post stimulus 3 s
    for w = 1:2
        ste_pre{w} = std(data_pre{w})./sqrt(length(data_pre{w})); 
        mean_pre{w} = mean(data_pre{w});
        ste_post{w} = std(data_post{w})./sqrt(length(data_post{w})); 
        mean_post{w} = mean(data_post{w});

        %subtracting mean for comparing with Clouter et al., 2017 and Wang et al.,
        %2018
        data_pren{w} = data_pre{w} - mean(data_pre{w},2);
        data_postn{w} = data_post{w} - mean(data_post{w},2);
        ste_pren{w} = std(data_pren{w})./sqrt(length(data_pren{w})); 
        mean_pren{w} = mean(data_pren{w});
        ste_postn{w} = std(data_postn{w})./sqrt(length(data_postn{w})); 
        mean_postn{w} = mean(data_postn{w});
    end
    
    recall_cueNP_plt = [mean_post{1} - ste_post{1}; mean_post{1}; mean_post{1} + ste_post{1}];
    recall_cueNP_blplt = [mean_pre{1} - ste_pre{1}; mean_pre{1}; mean_pre{1} + ste_pre{1}];

    recall_cueNP_pltn = [mean_postn{1} - ste_postn{1}; mean_postn{1}; mean_postn{1} + ste_postn{1}];
    recall_cueNP_blpltn = [mean_pren{1} - ste_pren{1}; mean_pren{1}; mean_pren{1} + ste_pren{1}];

    e_NP = (recall_cueNP_pltn(3,:)-recall_cueNP_pltn(1,:))/2; %error bars

    recall_cueP_plt = [mean_post{2} - ste_post{2}; mean_post{2}; mean_post{2} + ste_post{2}];
    e_P = (recall_cueP_plt(3,:)-recall_cueP_plt(1,:))/2; %error bars

    if length(phase_con) == 4
        y_postNP(num_freqs,:) = [recall_cueNP_plt(2,1) mean(recall_cueNP_plt(2,2:4))]; %average across aync condition
        y_bl(num_freqs,:) = [recall_cueNP_blplt(2,1) mean(recall_cueNP_blplt(2,2:4))]; %average across aync condition
        ebNP(num_freqs,:) = [(recall_cueNP_plt(3,1)-recall_cueNP_plt(1,1))/2 (mean(recall_cueNP_plt(3,2:4))-mean(recall_cueNP_plt(1,2:4)))/2]; %error bars also average across asyn condition
        eb_bl(num_freqs,:) = [(recall_cueNP_blplt(3,1)-recall_cueNP_blplt(1,1))/2 (mean(recall_cueNP_blplt(3,2:4))-mean(recall_cueNP_blplt(1,2:4)))/2]; %error bars also average across asyn condition
    else
        y_postNP(num_freqs,:) = recall_cueNP_plt(2,:);
        y_bl(num_freqs,:) = recall_cueNP_blplt(2,:);
        ebNP(num_freqs,:) = (recall_cueNP_plt(3,:)-recall_cueNP_plt(1,:))/2;
        eb_bl(num_freqs,:) = (recall_cueNP_blplt(3,:)-recall_cueNP_blplt(1,:))/2;
        y_postP(num_freqs,:) = recall_cueP_plt(2,:);
        ebP(num_freqs,:) = (recall_cueP_plt(3,:)-recall_cueP_plt(1,:))/2;
    end
    
    %% only plot this for figure 3C
    if freqs(num_freqs) == 4 && length(phase_con)==4
        left_color = [0 0 0]; right_color = [0 0.45 0.7];
        fig = figure;
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        box off; hold on; set(gca,'fontname','arial','fontsize',16)

        x_con = [1 2.5 4 5.5];

        % load Clouter et al data
        data = readtable('experiment1_clouteretal.csv');

        phase = data.PhaseOffset;
        avgpr = data.AvgCorrect;

        idx0 = phase==0;
        idx90 = phase==90;
        idx180 = phase==180;
        idx270 = phase==270;
        clouter(:,1) = 1:24;
        clouter(:,2:5) = [avgpr(idx0,1) avgpr(idx90,1) avgpr(idx180,1) avgpr(idx270,1)];

        clouter(:, 2:5) = clouter(:,2:5) - mean(clouter(:,2:5),2);

        n_clouter = length(clouter);
        semc =  std(clouter(:,2:5))/sqrt(n_clouter);    
        theta_clouter(1,:) = mean(clouter(:,2:5)) - semc; % STE
        theta_clouter(2,:) = mean(clouter(:,2:5));
        theta_clouter(3,:) = mean(clouter(:,2:5)) + semc; % STE
        theta_clouter = theta_clouter*100;
        
        e_clouter = (theta_clouter(3,:)-theta_clouter(1,:))/2; %error bars

        % load Wang et al data
        dataw = readtable('experiment_single_trial_wangetal.csv');
        idx0 = dataw.Phase_offset_condition==0;
        idx90 = dataw.Phase_offset_condition==90;
        idx180 = dataw.Phase_offset_condition==180;
        idx270 = dataw.Phase_offset_condition==270;

        wang = [dataw.Accuracy(idx0,1) dataw.Accuracy(idx90,1) dataw.Accuracy(idx180,1) dataw.Accuracy(idx270,1)];

        wang = wang - mean(wang,2);
        n_wang = length(wang);
        semcw =  std(wang(:,1:4))/sqrt(n_wang);
        theta_wang(1,:) = mean(wang(:,1:4)) - semcw; % STE
        theta_wang(2,:) = mean(wang(:,1:4));
        theta_wang(3,:) = mean(wang(:,1:4)) + semcw; % STE
        theta_wang = theta_wang*100;
        
        e_wang = (theta_wang(3,:)-theta_wang(1,:))/2; %error bars

        yyaxis left;
        % ylim([35 65]);
        errorbar(x_con,theta_clouter(2,:),e_clouter,'MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','o',...
            'Color',[0 0 0],'LineStyle','none',...
            'LineWidth',1,...
            'CapSize',12,...
            'Color',[0 0 0]);
        hold on
        errorbar(x_con + 0.3,theta_wang(2,:),e_wang,'MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','s',...
            'Color',[0 0 0],'LineStyle','none',...
            'LineWidth',1,...
            'CapSize',12,...
            'Color',[0 0 0]);
        ylabel('Accuracy (%)')
        % ylim([35 60]);
        xlim([0.2 7]);
        xticks([0 x_con+0.3 7]);
        xticklabels({'',['0',char(176)],['90',char(176)],['180',char(176)],['270',char(176)],''});
        xlabel('Phase offset condition')

        yyaxis right; 
        errorbar(x_con + 0.6,recall_cueNP_pltn(2,:),e_NP,'MarkerSize',8,'MarkerFaceColor',[0 0.45 0.7],'Marker','^',...
            'Color',[0 0.45 0.7],'LineStyle','none',...
            'LineWidth',1,...
            'CapSize',12,...
            'Color',[0 0.45 0.7]);
        ylabel('Synaptic Efficacy')
        % ylim([-0.2 0.4]);

        legend('Data from Clouter et al., 2017','Data from Wang et al., 2018','Simulated Hippocampal Weights',...
            'FontSize',12,'EdgeColor','No');
        saveas(fig,['analyse_figures/' savefilename '/mean_wc_clouter_wang.tiff']); close(fig);
    end
end

%% plot figure 5A
if length(freqs) == 3 && length(phase_con) == 4
    fig_freqspec = figure;
    box off; hold on; 
    set(gca,'fontname','arial','fontsize',16);
    errorbar(x_freqspec1, y_bl(1,:),eb_bl(1,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0.7 0.7 0.7],'Marker','o',...
        'Color',[0.7 0.7 0.7],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0.7 0.7 0.7]);
    errorbar(x_freqspec1, y_postNP(1,:),ebNP(1,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','o',...
        'Color',[0 0 0],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0 0 0]);


    errorbar(x_freqspec2, y_bl(2,:),eb_bl(2,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0.7 0.7 0.7],'Marker','^',...
        'Color',[0.7 0.7 0.7],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0.7 0.7 0.7]);
    errorbar(x_freqspec2, y_postNP(2,:),ebNP(2,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','^',...
        'Color',[0 0 0],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0 0 0]);

    errorbar(x_freqspec3, y_bl(3,:),eb_bl(3,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0.7 0.7 0.7],'Marker','d',...
        'Color',[0.7 0.7 0.7],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0.7 0.7 0.7]);
    errorbar(x_freqspec3, y_postNP(3,:),ebNP(3,:),'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','d',...
        'Color',[0 0 0],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0 0 0]);

    label = [{'','S', 'A','S', 'A','S', 'A'}; {'','  \delta','', '  \theta','', '  \alpha',''}];
    ticklabel = strtrim(sprintf('%s\\newline%s\n', label{:}));

    xlim([0.9 1.54]);
    ylim([0.1 0.7]);
    xlabel('stimulus frequency')
    ylabel('Synaptic Efficacy')
    xticks([0.9 1 1.04 1.2 1.24 1.4 1.44]); %xticklabels({'','\deltaS',' A','\thetaS',' A','\alphaS',' A'});
    xticklabels(ticklabel);
    saveas(fig_freqspec,['analyse_figures/' savefilename '/mean_wc_freq_spec.tiff']); close(fig_freqspec);
    
% plot figure 6
elseif length(phase_con) > 4
    freq.label = {'delta','theta','alpha','beta', 'low gamma', 'high gamma'};
    freq.value = [1.652 4 10.472 18.335 41.236 71.771];

    for f = 1:length(freqs)
        fpc(f) = figure;
        set(fpc(f),'position',[0        0        2000       1000]);
        hold on
        x_con = [-5 -4 -3 -2 -1 0 1 2];
        errorbar(x_con,y_postP(f,:),ebP(f,:),'MarkerSize',8,'MarkerFaceColor',[0 0.45 0.7],'Marker','o',...
            'Color',[0 0.45 0.7],'LineWidth',1,'CapSize',12,...
            'Color',[0 0.45 0.7]);
        errorbar(x_con,y_postNP(f,:),ebNP(f,:),'MarkerSize',8,'MarkerFaceColor',[0.8 0.4 0],'Marker','o',...
            'Color',[0.8 0.4 0],'LineWidth',1,'CapSize',12,...
            'Color',[0.8 0.4 0]);
        xlabel({'                      V leads A (A leads V)                                   A leads V (V leads A)'; 'Phase offset condition'})
        ylabel('Synaptic Efficacy')
        xticks(x_con); %]);
        xlim([-5.5 2.5]);
        ylim([0.1 0.65]);
        xticklabels({['225',char(176),'(135',char(176),')'],['180',char(176)], ['135',char(176),'(225',char(176),')'],['90',char(176),'(270',char(176),')'],...
            ['45',char(176),'(315',char(176),')'],['0',char(176)],['45',char(176),'(315',char(176),')'],['90',char(176),'(270',char(176),')']});
        legend('Visual -> Auditory', 'Auditory -> Visual','EdgeColor','none');
        title({'Mean hippocampal weight change between groups',['Stimulus frequency: ' freq.label{freq.value == freqs(f)}]})
        ax = gca; ax.FontSize = 16;
        ax.FontName = 'Arial';
        saveas(fpc(f),['analyse_figures/' savefilename '/meanWC_wider_phasecon_', freq.label{freq.value == freqs(f)},'.tiff']); close(fpc(f));
    end
    
end

%% comparison between theta 0 vs. 180 vs. no flicker
if no_flicker
    x_fnf = [1 1.5 2];
    data_post = zeros(trials, 1);
    data_pre = data_post;

    %% AVERAGE WEIGHT CHANGE PLOTS
    if(sum(cellfun(@sum,strfind(type,'WC')))>0)
        load([filename{3} '/variables.mat']); par = vars.par;
        load([vars.test_f '/Analysis/weight-change-data.mat']);
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

                    %% EXTRACT AND PLOT HIPPOCAMPAL WEIGHT CHANGE
                    % plot weight change between groups only
                    p = 2;
                    w=1; % aud -> visual
                    n=2; % DL
                    % weight change data over simulations
                    average_y = transpose(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2}]){k}{i,j});
                    ste_y = transpose(std(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j},0,2)/sqrt(trials));
                    recall_AL_cueNP(1,:) = average_y - ste_y;
                    recall_AL_cueNP(2,:) = average_y;
                    recall_AL_cueNP(3,:) = average_y + ste_y;
                    x = -par{k}{i,j}{1}.pre_stim_length+1 : par{k}{i,j}{1}.stim_length;
                    x = x/1000;

                    % post learning WC: average across post stimulus end - 250 ms
                    % baseline WC: average across prestimulus 250 ms to onset 
                    stim_onset = find(x==0);
                    stim_end = find(x==x(end));
                    data_post(:,1) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_end - 250+1 : stim_end,:),1));       
                    data_pre(:,1) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_onset-250 : stim_onset -1,:),1));  

                end
            end
        end
    end
    clear recall_AL_cueNP x
    

    %% average across post stimulus 1.5 s
    ste_pre = std(data_pre)./sqrt(length(data_pre)); 
    mean_pre = mean(data_pre);
    ste_post = std(data_post)./sqrt(length(data_post)); 
    mean_post = mean(data_post);

    recall_cueNP_plt = [mean_post - ste_post; mean_post; mean_post + ste_post];
    recall_cueNP_blplt = [mean_pre - ste_pre; mean_pre; mean_pre + ste_pre];

    y_postNP(1,3) = recall_cueNP_plt(2,1); 
    y_bl(1,3) = recall_cueNP_blplt(2,1); 
    ebNP(1,3) = (recall_cueNP_plt(3,1)-recall_cueNP_plt(1,1))/2; 
    eb_bl(1,3) = (recall_cueNP_blplt(3,1)-recall_cueNP_blplt(1,1))/2;

    
    %% plot figure 5B
    fig_fnf = figure;
    box off; hold on; 
    set(gca,'fontname','arial','fontsize',16);
    errorbar(x_fnf, y_bl,eb_bl,'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0.7 0.7 0.7],'Marker','o',...
    'Color',[0.7 0.7 0.7],'LineWidth',1,...
    'CapSize',12,...
    'Color',[0.7 0.7 0.7]);
    errorbar(x_fnf, y_postNP,ebNP,'LineStyle','none','MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','o',...
        'Color',[0 0 0],'LineWidth',1,...
        'CapSize',12,...
        'Color',[0 0 0]);
    ylabel('Synaptic Efficacy')
    xlim([0.5 2.5]);
    ylim([0.1 0.7]);
    xticks([1 1.5 2]); xticklabels({['0' char(176)],['180' char(176)],'No Flicker'});
    saveas(fig_fnf,['analyse_figures/' savefilename '/mean_wc_theta_vs_noflicker.tiff']); close(fig_fnf);
end

