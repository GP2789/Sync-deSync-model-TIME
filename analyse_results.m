function [ fstat_clouter, fstat_wang ] = analyse_results( filename, savefilename, type, trials, phase_con)
% Analyse simulation data in the directory specified by filename. Analyse
% by type ('FR', 'WC' for Activity, Weight Change, respectively).
% generate model comparisons results for comparing STDP model with data
% from clouter et al. and wang et al.
line_color = {[0.8 0.6 0.7; 0.8 0.6 0.7] [0 0 0; 0 0 0] [0.8 0.4 0; 0.8 0.4 0] [0.95 0.9 0.25; 0.95 0.9 0.25] [0.2 1 1; 0 1 1] [1 0.2 1; 1 0 1]...
    [1 1 0.2; 1 1 0] [0.85 0.2 0; 0.8500 0.3250 0.0980]};
con_label = {['0',char(176), 'phase offset'], ['90',char(176), ' phase offset'], ['180',char(176), ' phase offset'], ['270',char(176), ' phase offset']};
        
mkdir(['analyse_figures/' savefilename]);
for cond = 1:length(phase_con)
    %% INITIALISE DATA STRUCTURES
    load([filename{2,cond} '/variables.mat']); par = vars.par;
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
                    I_Hip_BL = zeros(1, sim_length); I_Hip_AL = zeros(1, sim_length);
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

                    %% PLOT DL ACTIVITY
                    if cond == 1
                       f_fr = figure;
                    end
                    hold on
                    r(cond) = fill([x fliplr(x)],[DL_P(1,:) fliplr(DL_P(3,:))], line_color{cond}(1,:),'edgecolor',line_color{cond}(1,:)); alpha(0.2); %alpha(a); % plot mean distribution
                    plot(x,DL_P(2,:),'color', line_color{cond}(2,:), 'linestyle', '-','linewidth',2); % plot average
                    if cond == 1
                        xlim([-2 3]);  title('Hippocampal Activity During Learning: visual neurons')
                        ylabel('Activity (Hz)'); 
                        xlabel('Time (s)');     
                    end
                    ax = gca; ax.FontSize = 16;

                    if cond == length(phase_con)
                       legend(r,con_label,'EdgeColor','none','Location','NorthWest');
                       saveas(f_fr,['analyse_figures/' savefilename '/hipp_activity_DL.tiff']); close(f_fr);
                    end
                    
                end
            end
        end
    end 
end

y_post = zeros(size(filename,1),4);
y_bl = y_post;
eb = y_post;
eb_bl = eb;
for num_models = 1:size(filename,1)
    data_post = zeros(trials, length(phase_con));
    data_pre = data_post;

    for cond = 1:length(phase_con)
        %% AVERAGE WEIGHT CHANGE PLOTS
        if(sum(cellfun(@sum,strfind(type,'WC')))>0)
            load([filename{num_models,cond} '/variables.mat']); par = vars.par;
            load([filename{num_models,cond} '/Analysis/weight-change-data.mat']);
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
%                         if length(phase_con) == 4
%                             if cond == 1
%                                 f_wc = figure;
%                             end
%                             hold on;
%                             rwc(cond) = fill([x fliplr(x)],[recall_AL_cueNP(1,:) fliplr(recall_AL_cueNP(3,:))], line_color{cond}(1,:),'edgecolor',line_color{cond}(1,:)); alpha(0.2); %alpha(a); % plot mean distribution
%                             plot(x,recall_AL_cueNP(2,:),'color', line_color{cond}(2,:), 'linestyle', '-','linewidth',2); % plot average
%                             if cond == 1
%                                 xlim([-2 3]); ylim([0 0.7]); title('Weight Change Between Groups: A --> V Synapses')
%                                 ylabel('Synaptic Efficacy'); 
%                                 xlabel('Time (s)');     
%                             end
%                             ax = gca; ax.FontSize = 16;
% 
%                             if cond == length(phase_con)
%                                fill([x(stim_end - 250+1) x(stim_end - 250+1) x(stim_end) x(stim_end)],[0 0.7 0.7 0],[0.8 0.8 0.8],'FaceAlpha',.4,'EdgeAlpha',0)
%                                legend(rwc,con_label,'EdgeColor','none','Location','NorthWest');
%                                saveas(f_wc,['analyse_figures/' savefilename '/weight_change_freq' num2str(freqs(num_freqs)) '.tiff']); close(f_wc);
%                             end
%                         end
                        % post learning WC: average across post stimulus end - 250 ms
                        % baseline WC: average across prestimulus 250 ms to onset 
                        stim_onset = find(x==0);
                        stim_end = find(x==x(end));
                        data_post(:,cond) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_end - 250+1 : stim_end,:),1));       
                        data_pre(:,cond) = transpose(mean(WC.([n_sim{n}]).([n_syn_g{w+(p-1)*2},'_raw']){k}{i,j}(stim_onset-250 : stim_onset -1,:),1));  

                    end
                end
            end
        end
        clear recall_AL_cueNP x
    end

    %% average across post stimulus 3 s
    ste_pre = std(data_pre)./sqrt(length(data_pre)); 
    mean_pre = mean(data_pre);
    ste_post = std(data_post)./sqrt(length(data_post)); 
    mean_post = mean(data_post);

    %subtracting mean for comparing with Clouter et al., 2017 and Wang et al.,
    %2018
    data_pren = data_pre - mean(data_pre,2);
    data_postn = data_post - mean(data_post,2);
    ste_pren = std(data_pren)./sqrt(length(data_pren)); 
    mean_pren = mean(data_pren);
    ste_postn = std(data_postn)./sqrt(length(data_postn)); 
    mean_postn = mean(data_postn);

    recall_cueNP_plt = [mean_post - ste_post; mean_post; mean_post + ste_post];
    recall_cueNP_blplt = [mean_pre - ste_pre; mean_pre; mean_pre + ste_pre];

    recall_cueNP_pltn = [mean_postn - ste_postn; mean_postn; mean_postn + ste_postn];
    recall_cueNP_blpltn = [mean_pren - ste_pren; mean_pren; mean_pren + ste_pren];

    e_NP = (recall_cueNP_pltn(3,:)-recall_cueNP_pltn(1,:))/2; %error bars
   
    y_post(num_models,:) = recall_cueNP_plt(2,:);
    y_bl(num_models,:) = recall_cueNP_blplt(2,:);
    eb(num_models,:) = (recall_cueNP_plt(3,:)-recall_cueNP_plt(1,:))/2;
    eb_bl(num_models,:) = (recall_cueNP_blplt(3,:)-recall_cueNP_blplt(1,:))/2;
end

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

n_clouter = length(clouter);
semc =  std(clouter(:,2:5))/sqrt(n_clouter);    
theta_clouter(1,:) = mean(clouter(:,2:5)) - semc; % STE
theta_clouter(2,:) = mean(clouter(:,2:5));
theta_clouter(3,:) = mean(clouter(:,2:5)) + semc; % STE
theta_clouter = theta_clouter * 100;

e_clouter = (theta_clouter(3,:)-theta_clouter(1,:))/2; %error bars

% load Wang et al data
dataw = readtable('experiment_single_trial_wangetal.csv');
idx0 = dataw.Phase_offset_condition==0;
idx90 = dataw.Phase_offset_condition==90;
idx180 = dataw.Phase_offset_condition==180;
idx270 = dataw.Phase_offset_condition==270;

wang = [dataw.Accuracy(idx0,1) dataw.Accuracy(idx90,1) dataw.Accuracy(idx180,1) dataw.Accuracy(idx270,1)];

n_wang = length(wang);
semcw =  std(wang(:,1:4))/sqrt(n_wang);
theta_wang(1,:) = mean(wang(:,1:4)) - semcw; % STE
theta_wang(2,:) = mean(wang(:,1:4));
theta_wang(3,:) = mean(wang(:,1:4)) + semcw; % STE
theta_wang = theta_wang * 100;

e_wang = (theta_wang(3,:)-theta_wang(1,:))/2; %error bars


%% model fit to Clouter et al.
x = 1:4;
% design matrix for full model
X1 = [ ones(4,1) y_post(1,:)' ];

% fit to Clouter et al. exp 1
b1 = (X1'*X1)\(X1'*theta_clouter(2,:)');

% design matrix for STDP only model
X2 = [ ones(4,1) y_post(2,:)' ];

% fit to Clouter et al. exp 1
b12 = (X2'*X2)\(X2'*theta_clouter(2,:)');

fig_clouter = figure;
box off; hold on; 
set(gca,'fontname','arial','fontsize',16);
errorbar(x,theta_clouter(2,:),e_clouter,'MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','o',...
    'Color',[0 0 0],'LineStyle','none',...
    'LineWidth',1,...
    'CapSize',12,...
    'Color',[0 0 0]);
ylabel('Accuracy (%)')
xlim([0.8 4.2]);
xticks([0 1 2 3 4 5]);
xticklabels({'',['0',char(176)],['90',char(176)],['180',char(176)],['270',char(176)],''});
xlabel('Phase offset condition')
plot(x, X1*b1,'Color',[0.35 0.7 0.9],'LineWidth',2);
plot(x, X2*b12,'Color',[0.93 0.69 0.13],'LineWidth',2);
legend('Data from Clouter et al., 2017',...
    ['Mean hippocampal weights predicted', newline, 'by the full model'],...
    ['Mean hippocampal weights predicted', newline, 'by the model with Theta-phase-learning only']...
    ,'FontSize',12,'EdgeColor','No','Color','none');
saveas(fig_clouter,['analyse_figures/' savefilename '/FullvsTheta_fit_clouter.tiff']); close(fig_clouter);

% model comparison
% empirical data mean
mc = theta_clouter(2,:);
% predicted by full model and STDP only model
m1 = X1*b1;
m2 = X2*b12;
% rsq_m1 = 1 - sum((mc - m1').^2)/sum((mc - mean(mc)).^2);
% rsq_m2 = 1 - sum((mc - m2').^2)/sum((mc - mean(mc)).^2);

% sum of squared errors for full model and STDP only
sse1 = sum((mc' - m1).^2);
sse2 = sum((mc' - m2).^2);

% statistic F value for comparison between the two models
fstat_clouter = (sse2 - sse1)/(sse1/3);


%% fit to Wang et al.
b2 = (X1'*X1)\(X1'*theta_wang(2,:)');

b22 = (X2'*X2)\(X2'*theta_wang(2,:)');

fig_wang = figure;
box off; hold on; 
set(gca,'fontname','arial','fontsize',16);
errorbar(x,theta_wang(2,:),e_wang,'MarkerSize',8,'MarkerFaceColor',[0 0 0],'Marker','o',...
    'Color',[0 0 0],'LineStyle','none',...
    'LineWidth',1,...
    'CapSize',12,...
    'Color',[0 0 0]);
ylabel('Accuracy (%)')
xlim([0.8 4.2]);
% ylim([35 60]);
xticks([0 1 2 3 4 5]);
xticklabels({'',['0',char(176)],['90',char(176)],['180',char(176)],['270',char(176)],''});
xlabel('Phase offset condition')
hold on
plot(x, X1*b2,'Color',[0.35 0.7 0.9],'LineWidth',2);
plot(x, X2*b22,'Color',[0.93 0.69 0.13],'LineWidth',2);
legend('Data from Wang et al., 2018',...
    ['Mean hippocampal weights predicted', newline, 'by the full model'],...
    ['Mean hippocampal weights predicted', newline, 'by the model with Theta-phase-learning only']...
    ,'FontSize',12,'EdgeColor','No','Color','none');
saveas(fig_wang,['analyse_figures/' savefilename '/FullvsTheta_fit_wang.tiff']); close(fig_wang);

% model comparison
% empirical data mean
mw = theta_wang(2,:);
% predicted by full model and STDP only model
m1 = X1*b2;
m2 = X2*b22;
% rsq_mw1 = 1 - sum((mw - m1').^2)/sum((mw - mean(mw)).^2);
% rsq_mw2 = 1 - sum((mw - m2').^2)/sum((mw - mean(mw)).^2);

% sum of squared errors for full model and STDP only
ssew1 = sum((mw' - m1).^2);
ssew2 = sum((mw' - m2).^2);

% statistic F value for comparison between the two models
fstat_wang = (ssew2 - ssew1)/(ssew1/3);