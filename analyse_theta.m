function [ filename ] = analyse_theta( filename, n_trials, simulate, condition, het_STDP)

%% initialise
clearvars global par; global par;
set_parameters;
par.n_Items = 2;
par.Hip_per_item = [5 5];
par.n_Hip = sum(par.Hip_per_item);
par.n_NC = 0;
par.network_size = par.n_NC + par.n_Hip;
par.nG = {'Hip'};
par.STDP_hetero = het_STDP;

par.Hip_r_phase = 0; par.NC_r_phase = 0;
par.stim_type = 'DC';         % set to DC input
par.stim_TS = 100;            % frequency of burst
if(contains(condition,'single'))
    burst_n = [0 1 2 3 4];      % number of shocks in burst
    stim_phase = {0, 1};      % theta phase of burst
    filename = [filename 'single-burst'];
    stim_n = {'trough','peak'};
elseif(contains(condition,'both'))
    burst_n = 4;              % number of shocks in burst
    stim_phase = {[0 1]};     % theta phase of burst
    filename = [filename 'multi-burst'];
    stim_n = {'both'};
end
par.stimulus_strength = 3; % strength of burst
par.pre_stim_length = 250;
par.stim_length = 500;
par.post_stim_length = 250;
par.sim_length = par.pre_stim_length + par.stim_length + par.post_stim_length;
par.Hip_amp = 0.01;
par.Hip_ADP = false;
par.rand_phase = false;
par.Hip_SW = 0;
par.Hip_SR = 0;
par.Hip_inter_syn_str = 0.01;
par.Hip_intra_syn_str = 0.01;
par.weight_max = 0.01;
par.stim_mod = false;
par.Hip_reset = false;

% set sim_length and sim_type

par.sim_type = 'HL';

stim_dW = nan(length(burst_n), length(stim_phase), n_trials);
unstim_dW = nan(length(burst_n), length(stim_phase), n_trials);

HET_stim       = cell(length(stim_phase), length(burst_n));
HET_unstim     = cell(length(stim_phase), length(burst_n));
LTP_stim       = cell(length(stim_phase), length(burst_n));
LTP_unstim     = cell(length(stim_phase), length(burst_n));
LTD_stim       = cell(length(stim_phase), length(burst_n));
LTD_unstim     = cell(length(stim_phase), length(burst_n));
wc_pos_stim   = cell(length(stim_phase), length(burst_n));
wc_pos_unstim = cell(length(stim_phase), length(burst_n));
wc_neg_stim   = cell(length(stim_phase), length(burst_n));
wc_neg_unstim = cell(length(stim_phase), length(burst_n));
stim_p_t       = cell(length(stim_phase), length(burst_n));
unstim_p_t     = cell(length(stim_phase), length(burst_n));

% make directories for simulation
filename = [filename '_' int2str(n_trials) 'T'];
if(exist(filename,'dir')~=7); mkdir(filename); end
if(exist([filename '/Data'],'dir')~=7); mkdir([filename '/Data']); end

%% simulate
%fprintf('simulating trials ... \n');
for t = 1:n_trials
    if(exist([filename '/T' int2str(t)],'dir')~=7); mkdir([filename '/T' int2str(t)]); end
    for i = 1:length(stim_phase)
        par.stim_phase = stim_phase{i};
        for j = 1:length(burst_n)
            par.burst_n = burst_n(j);
            T = 0:1/1000:par.sim_length/1000 - 1/1000;
            
            if(simulate)
                create_network(); % create network
                sim_stats = simulate_network(); % simulate network
            else
                %fprintf('\nloading trial %.0f/%.0f\n', t, n_trials)
                load([filename '/Data/T' int2str(t) '_' stim_n{i} '_' num2str(burst_n(j)) '-shocks.mat'])
            end

            stim_pathways = zeros(par.network_size, par.network_size);
            stim_pathways( 1:par.Hip_per_item(1), 1:par.Hip_per_item(1) ) = 1;
            unstim_pathways = zeros(par.network_size, par.network_size);
            unstim_pathways( par.Hip_per_item(1)+1:end, 1:par.Hip_per_item(1) ) = 1;
            %unstim_pathways( 1:par.Hip_per_item(1), par.Hip_per_item(1)+1:end ) = 1;
            %unstim_pathways( par.Hip_per_item(1)+1:end, par.Hip_per_item(1)+1:end ) = 1;
            stim_p = find(stim_pathways & sim_stats.weight_matrix_STDP);
            unstim_p = find(unstim_pathways & sim_stats.weight_matrix_STDP);

            for t2 = 1:length(T)
                wm = sim_stats.I_REC.HET(:,:,t2); 
                HET_stim{i,j}(t,t2) = mean(wm(stim_p))/(par.weight_max)*100; 
                HET_unstim{i,j}(t,t2) = mean(wm(unstim_p))/(par.weight_max)*100;
                wm = sim_stats.I_REC.LTP(:,:,t2); 
                LTP_stim{i,j}(t,t2) = mean(wm(stim_p))/(par.weight_max)*100; 
                LTP_unstim{i,j}(t,t2) = mean(wm(unstim_p))/(par.weight_max)*100;
                wm = sim_stats.I_REC.LTD(:,:,t2); 
                LTD_stim{i,j}(t,t2) = mean(wm(stim_p))/(par.weight_max)*100; 
                LTD_unstim{i,j}(t,t2) = mean(wm(unstim_p))/(par.weight_max)*100;
                
                wm = sim_stats.I_REC.wc_pos(:,:,t2); 
                wc_pos_stim{i,j}(t,t2) = mean(wm(stim_p)); 
                wc_pos_unstim{i,j}(t,t2) = mean(wm(unstim_p));
                wm = sim_stats.I_REC.wc_neg(:,:,t2); 
                wc_neg_stim{i,j}(t,t2) = mean(wm(stim_p)); 
                wc_neg_unstim{i,j}(t,t2) = mean(wm(unstim_p));
                
                wm = sim_stats.weight_matrix(:,:,t2); 
                stim_p_t{i,j}(t,t2) = mean(wm(stim_p)); 
                unstim_p_t{i,j}(t,t2) = mean(wm(unstim_p)); clear wm;
            end
            
            stim_p_t{i,j}(t,:) = (stim_p_t{i,j}(t,:) / mean(stim_p_t{i,j}(t,1:par.pre_stim_length))) * 100 - 100;
            unstim_p_t{i,j}(t,:) = (unstim_p_t{i,j}(t,:) / mean(unstim_p_t{i,j}(t,1:par.pre_stim_length))) * 100 - 100;
            
            save([filename '/Data/T' int2str(t) '_' stim_n{i} '_' num2str(burst_n(j)) '-shocks.mat'], 'sim_stats');
            
            stim_dW(j,i,t) = stim_p_t{i,j}(t,par.pre_stim_length + 400);
            unstim_dW(j,i,t) = unstim_p_t{i,j}(t,par.pre_stim_length + 400);
        end
    end
    fprintf('\b\b\b\b\b%3.0f%%\n', t/n_trials*100)
end

%% plot 1
left_color = [0 0 0]; right_color = [0 0 0];
if(contains(condition,'single')); x_lim = [0 0.5];
else; x_lim = [0 max(T)];
end

for i = 1:length(stim_phase)
    for j = 1:length(burst_n)
        load([filename '/Data/T' int2str(t) '_' stim_n{i} '_' num2str(burst_n(j)) '-shocks.mat'])
        
        fig = figure(1); 
        if(het_STDP); set(fig, 'position', [0 0 600 900],'defaultAxesColorOrder',[left_color; right_color]); 
        else; set(fig, 'position', [0 0 600 500],'defaultAxesColorOrder',[left_color; right_color]); 
        end
        if(het_STDP); subplot(4,1,1); else; subplot(3,1,1); end
        hold on; box on; set(gca,'fontsize',14,'fontname','Arial')
        title('Raster Plot & \theta-STDP'); 
        yyaxis right; clear l;
        %l(1) = plot(T, sim_stats.I_REC.HIP_THETA, 'b-', 'linewidth',1);
        l(1) = plot(T, 1-sim_stats.I_REC.HIP_THETA, 'color', right_color, 'linewidth',2);
        ylim([0 2]); yticks([0 1]); ylabel('\theta-amplitude')
        yyaxis left;
        l(2) = scatter(sim_stats.spike_detector(:,2)/1000,sim_stats.spike_detector(:,1),'k.');
        xticklabels({}); ylabel('neuron ID'); xlim(x_lim);
        ylim([-par.network_size, par.network_size]); yticks([0 par.network_size/2 par.network_size])
        legend(l, 'Hipp.-\theta', 'Spikes', ...
            'location','northwest','orientation','vertical');

        
        y_lim = ceil(max(abs([mean(HET_stim{i,j},1) mean(LTP_stim{i,j},1) mean(LTD_stim{i,j},1)...
            mean(HET_unstim{i,j},1) mean(LTP_unstim{i,j},1) mean(LTD_unstim{i,j},1)]))/10)*10;

        if(het_STDP); subplot(4,1,2); title('Stim. Plasticity'); 
        else; subplot(3,1,2); title('Plasticity');  end 
        hold on; box on; set(gca,'fontsize',14,'fontname','Arial')
         clear l;
        yyaxis left
        if(het_STDP); dpdt = mean(LTD_stim{i,j},1) + mean(LTP_stim{i,j},1) + mean(HET_stim{i,j},1); 
        else; dpdt = mean(LTD_stim{i,j},1) + mean(LTP_stim{i,j},1);
        end
        l(1) = plot(T,dpdt,'k-','linewidth',2);
        %plot(T,mean(LTP_stim{i,j},1),'k-','linewidth',2);
        %if(het_STDP); l(3) = plot(T,mean(HET_stim{i,j},1),'k-','linewidth',2); end
        plot([0 max(T)], [0 0],'k:','linewidth',1)
        xticklabels({}); xlim(x_lim);
        ylabel('d\rho/dt (%)'); ylim([-1 1] * ceil(max(abs(ylim))/10)*10);
        
        yyaxis right
        plot([0 max(T)], [1 1]*par.T_p, 'b--')
        plot([0 max(T)], [1 1]*-par.T_d, 'r--')
        l(2) = fill([T fliplr(T)], [mean(wc_pos_stim{i,j},1) zeros(1,length(T))],...
            [0 0 1],'edgecolor','none');
        l(3) = fill([T fliplr(T)], [-mean(wc_neg_stim{i,j},1) zeros(1,length(T))],...
            [1 0 0],'edgecolor','none');
        ylabel('potential plasticity'); alpha(0.35);ylim([-1 1] * ceil(max(abs(ylim))));
        yticks([-par.T_d 0 par.T_p]);
        legend(l, 'd\rho/dt', 'F_L_T_P','F_L_T_D', 'location','west','orientation','vertical');

        if(het_STDP)
            subplot(4,1,3); hold on; box on; set(gca,'fontsize',14,'fontname','Arial')
            title('No Stim. Plasticity');
            yyaxis left
            if(het_STDP); dpdt = mean(LTD_unstim{i,j},1) + mean(LTP_unstim{i,j},1) + mean(HET_unstim{i,j},1); 
            else; dpdt = mean(LTD_unstim{i,j},1) + mean(LTP_unstim{i,j},1);
            end
            l(1) = plot(T,dpdt,'k-','linewidth',2);
            if(het_STDP); l(3) = plot(T,mean(HET_unstim{i,j},1),'k-','linewidth',2); end
            plot([0 max(T)], [0 0],'k:','linewidth',1)
            xticklabels({}); xlim(x_lim);
            ylabel('d\rho/dt (%)'); ylim([-1 1] * ceil(max(abs(ylim))/10)*10);

            yyaxis right
            plot([0 max(T)], [1 1]*par.T_p, 'b--')
            plot([0 max(T)], [1 1]*-par.T_d, 'r--')
            l(2) = fill([T fliplr(T)], [mean(wc_pos_unstim{i,j},1) zeros(1,length(T))],...
                [0 0 1],'edgecolor','none');
            l(3) = fill([T fliplr(T)], [-mean(wc_neg_unstim{i,j},1) zeros(1,length(T))],...
                [1 0 0],'edgecolor','none');
            ylabel('potential plasticity'); alpha(0.35); ylim([-1 1] * max(abs(ylim)));
            legend(l, 'd\rho/dt', 'F_L_T_P','F_L_T_D', 'location','west','orientation','vertical');
        end

        if(het_STDP); subplot(4,1,4); else; subplot(3,1,3); end
        hold on; box on;set(gca,'fontsize',14,'fontname','Arial')
        title('Change in Synaptic Efficacy'); clear l; 
        l(1) = plot(T,mean(stim_p_t{i,j},1),'k-','linewidth',2);
        if(het_STDP); l(2) = plot(T,mean(unstim_p_t{i,j},1),'k--','linewidth',2); end
        plot([0 max(T)],[0 0],'k:')
        ylabel('\rho (% change)'); ylim([-50 100]);%ylim([-1 1] * ceil(max(abs(ylim))/10)*10);
        if(het_STDP); legend(l,'Stim.', 'No Stim.', 'location','west','orientation','vertical');
        end
        xlabel('time (s)'); xlim(x_lim);
        
        saveas(fig,[filename '/weight_change_' stim_n{i} '_' num2str(burst_n(j)) '-shocks.tiff']); close(fig);
        
        fig = figure(); set(fig,'position',[0 0 500 300])
        hold on; box on; set(gca,'fontsize',14,'fontname','Arial')
        xlabel('\rho (% change)'); %title('Summary of long-term synaptic plasticity during \theta')
        burst_at_trough_stim   = stim_p_t{i,j}(:,par.pre_stim_length+200);
        burst_at_trough_unstim = -unstim_p_t{i,j}(:,par.pre_stim_length+200);
        burst_at_peak_stim     = (burst_at_trough_stim - stim_p_t{i,j}(:,end))./burst_at_trough_stim*100;
        burst_at_peak_unstim   = -unstim_p_t{i,j}(:,end) - burst_at_trough_unstim;
        
        boxplot(fliplr([burst_at_trough_stim, burst_at_trough_unstim, ...
            burst_at_peak_stim, burst_at_peak_unstim]), 'labels', ...
            fliplr({'Stim. LTP', 'No Stim. Het. LTD', 'Stim. LTD', 'No Stim. Het. LTD'}),...
            'colors','k','orientation','horizontal','symbol','');
        ax = gca(); bph = ax.Children; bpchil = bph.Children;
        xlim([floor(min([bpchil(3*4:4*4).XData])/10)*10 ceil(max([bpchil(4*4:5*4).XData])/10)*10]);
        %plot([0 0], ylim, 'k:'); 
        plot(xlim,[1 1]*diff(ylim)/2+.5,'k:')
        try
            xticks(unique(sort(round([median(burst_at_trough_stim) median(burst_at_trough_unstim)...
                median(burst_at_peak_stim) median(burst_at_peak_unstim) 0]))))
        catch
        end
        %ytickangle(45);
        saveas(fig,[filename '/weight_change_theta_' stim_n{i} '_' num2str(burst_n(j)) '-shocks.tiff']); close(fig);
    end
end

%% plot 2
if(contains(condition,'single'))
    fig = figure(); set(fig,'position',[0 0 450 400])

    zero_shock = squeeze(stim_dW(burst_n==0,:,:)); zero_shock = zero_shock(:);
    peak_shocks = squeeze(stim_dW(burst_n~=0,2,:));
    trough_shocks = squeeze(stim_dW(burst_n~=0,1,:));

    peak_m = mean(peak_shocks,2);
    trough_m = mean(trough_shocks,2);
    zero_m = mean(zero_shock);
    peak_sde = std(peak_shocks,[],2) / sqrt(n_trials);
    trough_sde = std(trough_shocks,[],2) / sqrt(n_trials);
    zero_sde = std(zero_shock) / sqrt(n_trials);
    m = [flipud(peak_m); zero_m; trough_m];
    sde = [flipud(peak_sde); zero_sde; trough_sde];

    b1 = bar(m,'facecolor','w','edgecolor','k','linewidth',2);
    b1.FaceColor = 'flat';
    b1.CData(1:length(burst_n),:) = repmat([0 0 0],[length(burst_n),1]);
    set(gca,'XTickLabel',[cellfun(@num2str,num2cell(fliplr(burst_n)),'uniformoutput',false) ...
        cellfun(@num2str,num2cell(burst_n(2:end)),'uniformoutput',false)], 'FontSize',14,'FontName','Arial')
    hold on;
    errorbar(1:length(m), m, sde, sde,'color','k','linestyle','none','linewidth',2)
    ylim([-50 100]); yticks(-50:25:100); yticklabels({'-50','','0','','50','','100'});
    ylabel('% Weight Change'); xlabel('Number of Shocks in Burst');
   %title('Dependence of synaptic modifications on number of shocks in a single burst');
    saveas(fig,[filename '/weight_change.tiff']); close(fig);
end

end

