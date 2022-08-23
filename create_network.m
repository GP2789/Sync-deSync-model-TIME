function [  ] = create_network( )
% Create and connect neural network based on parameters specified in set_parameters.

%% DECLARE GLOBAL VARIABLES
global par V_m I ADP osc MUA delay ref_period weight_matrix weight_matrix_STDP wc_pos wc_neg;

%% INITIALISE DATA STRUCTURES
% neuron data storage
V_m = zeros(par.network_size, par.sim_length); % keeps track of voltages through time
par.V_ref = zeros(par.network_size, 1); % stores refractory period data for neurons
par.tau_syn = ones(par.network_size, par.network_size); % stores synaptic time constants
ref_period = ones(par.network_size, par.sim_length); % tracks current refractory period status through time
I.SYN = zeros(par.network_size, par.sim_length); % input for spikes
I.OSC = zeros(par.network_size, par.sim_length); % input for static oscillations
I.ADP = zeros(par.network_size, par.sim_length); % input for ADP functions
I.DC = zeros(par.network_size, par.sim_length); % input for DC functions
I.L = zeros(par.network_size, par.sim_length); % input for leak functions
MUA = zeros(par.network_size, par.sim_length); % tracks spike events through time

% synapse data storage
weight_matrix = zeros(par.network_size, par.network_size, par.sim_length); % tracks weights through time
wc_pos = zeros(par.network_size, par.network_size, par.sim_length); % tracks decaying + weight change through time
wc_neg = zeros(par.network_size, par.network_size, par.sim_length); % tracks decaying - weight change through time
weight_matrix_STDP = zeros(par.network_size, par.network_size); % informs on whether synapses have learning enabled (binary)
delay = zeros(par.network_size, par.network_size); % stores delay information on synapses

% set synaptic time constants
par.tau_syn(1:par.n_NC,1:par.n_NC) = par.NC_intra_tau_syn;
par.tau_syn(1:par.n_NC,1+par.n_NC:par.network_size) = par.NC_Hip_intra_tau_syn;
par.tau_syn(1+par.n_NC:par.network_size,1:par.n_NC) = par.Hip_NC_intra_tau_syn;
par.tau_syn(1+par.n_NC:par.network_size,1+par.n_NC:par.network_size) = par.Hip_intra_tau_syn;

V_m(:,1) = par.V_m_initial; % initial voltage starting value

t = 0:1/1000:par.sim_length/1000-0.001;
c = 0; nG = par.nG;

%% CREATE AND CONNECT NEURONS
% loops through neuron groups to extract and set intrinsic and network properties defined in set_parameters.m
for g=1:length(nG)
    %% create & connect oscillation sine wave
    if( par.rand_phase ) % random phase option
        par.([nG{g} '_r_phase']) =  round(rand()*(1000/par.([nG{g} '_freq'])))/1000;
    else
        par.([nG{g} '_r_phase']) = 0;
    end
    % create static oscillations
    osc.(nG{g})  = cos(2*pi*par.([nG{g} '_freq'])*(t+par.([nG{g} '_r_phase'])))*par.([nG{g} '_amp']);
    if(strcmp(par.nG{g},'Hip'))
       if par.rand_EC_phase
           EC_offset = randn / 3 * par.EC_SD; % random normal sampling at defined SD
           EC_offset = (EC_offset * (1000/par.Hip_freq/2)) / 1000; % set to new SD around width of theta peak
       else
           EC_offset = 0;
       end
       osc.EC  = cos(2*pi*par.Hip_freq*(t+EC_offset+par.Hip_r_phase)); % create EC oscillation               
    end
    % phase reset oscillation to stimulus input
    if (par.([nG{g} '_reset']) && strcmp(par.sim_type,'DL')==1) % isempty(strfind(par.sim_type,'idle'))) 
        t_reset = t(par.pre_stim_length+1:end);
        ph_reset = cos(2*pi*par.([nG{g} '_freq'])*( t_reset+par.reset_ph*250/360/1000 ))*par.([nG{g} '_amp']);
        osc.(nG{g})(par.pre_stim_length+1:end) = ph_reset;
        if(strcmp(par.nG{g},'Hip'))
            ph_reset = cos(2*pi*par.Hip_freq*( t_reset + EC_offset + par.reset_ph*250/360/1000 ));
            osc.EC(par.pre_stim_length+1:end) = ph_reset; 
        end
        % save EC_offset DL
        par.([nG{g} '_EC_offset']) = EC_offset;
    
%     elseif (par.([nG{g} '_reset']) && strcmp(par.sim_type,'NP-AL')==1) % isempty(strfind(par.sim_type,'idle'))) 
%         t_reset = t(par.pre_stim_length+1:end);
%         ph_reset = cos(2*pi*par.([nG{g} '_freq'])*( t_reset+par.stim_G2_ph*250/360/1000 ))*par.([nG{g} '_amp']);
%         osc.(nG{g})(par.pre_stim_length+1:end) = ph_reset;
%     elseif (par.([nG{g} '_reset']) && strcmp(par.sim_type,'P-AL')==1) % isempty(strfind(par.sim_type,'idle'))) 
%         t_reset = t(par.pre_stim_length+1:end);
%         ph_reset = cos(2*pi*par.([nG{g} '_freq'])*( t_reset+par.stim_G1_ph*250/360/1000 ))*par.([nG{g} '_amp']);
%         osc.(nG{g})(par.pre_stim_length+1:end) = ph_reset;
    end
    if(strcmp(par.nG{g},'Hip'))
        osc.EC = 1-(osc.EC/max(osc.EC) + 1)/2; % flip phase of EC by 180
        osc.EC = (osc.EC + (1-par.EC_W)) / ( 1 + (1-par.EC_W) ); % modulate by EC_W variable (lower bound)
    end
    % add oscillation to I input for selected neurons
    I.OSC(c+1:c+par.(['n_' nG{g}]),:) = repmat(osc.(nG{g}), [par.(['n_' nG{g}]) 1]);
    
    %% create ADP function: equation from Jensen/Lisman/Idiart 96
    if(par.([nG{g} '_ADP']))
        adp = 1:1:(1000/par.([nG{g} '_ADP_freq']));
        adp =  adp.*exp(1-adp/(1000/par.([nG{g} '_ADP_freq'])));
        ADP.(nG{g}) = (adp/length(adp)) * par.([nG{g} '_ADP_amp']);
    end
    
    %% create & connect background noise
    % create EPSP for spike events
    psp = compute_psp(par.([nG{g} '_ST']))';
    for i=c+1:c+par.(['n_' nG{g}])
        % create new spike train for each neuron
        spikes = spike_train(par.([nG{g} '_SR']),par.sim_length,par.([nG{g} '_SR'])/100);
        % add EPSP * spike events to I input of selected neurons
        for s = 1:length(spikes)
            I.SYN(i,s:min((s-1)+length(psp),par.sim_length)) = I.SYN(i,s:min((s-1)+length(psp),par.sim_length)) + ...
                spikes(s)*psp(1:min(length(psp),par.sim_length-s+1)) * par.([nG{g} '_SW']);
        end
    end
    
    %% connect network groups together
    for g2 = 1:length(nG)
        c1 = c;
        for i = 1:par.n_Items
            c2 = (g2-1)*par.(['n_' nG{1}]);
            for i2 = 1:par.n_Items
                if(g ~= g2); conn_n = [nG{g2} '_']; else; conn_n = []; end
                if(i == i2); conn_n = [conn_n 'intra']; else; conn_n = [conn_n 'inter']; end
                % find X->Y neuron IDs
                from = c1+1:c1+par.([nG{g} '_per_item']);
                to = c2+1:c2+par.([nG{g2} '_per_item']);
                % set off/on synapses based on X->Y connectivity parameters
                % allow +/-5%  on G1 -> G2 percentage connected
                syn = zeros(length(from), length(to));
                while(abs(sum(syn,'all')/numel(syn) - par.([nG{g} '_' conn_n '_R_conn'])/100) >= 0.05)
                    syn = double(par.([nG{g} '_' conn_n '_R_conn'])/100 > rand(length(from),length(to)));
                end
                % enable STDP on selected synapses
                if( par.([nG{g} '_STDP']) && par.([nG{g2} '_STDP']) ) % STDP enabled?
                    weight_matrix_STDP(from,to) = syn;
                    if( ~strcmp(par.sim_type, 'HL')) % not testing HL fig 7 findings?
                        normal_map = false;
                        str = par.weight_max;
                        if(par.([nG{g} '_' conn_n '_syn_str']) == 0) % low weights
                            % set weights to cluster about the attractor state
                            str = (1-par.T_h) * str + par.W_sd * str * (randn(size(syn))/3);
                        else                                         % high weights
                            % set weights to cluster about 1 minus the attractor state
                            str = par.T_h * str + par.W_sd * str * (randn(size(syn))/3);
                        end
                        str(str > par.weight_max) = par.weight_max;
                    else; normal_map = true;
                    end
                else; normal_map = true;
                end
                if(normal_map)
                    % sets weights as distributed around pre-defined
                    % proportional value of pre-defined weight value
                    str = par.([nG{g} '_' conn_n '_syn_str']);
                    str = par.W_m * str + par.W_sd * str * (randn(size(syn))/3);
                    str(str > par.([nG{g} '_' conn_n '_syn_str'])) = par.([nG{g} '_' conn_n '_syn_str']);
                end
                str(str < 0) = 0;
                
                % input weighted synapses and delays for X->Y
                weight_matrix(from,to,:) = repmat(syn .* str,[1 1 par.sim_length]);
                delay(from,to) = syn * par.([nG{g} '_' conn_n  '_delay']);
                c2 = c2 + par.([nG{g2} '_per_item']); % counting variable
            end
            c1 = c1 + par.([nG{g} '_per_item']); % counting variable
        end
    end
    c = c + par.(['n_' nG{g}]); % counting variable
end
%% no self connections
for i=1:size(weight_matrix,1); for j=1:size(weight_matrix,2); if(i==j); weight_matrix(i,j,:) = 0; weight_matrix_STDP(i,j) = 0; delay(i,j) = 0;end; end; end

%% create & connect stimulus input
if(isempty(strfind(par.sim_type,'idle'))==1)
    if(strfind(par.sim_type,'P')==1) % P stimulation
        n1 = 1 : par.NC_per_item; n2 = [];
    elseif(strfind(par.sim_type,'NP')==1) % NP stimulation
        n1 = par.NC_per_item+1 : par.NC_per_item*2; n2 = [];
    elseif(strfind(par.sim_type,'DL')==1)
        n1 = 1 : par.NC_per_item; % learning (P & NP stimulation)
        n2 = par.NC_per_item+1 : par.NC_per_item*2; % learning (P & NP stimulation)
    else
        n1 = par.n_NC + 1 : par.n_NC + par.Hip_per_item; n2 = [];
    end
    % set oscillatory modulation of inputs
    input_1 = zeros(1, par.sim_length);
    if(par.stim_mod)
        if par.stim_range
            mod_1 = cos(2*pi*par.stim_F*(t+par.stim_G1_ph*(1000/par.stim_F)/360/1000));
            mod_2 = cos(2*pi*par.stim_F*(t+par.stim_G2_ph*(1000/par.stim_F)/360/1000));
        else
            mod_1 = (cos(2*pi*par.stim_F*(t+par.stim_G1_ph*(1000/par.stim_F)/360/1000)) + 1)/2;
            mod_2 = (cos(2*pi*par.stim_F*(t+par.stim_G2_ph*(1000/par.stim_F)/360/1000)) + 1)/2;
        end
    else
        mod_1 = ones(size(input_1)); mod_2 = mod_1;
    end
    
    if(strcmp(par.stim_type,'DC')) %% DC constant (or burst) input
        if(~isempty(par.burst_n)) % for burst input
            theta = (osc.Hip(par.pre_stim_length:end)/max(osc.Hip) + 1)/2;
            stim = [];
            for i = 1:par.burst_n
                stim = [stim ones(1, 500/par.stim_TS) zeros(1,500/par.stim_TS)];
            end
            %stim = ones(1, par.stim_length);
            stim = stim * par.stimulus_strength;
            for b = 1:length(par.stim_phase)
                if(b==1)
                    theta_i(b) = find(theta == par.stim_phase(b), 1, 'first');
                else
                    theta_t = find(theta == par.stim_phase(b));
                    theta_t = theta_t(theta_t > theta_i(b-1));
                    theta_i(b) = theta_t(2);
                end
                input_1(par.pre_stim_length + theta_i(b) - ceil(length(stim)/2) + 1 : ...
                    par.pre_stim_length + theta_i(b) + floor(length(stim)/2)) = stim;
            end
        else % for constant input
            input_1(par.pre_stim_length + 1 : par.pre_stim_length + par.stim_length) = par.stimulus_strength;
        end
        input_2 = input_1 .* mod_2;
        input_1 = input_1 .* mod_1;
        
        I.DC(n1,:) = repmat(input_1, [length(n1),1]);
        I.DC(n2,:) = repmat(input_2, [length(n2),1]);
    else
        T = 1:1:par.stim_length;
        % create stimulus alpha function
        stim = [zeros(1, par.pre_stim_length) (exp(1)*T./par.stim_TS).*exp(-T./par.stim_TS).*heaviside(T)];
        psp = compute_psp(par.stimulus_ST)';
        
        % create spike train based on parameters
        input_1(par.pre_stim_length+1:end) = ...
            spike_train(par.stimulus_SR, par.stim_length, ceil(par.stimulus_SR/100));
        input_1 = input_1 .* stim * par.stimulus_SW;  % modulate spike train by alpha function & weight
%         input_1 = input_1 * par.stimulus_SW;
        input_2 = input_1 .* mod_2;                   % modulate stim 2 by oscillation (if selected)
        input_1 = input_1 .* mod_1;                   % modulate stim 1 by oscillation (if selected)
        
        % add generated spike events to future input
        for s = 1:length(input_1)
            I.SYN(n1,s:min((s-1)+length(psp),par.sim_length)) = I.SYN(n1,s:min((s-1)+length(psp),par.sim_length)) + ...
                repmat(input_1(s)*psp(1:min(length(psp),par.sim_length-s+1)),[length(n1) 1]);
            I.SYN(n2,s:min((s-1)+length(psp),par.sim_length)) = I.SYN(n2,s:min((s-1)+length(psp),par.sim_length)) + ...
                repmat(input_2(s)*psp(1:min(length(psp),par.sim_length-s+1)),[length(n2) 1]);
        end
    end
    
end

end

