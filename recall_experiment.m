function [  ] = recall_experiment( filename )
% Perform learning and recall experiment by creating and simulating
% networks specified by set_parameters. Give a filename for directory location (string).

%% DECLARATIONS
% initiate global variables and set parameters
global weight_matrix weight_matrix_STDP par;
%clearvars -global par; 
%set_parameters; global par; 

% make directories for simulation
if(exist(filename,'dir')~=7); mkdir(filename); end
if(exist([filename '/Data'],'dir')~=7); mkdir([filename '/Data']); end

%% MAIN SIMULATION

% random phase for static oscillations
par.Hip_r_phase_n = zeros(1, length(par.sim_order));
par.NC_r_phase_n = zeros(1, length(par.sim_order));

if( par.rand_phase )
    deg = 0:30:360; % 0-360 in 30deg intervals
    par.Hip_r_phase = round(deg(ceil(rand()*length(deg)))*(250/360))/1000; % choose random Theta phase
    par.NC_r_phase = round(deg(ceil(rand()*length(deg)))*(100/360))/1000; % choose random Alpha phase
else % set no random phases for static oscillations
    par.Hip_r_phase = 0; par.NC_r_phase = 0;
end
par.recall = false;
for n = 1:length(par.sim_order) % loop through each experiment phase
    fprintf(['\nsimulating ' par.sim_order{n} '...\n' ])
    % set sim_length and sim_type
    if(strcmp(par.sim_order{n},'idling')==1)
        par.sim_length = par.idling_length;
    else
        par.sim_length = par.pre_stim_length + par.stim_length;
    end
    par.sim_type = par.sim_order_n{n};
    
    % no STDP learning after learning
    if(~isempty(strfind(par.sim_order_n{n},'AL'))==1)
        par.Hip_STDP = false;
        par.recall = true;
        % create network
        create_network(); 
        par.Hip_r_phase_n(n) = par.Hip_r_phase; par.NC_r_phase_n(n) = par.NC_r_phase;
        if(n>1) % carry forward weight matrix from previous phase but STDP matrix is not carried from previous phase if no STDP (= 0)
            weight_matrix = repmat(carried_wm, [1 1 par.sim_length]); 
        end
    else
%     end
        par.Hip_STDP = true;
        par.recall = false;
%     end
        % create network
        create_network(); 
        par.Hip_r_phase_n(n) = par.Hip_r_phase; par.NC_r_phase_n(n) = par.NC_r_phase; 
        if(n>1) % carry forward weight matrix from previous phase
            weight_matrix = repmat(carried_wm, [1 1 par.sim_length]); 
            weight_matrix_STDP = carried_wm_STDP;
        else % store weight matrix from first phase
            carried_wm_STDP = weight_matrix_STDP;
        end
    end
    [data.sim_stats] = simulate_network(); % simulate network
    sim_length_total(n) = par.sim_length;
    carried_wm = weight_matrix(:,:,end); 
    save([filename '/Data/' par.sim_order_n{n} '.mat'],'data')
    
end
par.sim_length_total = sim_length_total;
save([filename '/Data/parameters.mat'],'par')

end


