function [ ] = STDP( spikes, t, theta )
%% INITIALISE
global weight_matrix_STDP I_REC par weight_matrix wc_pos wc_neg;

weight_matrix(:, :, t) = weight_matrix(:, :, t-1);         % update weight matrix
p_w = weight_matrix(:, :, t) / par.weight_max;             % normalise weight matrix (0 <= p_w <= 1)

%% CALCIUM DECAY
wc_pos(:, :, t) = wc_pos(:, :, t-1) ...                    % potentially contributing synapses
    - (wc_pos(:, :, t-1)./par.T_Ca);                       % decay at a certain exponential rate

wc_neg(:, :, t) = wc_neg(:, :, t-1) ...                    % potentially competing synapses
    - (wc_neg(:, :, t-1)./par.T_Ca);                       % decay at a certain exponential rate

%% SPIKE INDUCED PLASTICITY
if( ~isempty(spikes) )
    %% PRE-POST SYNAPTIC CALCIUM CHANGES MODULATED BY THETA
    wc_pos(spikes, :, t) = (wc_pos(spikes, :, t) ...        % potentially contibuting synapses
        + par.a_pos * theta) ...                            % addition is modulated by theta
        .* weight_matrix_STDP(spikes, :);                   % select only plastic synapses
    
%     wc_neg(:, spikes, t) = (wc_neg(:, spikes, t) ...    % potentially competing synapses
%             + par.a_neg) ...                        % addition is modulated by theta
%             .* weight_matrix_STDP(:, spikes);  
%         
%     wc_pos(:, spikes, t) = (wc_pos(:, spikes, t) ...        % potentially contibuting synapses
%         + par.a_pos) ...                            % addition is modulated by theta
%         .* weight_matrix_STDP(:, spikes);                   % select only plastic synapses
    
    if( par.theta_LTD )                                     % LTD MODULATED BY -THETA (SdS V2)
        wc_neg(:, spikes, t) = (wc_neg(:, spikes, t) ...    % potentially competing synapses
            + par.a_neg * (1-theta)) ...                    % addition is modulated by -theta
            .* weight_matrix_STDP(:, spikes);               % select only plastic synapses
    else                                                    % ALL STDP MODULATED BY +THETA (ORIGINAL SdS)
        wc_neg(:, spikes, t) = (wc_neg(:, spikes, t) ...    % potentially competing synapses
            + par.a_neg * theta) ...                        % addition is modulated by theta
            .* weight_matrix_STDP(:, spikes);               % select only plastic synapses
    end
    
    %% SPIKE-TIMING-DEPENDENT-PLASTICITY (STDP)
    if( par.STDP_th ) % FOR THRESHOLD INDUCED PLASTICITY
        %% HOMO-SYNAPTIC LTP (THRESHOLD INDUCED)
        LTP_criteria = heaviside( wc_pos(:, spikes, t) - par.T_p );
        if( any(LTP_criteria > 0, 'all') )
            dW = (par.G_p .* ( 1-p_w(:, spikes) ) ...               % rate of change (capped for already high weights)
                .* weight_matrix_STDP(:, spikes) ...                % select plastic synapses only
                .*  LTP_criteria )...                               % where calcium > threshold
                .* ( wc_pos(:, spikes ,t) - par.T_p );              % modulated by amount over threshold

            % according to a time constant
            p_w(:, spikes) = p_w(:, spikes) + dW;                   % add weight change
            LTP = dW;
            I_REC.LTP(:, spikes, t) = dW * par.weight_max;          % for recording purposes
            clear dW;
        end

        %% HOMO-SYNAPTIC LTD (THRESHOLD INDUCED)
        LTD_criteria = heaviside( wc_neg(spikes, :, t) - par.T_d );
        if( any(LTD_criteria > 0,'all') )
            dW = (par.G_p/2 .* p_w(spikes, :) ...                    % rate of change (capped for already low weights)
                .* weight_matrix_STDP(spikes, :) ...                 % select plastic synapses only
                .*  LTD_criteria ) ...                               % where calcium > threshold
                .* ( wc_neg(spikes, :, t) - par.T_d );               % modulated by amount over threshold
            % according to a time constant
            p_w(spikes, :) = p_w(spikes, :) - dW;                    % subtract weight change
            I_REC.LTD(spikes, :, t) = - dW * par.weight_max;         % for recording purposes
            clear dW;
        end
        
        %% HETERO-SYNAPTIC PLASTICITY (THRESHOLD INDUCED)
        if( par.STDP_hetero && any(LTP_criteria > 0, 'all') )
            pW = (exp(1) * (p_w(:, spikes) - (1-par.T_h)).^2 + 0.1) ...  % calculate probability of weight change
                <= rand(size(p_w(:, spikes)));                           % select random synapses to undergo changes
            pW = pW .* weight_matrix_STDP(:, spikes) ...                 % select only plastic synapses
                .* any(LTP_criteria, 1)...                               % for all other post-synaptic synapses 
                .* heaviside( par.T_p - wc_pos(:, spikes, t) );          % not including the currently potentiating synapses

            dW = ((1 ./ (1 + exp( (p_w(:, spikes) - 0.5)*exp(1) ))) ...  % calculate synaptic change
                - par.T_h) ...                                           % gravitating towards an attractor state (0<->1)
                * par.G_h .* pW;                                         % at a rate of change on selected synapses
                
            dW = dW .* (sum(LTP,1) ./ -sum(dW,1)) * par.P_h;              % modulated by amount of post-synaptic LTP
            dW(isnan(dW)) = 0;
            
            %dN = ((randn(size(p_w))*3) * 0.001 ) .* pW;                 % add noise to synaptic changes
            
            p_w(:,spikes) = p_w(:,spikes) + dW;% - dN;                   % add weight change
            I_REC.HET(:,spikes,t) = dW * par.weight_max;                 % for recording purposes
            clear dW dN LTP;
        end
        clear LTP_criteria LTD_criteria;
    else % FOR THE ORIGINAL STDP RULE
        %% HOMO-SYNAPTIC LTP (ORIGINAL)
        dW = wc_pos(:, spikes, t) .* weight_matrix_STDP(:, spikes);  % simple addition of potential calcium
        p_w(:, spikes) = p_w(:, spikes) + dW;                        % add weight change
        I_REC.LTP(:, spikes, t) = dW * par.weight_max;               % for recording purposes
        clear dW;
        
        %% HOMO-SYNAPTIC LTD (ORIGINAL)
        dW = wc_neg(spikes, :, t) .* weight_matrix_STDP(spikes, :);  % simple addition of potential calcium
%         dW = wc_pos(spikes, :, t) .* weight_matrix_STDP(spikes, :);  % simple addition of potential calcium
%         p_w(spikes, :) = p_w(spikes, :) + dW;                        % subtract weight change
        p_w(spikes, :) = p_w(spikes, :) - dW;                        % subtract weight change
        I_REC.LTD(spikes, :, t) = - dW * par.weight_max;               % for recording purposes
        clear dW;
    end
end

%% NON-SPECIFIC THETA-MODULATED WEIGHT DECAY
if( par.weight_decay )
    dW = exppdf(p_w, par.T_weight_decay)*0.00007 ...    % rate of decay ( higher for smaller weights )
        .* weight_matrix_STDP ...                       % select plastic synapses only
        .* (1-theta);                                   % modulate by reverse phase of theta
    p_w = p_w - dW;                                     % subtract weight change
    I_REC.STDP_decay(:, :, t) = dW;                     % for recording purposes
    clear dW;
end

%% piecewise linear bounding function (no weights below 0 or above max weight)
p_w( p_w < 0 ) = 0;
p_w( p_w > 1 ) = 1;

%% reassign weights
weight_matrix(:, :, t) = p_w * par.weight_max;
clear p_w;

end

