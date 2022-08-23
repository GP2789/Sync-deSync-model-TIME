function [ ] = make_figures( dir_name )
%% DESCRIPTION
% Simply run this function specificying a directory name, where all the
% data and figures will be reproduced. The functions below run code for
% reproducing the Huerta & Lisman (1995) findings, firstly with a further
% hetero-synaptic plasticity function and secondly without it.
% Hetero-synaptic plasticity is not essential to reproducing the main
% findings, but might be useful to those who wish to implement this
% learning rule on a larger data-set in order to maintain network
% equilibrium in a bio-pyhsically plausible way. Each set of simulations
% first runs a "single burst," whereby the hippocampal populations are
% stimulated with a burst containing a variable number of spikes (1-4), either at
% the peak or trough of the ongoing theta oscillation. Then we run a
% "double burst," whereby a burst of 4 spikes occurs at the trough and then
% the peak of the ongoing theta oscillation. These simulations assess how
% plasticity is designed to be dependent on the ongoing theta oscillation in this model -
% replicating the procedure used in Huerta & Lisman (1995) which first observed this phenomenon. 

%% FUNCTIONS TO REPRODUCE SIMULATIONS
tic
% Here a set of simu
fprintf('\nsimulations with hetero-synaptic plasticity:\n');
fprintf('simulating a single burst ... \n');
analyse_theta( [dir_name '/het-STDP_'], 25, true, 'single', true);
fprintf('simulating a double burst ... \n');
analyse_theta( [dir_name '/het-STDP_'], 100,true, 'both', true);

fprintf('\nsimulations without hetero-synaptic plasticity:\n');
fprintf('simulating a single burst ... \n');
analyse_theta( [dir_name '/'], 25, true, 'single', false);
fprintf('simulating a double burst ... \n');
analyse_theta( [dir_name '/'], 100, true, 'both', false);

S = toc;
M = floor(S/60);
S = S - M*60;
fprintf('\ncompleted in %.0f mins %.0f secs\n\n',M,S);

end

