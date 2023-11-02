% Richard Hu
% rhu16
% AS.250.302
clear; close all; clc
rng(1);

%% 1
rate_vs_conc = load('rate_vs_conc.dat');

substrateConc = rate_vs_conc(:,1);
reactionRate = rate_vs_conc(:,2);

beta0 = [2 1];
[beta, R, J, CovB, MSE, ErrorModelInfo] = nlinfit(substrateConc, reactionRate, @MM, beta0);
beta_variances = [CovB(1,1) CovB(2,2)];

reactionRate_model = MM(beta, substrateConc);

hold on
scatter(substrateConc, reactionRate, 'filled', 'DisplayName', 'data')
plot(substrateConc, reactionRate_model, 'DisplayName', 'MM fit')
legend
xlabel('Substrate Concentration')
ylabel('Reaction Rate')
title('1c) Michaelis-Menten Reaction Rate Model')
hold off
saveas(gcf, "problem1.eps", "psc2")

fprintf("1d)variances\n");
fprintf("\tV_max: %d,\tK_m: %d\n", beta_variances(1), beta_variances(2));

%% 2
fprintf("\n")
close all;

clust_data = load("clust_data.dat");
dim_labels = ['x'; 'y'];
dim_targets = [1 2];

clust_idx = kmeans(clust_data, 2);

plot_kmeans_2D(clust_data, dim_labels, dim_targets, clust_idx)
saveas(gcf, "problem2.eps", "psc2")

%% 3
fprintf("\n")
close all;
Xvstime = load('Xvstime.dat');
Yvstime = load('Yvstime.dat');
Zvstime = load('Zvstime.dat');

bead_pairs = [1 10; 1 4; 1 5; 2 6; 4 7; 5 10; 5 9];
num_clusters = 6;
[transition_matrix_smallstep, transition_probabilities_smallstep, state_probabilities_smallstep, features_cluster_means_smallstep] = ...
    markovState_kmeans(Xvstime, Yvstime, Zvstime, bead_pairs, num_clusters, 50);
fprintf("3c)\n");
display(transition_probabilities_smallstep)
fprintf("The highest probabilities are the diagonal, indicating that\n");
fprintf("the state ends up in the same cluster. The barriers for transition,\n");
fprintf("may be relatively high---or the states just too changed slowly.\n");
fprintf("\n");
fprintf("The lowest probabilities on my seed and kmeans were (1->3), (1->4)\n,");
fprintf("(3->1), (4->1), (4->6), and (6->4).\n");
fprintf("For a pair of clusters with low probabilities in both directions\n");
fprintf("(all 3 pairs), their conformations are distant. Tentative corroboration\n");
fprintf("by seeing the distance between clusters:\n");
display(features_cluster_means_smallstep)
fprintf("(1,4) is a distant pair cluster-mean-wise:\n");
display(distance(features_cluster_means_smallstep(1,:), features_cluster_means_smallstep(4,:)))
fprintf("(1,5) is more typical:\n");
display(distance(features_cluster_means_smallstep(1,:), features_cluster_means_smallstep(5,:)))

fprintf("3d) Probabilities of each state\n");
for i = 1:numel(state_probabilities_smallstep)
    fprintf("\tCluster %d: %f\n", i, state_probabilities_smallstep(i));
end

[transition_matrix_bigstep, transition_probabilities_bigstep, state_probabilities_bigstep, features_cluster_means_bigstep] = ...
    markovState_kmeans(Xvstime, Yvstime, Zvstime, bead_pairs, num_clusters, 100);
display(transition_probabilities_bigstep)

fprintf("Perfect markovian: ptrans_100 - ptrans_50*ptrans_50 = 0\n")
fprintf("Fractional error from timestep of 100:\n")
display((transition_probabilities_bigstep - transition_probabilities_smallstep*transition_probabilities_smallstep)/transition_probabilities_bigstep)
fprintf("Close but not exact\n")

%% 4
fprintf("\n")
sequence_table = readtable("onesequence_-22.79.dat");
hpChain = sequence_table.Var1;
chain_coords = sequence_table{:, 2:4};
step_size = 0.05;
step_factor = 0.999;
tolerance = 1e-6;
[coords_minE, F_minE, E_minE] = ...
    steepest_descent_HP(hpChain, chain_coords, step_size, step_factor, tolerance);
fprintf("4b) Convergence criteria: lambda=%f, a=%f, tolerance=%f\n",...
    step_size, step_factor, tolerance);
fprintf("Note: Using the other parameters, a relatively large variation of lambdas\n");
fprintf("still converged to the same minimum (tried within 0.1, 5)\n");
fprintf("4c) Energy at minima:%f\n", E_minE);

plot3(coords_minE(:,1), coords_minE(:,2), coords_minE(:,3))
title("4c) Configuration at energy minima")
xlabel("x")
ylabel("y")
zlabel("z")

saveas(gcf, "problem4.eps", "psc2")
%%
close all;