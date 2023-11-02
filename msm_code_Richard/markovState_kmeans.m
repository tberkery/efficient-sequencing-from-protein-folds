function[transition_matrix, transition_probabilities, state_probabilities, features_cluster_means]=...
    markovState_kmeans(Xvstime, Yvstime, Zvstime, bead_pairs, num_clusters, timestep)
    [num_features, group_size] = size(bead_pairs);
    if (group_size ~= 2)
        error("Error: Invalid pair format (should follow [x1 x2; y1 y2; ...])");
    end
    [num_beads, duration] = size(Xvstime);
    
    % find feature distanes for all timepoints
    features_vs_time = zeros(num_features, duration);
    for time = 1:duration
        for pair_idx = 1:num_features
            feature_pair = bead_pairs(pair_idx, :);
            
            % get bead coordinates
            bead1 = [Xvstime(feature_pair(1), time), ...
                Yvstime(feature_pair(1), time), ...
                Zvstime(feature_pair(1), time)];
            bead2 = [Xvstime(feature_pair(2), time), ...
                Yvstime(feature_pair(2), time), ...
                Zvstime(feature_pair(2), time)];
            % get scalar distance of the bead pair
            features_vs_time(pair_idx, time) = distance(bead1, ...
                bead2);
        end
    end
    
    % cluster
    [features_cluster_idx, features_cluster_means] = kmeans(transpose(features_vs_time), num_clusters, 'MaxIter', 1000);
    
    % make transition matrix
    transition_matrix = zeros(num_clusters);
    % loop (t1, t2=t1+duration) via t_i=1:duration
    for t_ref = 1:(duration-timestep)
        % t_ref: time of reference state; t_new: time of transitioned state
        t_new = t_ref + timestep;
        transition_matrix(features_cluster_idx(t_ref), features_cluster_idx(t_new)) = ...
            transition_matrix(features_cluster_idx(t_ref), features_cluster_idx(t_new)) + 1;
    end
    state_probabilities = probability_from_transitions(transition_matrix);
    transition_probabilities = norm_rows(transition_matrix);
end