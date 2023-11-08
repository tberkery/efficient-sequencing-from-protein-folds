% assign a data point to the nearest cluster, using its distance features
% as input
%
% Parameters:
%   feature_i - a 1xn feature/distance vector, to be assigned to a cluster
%   clust_means - the mean feature/distance vectors (each 1xn) of k
%   clusters, stored in dimensions kxn
%
% Return:
%   clust_idx - the index of the closest cluster
%   clust_dist - feature/distance vector describing difference between
%   clust_idx and the closest mean: feature_i - clust_means(clust_idx,:)
function[clust_idx, clust_diff] = classify_feature(feature_x, clust_means)
    [num_x, dim_x] = size(feature_x);
    [k, dim_means] = size(clust_means);
    if (num_x ~= 1)
        error("feature_x should be a 1xn vector describing one feature");
    end
    if (dim_x ~= dim_means)
        error("Features provided between feature_i and clust_means are incompatible. Check dimensions");
    end

    % find min cluster
    clust_dist = Inf;
    for (idx = 1:k)
        dist_i = distance_features(feature_x, clust_means(idx, :));
        if (dist_i < clust_dist)
            clust_dist = dist_i;
            clust_diff = clust_means(idx, :);
            clust_idx = idx;
        end
    end
end