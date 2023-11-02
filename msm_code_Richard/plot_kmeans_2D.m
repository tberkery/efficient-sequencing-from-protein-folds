%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot clusters and means of k-means clustering using data, cluster indice.
% 
% 
% Parameters:
%   clust_data - cluster data, 1 point per row
%   dim_labels - dimension labels for clust_data (e.g: ['x'; 'y'])
%   clust_idx - k-means clusters indexed 1 to k, 1 per point in clust_data
% 
% Return:
%   Plot of k-means clusters, their means. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=plot_kmeans_2D(clust_data, dim_labels, dim_targets, clust_idx)
    [num_points, num_dimensions] = size(clust_data);
    k = numel(dim_labels);
    unique_clusters = 1:k;
    
    % loosey goosey error handling
    if (num_dimensions ~= numel(dim_labels))
        error('Dimension labels do not match clust_data (number of columns)');
    end
    if (numel(dim_labels) ~= numel(unique(clust_idx)))
        error('Dimension labels do not match cluster classifications (clust_idx)');
    end
    if (num_points ~= numel(clust_idx))
        error('Clusters (clust_idx) do not represent cluster data (clust_data)');
    end
    if (numel(dim_targets) ~= 2)
        error('Must specify exactly two dimensions for a 2D plot')
    end
    if (~all(unique_clusters == transpose(sort(unique(clust_idx)))))
        error('Cluster indices must be [1, 2, ..., k]');
    end
    for dim_target = dim_targets
        if (dim_target < 1 || dim_target > k)
            error('Cluster targets must be within [1, k]');
        end
    end
    
    hold on
    all_means = [];
    for i = unique_clusters
        % row index of clust_data AND clust_idx, matching cluster i
        curr_idx = (clust_idx == i);
        % filter for current cluster
        curr_clust_data = clust_data(curr_idx, :);
        cluster_label = strcat('Cluster: ', num2str(i));
        
        % save mean of cluster data (per dimension)
        curr_means = mean(curr_clust_data, 1);
        all_means = [all_means; curr_means];
        
        % plot cluster data
        scatter(curr_clust_data(:, dim_targets(1)), curr_clust_data(:, dim_targets(2)), ...
            10, 'DisplayName', cluster_label);
    end
    scatter(all_means(:, dim_targets(1)), all_means(:, dim_targets(2)), ...
        50, 'black', 'filled', 'DisplayName', 'Means');
    legend;
    xlabel(dim_labels(1));
    ylabel(dim_labels(2));
    title(strcat('K-Means Clusters by: (', dim_labels(1), ',', dim_labels(2), ')'));
    hold off;
end