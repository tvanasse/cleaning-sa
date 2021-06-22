function clusters = snpm_find_clusters_graphalgs(t_values, threshold, sparse_adjacency_matrix)
    if threshold > 0
        significant_vertices_logical = t_values > threshold;
    else
        significant_vertices_logical = t_values < threshold;
    end
    
    significant_vertices                = find(significant_vertices_logical);
    significant_sparse_adjacency_matrix = sparse_adjacency_matrix(significant_vertices_logical, significant_vertices_logical);
    
    % strongly connected component - every vertexd is reachable from every other vertex
    [num_clusters cluster_nums] = graphalgs('scc', 0, 0, significant_sparse_adjacency_matrix);
    
    clusters = cell(1, num_clusters);
    
    for cluster_num = 1:num_clusters
        clusters{cluster_num} = significant_vertices(cluster_nums == cluster_num);
    end
end
