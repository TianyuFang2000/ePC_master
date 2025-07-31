function permuted_edge_labels = permute_edge_labels(edge_labels, Nv)
    % permute_edge_labels reshapes the edge labels into the upper triangle of
    % an Nv x Nv matrix, permutes the rows/columns, and extracts the permuted
    % edge labels from the upper triangle.
    %
    % INPUT:
    % edge_labels - A vector of length Ne = Nv * (Nv - 1) / 2, representing edge labels.
    % Nv          - The number of nodes (vertices) in the network.
    %
    % OUTPUT:
    % permuted_edge_labels - A vector containing the permuted edge labels from the upper triangle.

    % Number of edges
    Ne = length(edge_labels);
    
    % Number of nodes
    % Nv is given as input, ensure consistency with Ne
    expected_Ne = Nv * (Nv - 1) / 2;
    if Ne ~= expected_Ne
        error('The length of edge_labels does not match the expected number of edges for Nv nodes.');
    end
    
    % Step 1: Reshape edge labels into the upper triangle of a matrix
    % Initialize the adjacency matrix (Nv x Nv) as zeros
    adj_matrix = zeros(Nv);
    
    % Fill the upper triangular part with the edge labels
    mask=triu(true(Nv),1);
    adj_matrix(mask) = edge_labels;
    adj_matrix = adj_matrix+adj_matrix';   
    
    % Step 2: Randomly permute the rows and columns of the matrix
    perm_order = randperm(Nv);  % Generate a random permutation of row/column indices
    
    % Apply the permutation to the matrix
    permuted_matrix = adj_matrix(perm_order, perm_order);
    
    % Step 3: Extract the permuted edge labels from the upper triangle
    permuted_edge_labels = permuted_matrix(mask);
end
