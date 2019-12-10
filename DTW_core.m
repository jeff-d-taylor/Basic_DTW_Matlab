% Simplified, modular DTW function for calculating the similarity between
% reference and query datasets

% inputs: 
%    reference_timestamps: column of timestamps for reference timeseries data 
%    query_timestamps: column of timestamps for query timeseries data
%    reference_comparison_data: column of reference data for calculating similarity
%    query_comparison_data: column of query data for calculating similarity

% outputs:
%   path: warp path, which provides the estimated match points between reference and query datasets
%   cost and cumulative cost matrices

% Step 1: Calculate cost matrix
% Step 2: Calculate cumulative cost matrix
% Step 3: Produce the warp path


function [path,cost,cumulative_cost] = DTW_core(reference_timestamps,query_timestamps,reference_comparison_data,query_comparison_data)

% Step 1: Calculate cost matrix (using time and distance constraints)
    % cost_matrix_constain function inputs: (reference timestamps, query timestamps, reference comparison values, query comparison values)
    cost = cost_matrix_calculator(reference_timestamps,query_timestamps,reference_comparison_data,query_comparison_data);
    % output is the cost matrix

% calculate the cumulative cost matrix
cumulative_cost = cumulative_cost_calculator(cost);

% find the shortest path
path = shortest_path_calculator(cumulative_cost);

% internal functions for DTW:

    function cost = cost_matrix_calculator(reference_time,query_time,reference_values,query_values)
        
        % This version of the cost function calculates the cost matrix for matching
        % all points in one time series to all points in another time series.
        
        % input: reference time series timestamps, query time series timestamps
        %        reference_values = column of reference data for calculating similarity
        %        query_values = column of query data for calculating similarity
        
        % output: cost matrix, of size [size(reference_time,1) X size(query_time,1)]
        
        % preallocation to improve speed
        cost = zeros(size(reference_time,1),size(query_time,1));
        
        for i = 1:size(reference_time,1) % for each time step in reference timeseries
            for j = 1:size(query_time,1) % for each time step in query timeseries
                    cost(i,j) = abs(reference_values(i,1) - query_values(j,1));
            end
        end
    end

    function cumulative_cost = cumulative_cost_calculator(cost_matrix)
        % This function calculates the cumulative cost matrix. It uses a specific search pattern to do so.
        
        % input: cost_matrix = cost of matching two timeseries points in their respective datasets
        
        % temporary variables:
        %     minimum_step = the least-cost required to arrive at the previous step before the current location in the matrix
        
        % output:
        %     cumulative_cost = matrix containing the cumulative cost of matching all points in one timeseries dataset
        %                       to all of the other points in the other timeseries dataset
        
        cumulative_cost = zeros(size(cost_matrix,1),size(cost_matrix,2)); % pre-allocate to improve speed
        
        % at the very beginning (1,1), the initial starting cost is zero
        for j = 2:size(cost_matrix,2) % working with the first row of the matrix
            cumulative_cost(1,j) = cumulative_cost(1,j-1) + cost_matrix(1,j);
        end
        for i = 2:size(cost_matrix,1) % working with the first column of the matrix
            cumulative_cost(i,1) = cumulative_cost(i-1,1)+cost_matrix(i,1);
        end
        
        for i = 2:size(cost_matrix,1) % for each reference timestamp
            for j = 2:size(cost_matrix,2) % for each query timestamp
                % calculate the minimum cost amongst three options for movement in the matrix (back within this column, back within this row, and back diagonally
                % use the least cost location to calculate the least cost for arriving at the current location in the matrix
                cumulative_cost(i,j) = cost_matrix(i,j) + min(min(cumulative_cost(i-1,j),cumulative_cost(i,j-1)),cumulative_cost(i-1,j-1));
            end
        end
    end

    function [path_sort] = shortest_path_calculator(cumulative_cost)
        
        % This function finds the shortest path through the cumulative cost matrix.
        
        % input: cumulative_cost = cumulative cost matrix
        
        % temporary variables: i = reference index, j = query index, k = path index
        %                      path = structure of reference and query matching, from end to beginning
        
        % output: path_sort = structure containing the reference(x) and query(y) incidces for optimal matching, from beginning to end
        
        i = size(cumulative_cost,1); % number of timeseries points in reference timeseries
        % i is also used to keep track of the current location in the matrix
        j = size(cumulative_cost,2); % number of timeseries points in query timeseries
        % j is also used to keep track of the current location in the matrix
        k = 1; % initialize the first step in the shortest path
        % start building the shortest path from the end of the matrix, starting
        % with the last points in the matrix
        path.x(k,1) = size(cumulative_cost,1); % i and x = reference
        path.y(k,1) = size(cumulative_cost,2); % j and y = query
        
        % In the cumulative cost matrix(i,j), the reference location stays the same
        % within the row, and the query location stays the same within the column
        
        % Find shortest path
        while i>=1 && j>=1 % as long as we don't reach the beginning of the matrix
            k = k+1; % move forward a step in the path
            
            if i == 1 && j == 1 % if we reach the beginning of the matrix
                break % stop moving through the matrix
            elseif i == 1 % if we reach the first row, stay in the same row
                path.x(k,1) = i; path.y(k,1) = j-1; j = j-1; % only increment to next column
                continue
            elseif j == 1 % if we reach the first column, stay in the same column
                path.x(k,1) = i-1; path.y(k,1) = j; i = i-1; % only increment to next row
                continue
            end
            
            if i>1 && j>1 % if we haven't reached the end of the matrix yet
                % calculate the next step cost - the minimum cost to reach the next
                % point in the matrix, moving vertically, horizontally, or diagonally
                next_step_cost = min(min(cumulative_cost(i-1,j),cumulative_cost(i,j-1)),cumulative_cost(i-1,j-1));
                
                % check which step to take next, store the next step in the
                % shortest path, and increment to the next step location in i & j
                if next_step_cost == cumulative_cost(i-1,j-1)
                    % if all steps cost the same, default step is diagonally
                    path.x(k,1) = i-1; path.y(k,1) = j-1; i = i-1; j = j-1;
                elseif next_step_cost == cumulative_cost(i,j-1)
                    path.x(k,1) = i; path.y(k,1) = j-1; j = j-1;
                else
                    path.x(k,1) = i-1; path.y(k,1) = j; i = i-1;
                end
            end
        end
        % re-sort the datasets to set the path from beginning to end
        path_sort.x = sort(path.x);  path_sort.y = sort(path.y);
    end

end


