% This function plots the match solution from DTW warp path.

% input: reference_time = column of timestamp data for leader vehicle
%        query_time = column of timestamp data for follower vehicle
%        reference_values = column of data to be plotted for reference dataset
%        query_values = column of data to be plotted for query dataset
%        warp path = structure (x and y) which identifies the match points between the two timeseries datasets
%        cost = cost of matching two timeseries points in their respective datasets

function DTW_matching_plot(reference_time,query_time,reference_values,query_values,path)

plot(query_time,query_values,'k'); hold on % plot query timeseries data
plot(reference_time,reference_values,'k'); hold on % plot reference timeseries data

% plot best matching points
% x = reference data, y = query data
line_color = 'k'; % default line color is black
for i = 1:size(path.x,1) % for each match point in the match solution
    % only plot for values of which we are relatively sure of their value
    plot([reference_time(path.x(i,1),1),query_time(path.y(i,1),1)],[reference_values(path.x(i,1),1),query_values(path.y(i,1),1)],line_color);
    hold on
end
