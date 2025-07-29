function plot_gene_list(weight)

figure;
% Original vector X
X = weight;  

% Desired length of the resampled vector
desiredLength = 150;

% Original indices based on the length of X
originalIndices = linspace(1, length(X), length(X));  

% New indices for the desired resampled vector length
newIndices = linspace(1, length(X), desiredLength);

% Perform linear interpolation to resample the vector
X_resampled = interp1(originalIndices, X, newIndices, 'linear');

% X_resampled = [Weighted(1:5); Weighted(end-4:end) ];
hBar = bar(flip(X_resampled),'EdgeColor','none');  % Create the bar chart
numBars = length(X_resampled);  % Number of bars
blueToLightBlueLength = length(find(X_resampled<0));  % Length of the blue-to-light-blue transition
lightRedToRedLength = length(find(X_resampled>0));   % Length of the light-red-to-red transition

% Define specific RGB values for each color in the gradient
blue = [0 0.3 0.6];        % Blue
lightBlue = [0.9 0.95 0.98];  % Light Blue
lightRed = [1, 0.9, 0.9];  % Light Red
red = [1, 0, 0];        % Red

% Create the first gradient: blue to light blue
blueToLightBlue = interp1([0, 1], [blue; lightBlue], linspace(0, 1, blueToLightBlueLength), 'linear');

% Create the second gradient: light red to red
lightRedToRed = interp1([0, 1], [lightRed; red], linspace(0, 1, lightRedToRedLength), 'linear');

% Combine both gradients to form the full colormap
customColormap = [blueToLightBlue; lightRedToRed];

% Color the bars using the colormap
colors = flip(customColormap);  % Get the colormap as an N x 3 matrix
hBar.FaceColor = 'flat';  % Allow individual coloring of bars
hBar.CData = colors(1:numBars, :);  % Assign a color to each bar

% Set Y-axis labels
set(gca, 'LineWidth',1.25,'TickDir','out','box','off','XTick',[],'FontSize',18);
% Add labels
xlabel('Genes');
ylabel('Pearsons r');