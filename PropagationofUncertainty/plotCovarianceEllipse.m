function [Data1, Data2] = plotCovarianceEllipse(data)
    % plotCovarianceEllipse - Plots the covariance ellipse for 2D data.
    %
    % Usage:
    %   plotCovarianceEllipse(data)
    %
    % Input:
    %   data - A matrix with two columns representing 2D data points, 
    %          where each row is a data point [x, y].

    % Calculate the covariance matrix and mean of the data
    cov_matrix = cov(data);
    mean_x = mean(data(:, 1));
    mean_y = mean(data(:, 2));

    % Eigenvalues and eigenvectors of the covariance matrix
    [eigvec, eigval] = eig(cov_matrix);

    % Sort eigenvalues and eigenvectors in descending order
    [eigval_sorted, idx] = sort(diag(eigval), 'descend');
    eigvec = eigvec(:, idx);

    % Calculate ellipse angle, width, and height
    angle  = atan2(eigvec(2,1), eigvec(1,1));
    width  = 2 * sqrt(eigval_sorted(1));  % Major axis (2*std deviation)
    height = 2 * sqrt(eigval_sorted(2)); % Minor axis (2*std deviation)

    % % Plot data points
    % figure;
    % scatter(data(:, 1), data(:, 2), 'b', 'filled');
    % hold on;

    % Plot the covariance ellipse
    t = linspace(0, 2 * pi, 100);
    ellipse_x = width/2 * cos(t);
    ellipse_y = height/2 * sin(t);
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    ellipse_coords = R * [ellipse_x; ellipse_y];
    Data1 = mean_x + ellipse_coords(1,:);
    Data2 = mean_y + ellipse_coords(2,:);
    % plot(mean_x + ellipse_coords(1,:), mean_y + ellipse_coords(2,:), 'r', 'LineWidth', 2);

    % % Configure the plot
    % axis equal;
    % xlabel('X');
    % ylabel('Y');
    % title('Covariance Ellipse');
    % legend('Data Points', 'Covariance Ellipse');
    % grid on;
    % hold off;
end
