clc;
clear;
% Define the grid
h = 1/3;
x = 0:h:1;
n = length(x);

% Initialize the matrix A and vector b
A = zeros(n,n);


b = zeros(n,1);


% Set up the system Au = b
for i = 2:n-1
    A(i, i-1:i+1) = [1 -2*(1+(h.^2)) 1];
    b(i) = (4*x(i)^2 - 2*x(i) - 4)*(h.^2);
end

% Apply Neumann boundary conditions
A(1,1) = 1;
b(1) = 0;
A(n, n-1) = 2;
A(n, n) = -2*(1+(h.^2));
b(n) = -2*h^2 + 2;


disp( 'A matrix = ' );
disp(A);
disp ( 'b matrix = ');
disp(b);

% Solve the system using Gauss-Jordan elimination
for i = 1:n
    pivot = A(i,i);
    A(i,:) = A(i,:) / pivot;
    b(i) = b(i) / pivot;
    for j = [1:i-1, i+1:n]
        ratio = A(j,i);
        A(j,:) = A(j,:) - ratio * A(i,:);
        b(j) = b(j) - ratio * b(i);
    end
end


% Print the computed values
disp('Computed values = ');
disp(b);

% Compute the exact solution
exact = -2*x.^2 + x;

% Plot the computed and exact solutions
figure;
plot(x, b, 'o-', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Computed');
hold on;
plot(x, exact, 'o--', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Exact');
hold off;

% Add grid lines
grid on;

% Add labels and title
xlabel('x', 'FontSize', 14);
ylabel('u', 'FontSize', 14);
title('Comparison of Computed and Exact Solutions(Neumann) ', 'FontSize', 16);

% Add a legend
legend('Location', 'best', 'FontSize', 12);

% Change the axes line width
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Change the figure color
set(gcf, 'Color', [0.9, 0.9, 0.9]);
