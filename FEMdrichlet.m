% Clear the workspace
clc;
clear;

% Define symbolic variable
syms x;

% Define the function and its exact solution
f = @(x) 4.*x.^2 - 2.*x - 4;
exact_sol = @(x) -2.*x.^2 + x;

% Given limits and elements
a = 0;
b = 1;
n = 3;

% Dirichlet boundary conditions
ua = 0;
ub = -1;

% Calculate step size (h)
h = (b - a) / n;

% Generate node points
X = a:h:b;

% Initialize matrices for the finite element method
K = zeros(n+1, n+1); % Stiffness matrix
C = zeros(n+1, n+1); % Coefficient matrix
F = zeros(n+1, 1);   % Load vector
G = zeros(n+1, 1);   % Boundary condition vector

% Define coefficients for the matrices
K11 = (1/h) - 2*(-2*h/6);
K12 = (-1/h) - (-2*h/6);
C11 = 2*h/6;
C12 = h/6;

% Populate the matrices
K(1, 1) = 1;
K(n+1, n+1) = 1;

% Populate the interior nodes of the matrices
for i = 2:n
    K(i, i-1) = K12;
    K(i, i) = 2*K11;
    K(i, i+1) = K12;

    C(i, i-1) = C12;
    C(i, i) = 2*C11;
    C(i, i+1) = C12;
end

% Populate the load vector F
for i = 1:n+1
    F(i) = f(X(i));
end

% Apply Dirichlet boundary conditions
b_N = -C * F + G;
b_N(1) = ua;
b_N(n+1) = ub;

% Perform LU decomposition without using inbuilt functions
U = Gauss_Elimination(K,b_N);

% Display the computed values
disp('Computed Values:');
disp(U);

% Plot the calculated solution
figure;
plot(X, U);
hold on;

% Plot the exact solution
fplot(exact_sol, [a b]);

% Labeling and titling the plot
xlabel('x');
ylabel('u(x)');
legend('Calculated Solution', 'Exact Solution', 'Location', 'Best');
title("Solution of Dirichlet B.C. using FEM");
grid on;
hold off;

% Define the Gauss elimination function
function [X] = Gauss_Elimination(A,b)

[n, ~] = size(A);
X = zeros(n,1);

for i = 1:n
    for j = i+1:n
        factor = A(j,i)/A(i,i);
        A(j,:) = A(j,:) - factor*A(i,:);
        b(j) = b(j) - factor*b(i);
    end
end

X(n) = b(n)/A(n,n);

for i = n-1:-1:1
    sum = b(i);
    for j = i+1:n
        sum = sum - A(i,j)*X(j);
    end
    X(i) = sum/A(i,i);
end

end