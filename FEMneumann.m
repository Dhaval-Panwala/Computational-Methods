% Clearing the command window and workspace
clc
clear

% Symbolic variable declaration
syms x;

% Defining the function f(x) and its exact solution
f = @(x) 4*x.^2 - 2*x - 4;
exact_sol = @(x) -2*x.^2 + x;

% Given limits and elements for the finite element method
a = 0;
b = 1;
n = 3;  % Number of elements

% Neumann boundary condition values
ua = 0;       % u at x = 0
dudx = -3;    % du/dx at x = 1

% Calculating the element size (h)
delx = (b - a) / n;
X = a:delx:b;  % Generating nodes

% Initializing matrices K, C, F, and G
K = zeros(n+1, n+1);
C = zeros(n+1, n+1);
F = zeros(n+1, 1);
G = zeros(n+1, 1);

K = zeros(n+1, n+1);
K11 = (1/delx)-2*(-2*delx/6);
K12 = (-1/delx)-(-2*delx/6);
K(1, 1) = 1;
K(n+1, n+1) = K11; 
K(n+1,n) = K12;

for i = 2:n
    K(i, i-1) = K12;
    K(i, i) = 2*K11;
    K(i, i+1) = K12;
end  

% Initialize the C matrix for the system g=-CF
C = zeros(n+1, n+1);
C11 = 2*delx/6;
C12 = delx/6;
C(1, 1) = C11;
C(1, 2) = C12;
C(n+1, n+1) = C11;
C(n+1, n) = C12;

for i = 2:n
    C(i, i-1) = C12;
    C(i, i) = 2*C11;
    C(i, i+1) = C12;
end


% Calculating the load vector F
for i = 1:n+1
    F(i) = f(X(i));
end

% Neumann boundary condition
G(n+1) = dudx;

% Combining the matrices to form the system of equations Ku = b_N
b_N = -C * F + G;
b_N(1) = ua;

% Solving the system using Gaussian Elimination
U = myGaussElimination(K, b_N);

% Displaying the computed values
disp('U Values:');
disp(U);

% Plotting the calculated solution along with the exact solution
figure;
plot(X, U, 'o-', 'LineWidth', 2);
hold on;
fplot(exact_sol, [a b], 'r--', 'LineWidth', 2);
xlabel('x');
ylabel('u(x)');
legend('Calculated Solution', 'Exact Solution', 'Location', 'Best');
title('Solution of Neumann B.C. using FEM');
grid on;
hold off;

% Custom Gaussian Elimination function
function X = myGaussElimination(A, b)
    n = length(b);
    X = zeros(n, 1);

    % Forward elimination
    for k = 1:n-1
        for i = k+1:n
            factor = A(i, k) / A(k, k);
            A(i, k:n) = A(i, k:n) - factor * A(k, k:n);
            b(i) = b(i) - factor * b(k);
        end
    end

    % Backward substitution
    X(n) = b(n) / A(n, n);
    for i = n-1:-1:1
        X(i) = (b(i) - A(i, i+1:n) * X(i+1:n)) / A(i,i);
    end
end