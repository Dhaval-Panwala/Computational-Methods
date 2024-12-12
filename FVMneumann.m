clc
clear
% defining function Given in question
syms x;
f = @(x) (4.*x.^2) - (2.*x) - 4;
exactSol = @(x) (-2.*x.^2) + x;

% upper and lower limit of x and no. of elements
a = 0;
b = 1;
n = 3;

% Neumann boundary condition ( Given ) 
ua = 0;
dudx = -3;


% Defining step size
h = (b-a)/n;
X = a:h:b;

% defining A and making A matrix
A = zeros(n+1, n+1);
A(1,1) = 1;
A2 = -1*((2/(h^2))-(-2));
A1 = 1/(h^2);
A(n+1,n) = 2*A1;
A(n+1,n+1) = A2;

% defining b and making b matrix
b = zeros(n+1, 1);
b(1) = ua;
b(n+1) = f(X(n+1)) - 2*dudx./h;    


for i = 2:n
    A(i, i-1) = A1;
    A(i, i) = A2;
    A(i, i+1) = A1;  
    b(i)= f(X(i));     
end  

%Solving AU = b using Gauss Elimination
U = Guass_Elimination(A,b);

% displaying results
disp ("Computed values of U");
disp(U);   

% Plotting the comparison betwween computed and exact solution
figure;
plot(X, U); 
hold on;
fplot(exactSol, [0 1]);
xlabel("x");
ylabel("u(x)");
title("Solution of Neumann B.C. using FVM(via FDM)");
legend('Computed Solution', 'Exact solution');
grid on;

% Gauss Elimination function
function [X] = Guass_Elimination(A,b)

    [n, ~] = size(A); 

    % Forward Elimination
    for i = 1:n
        for j = i+1:n
            fac = A(j, i)/A(i, i); % Define factor
            A(j, i) = 0; % Make the A(i, j) term 0 for j<i
            for k = i+1:n
                A(j, k) = A(j, k) - fac*A(i, k); % Define Row Operations for A Matrix
            end
            b(j) = b(j) - fac*b(i); % Define Row Operations for b Matrix
        end
    end

    % Backward Substitution
    X(n) = b(n)/A(n,n); % Define X_n

    for i = n-1:-1:1
        sum = b(i); % Define variable to calculate X_n-1, X_n-2,....X_1
        for j = i+1:n
            sum = sum - A(i,j)*X(j); % Define calculation for X_n-1, X_n-2,....X_1
        end
        X(i) = sum/A(i,i); % Define X(i)    
    end

end
