% Ypologistika Mathhmatika I
% Project 3
% Martha Papadopoulou
% AEM 4438

clear; clc; close all;

syms x y %define variables symbolically
f = x^2 + y^3; % define integrand

%define integral bounds
x_lower = 0;
x_upper = 1;
y_lower = x;
y_upper = 2*x;

%% 3D plot of integrand within bounds
x_values = linspace(x_lower, x_upper, 100);
y_values = linspace(subs(y_lower, x, x_lower), subs(y_upper, x, x_upper), 100);
[X, Y] = meshgrid(double(x_values), double(y_values)); % create a grid of values
mask = (Y >= X) & (Y <= 2*X); % restriction for y

F = double(X.^2 + Y.^3);
F(~mask) = NaN;

figure(1);
surf(X, Y, F); % 3d plot of f(x,y)
colormap(spring); % colour of surface 
title('3D Plot of f(x, y) = x^2 + y^3 within integral bounds');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');

%% (a) solve analytically
g = int(f, y, y_lower, y_upper);
result_a = int(g, x, x_lower, x_upper);
result_a = double(result_a);
fprintf('The analytical result of the double integral is: %d\n', result_a);

%% (b) solve numerically with Simpson
result_b = doubleIntSimpson(f, x_lower, x_upper, y_lower, y_upper);
accuracy_b = abs(result_a - result_b);
fprintf('The numerical result of the double integral, using Simpson, is: %f\n', result_b);
fprintf('Accuracy of Simpsons rule is %e\n', accuracy_b);

%% (c) solve numerically with Monte Carlo
n = 100; % number of iterations
results_c = zeros(n,1); % initialize results vector

for i = 1:n
    results_c(i) = doubleIntMonteCarlo(f, x_lower, x_upper, y_lower, y_upper);
end
result_c = mean(results_c); % mean value of results 
variance_c = var(results_c); % variance of results
fprintf('The numerical result of the double integral, using Monte Carlo integration %d times, is %f\n', n, result_c);
fprintf('The variance of the method is %f\n', variance_c);

% scatter plot of results for each iteration
figure(2);
scatter((1:n),results_c);
xlabel('Iteration number n');
ylabel('Result of integral I');
title('Scatter plot of the results of Monte Carlo integration for every iteration');
grid on;
