% Ypologistika Mathhmatika I
% Project 1
% Martha Papadopoulou

clear
clc

figure(1);
fplot(@(x) exp(x)+x-2);
xlabel('x');
ylabel('f(x)');
title('Plot of f(x) = e^{x} + x - 2');
grid on;
hold on;

%From the plot above we can see that the root is between 0<x<1.
%So for simplicity we will use those as are starting points.

n = 6; %significant figures
es = 0.5*10^(2-n); %prespecified percent tolerance
f = @(x) exp(x)+x-2;

%Bisection Method
xl_b = 0; %lower guess
xu_b = 1; %upper guess
ea_b = 100; %starting value to enter the loop
ea_values_b = []; %approximation error table
count_b = 0; 

while ea_b>es && count_b <= 1000
    
    xr_b = (xu_b+xl_b)/2;
    
    if f(xl_b)*f(xr_b)<0
        xu_b = xr_b;
    elseif f(xl_b)*f(xr_b)>0
        xl_b = xr_b;
    end
    
    ea_b = abs((xu_b-xl_b)/(xu_b+xl_b))*100;
    ea_values_b = [ea_values_b, ea_b];
    count_b = count_b + 1;
    
end

root_b = chop(xr_b, n);
    
if count_b == 1001
    disp('Not able to find a root within 1000 iterations.')
else
    figure(2);
    plot(1:count_b,ea_values_b, '-o');
    xlabel('Number of repetitions')
    ylabel('Approximate error \epsilon_{a}%')
    title('Convergence of Bisection Method')
    grid on;

    fprintf('The root using the Bisection Method for %d significant figures is: %.6f\n', n, root_b);
    disp(['Number of repetitions using the Bisection Method: ' num2str(count_b)]);
end

%False Position Method (Regula Falsi)
xl_r = 0; %lower guess
xu_r = 1; %upper guess
ea_r = 100; %starting value to enter the loop
ea_values_r = []; %approximation error table
count_r = 0;
xur_count = 0; %counts how many times xu_r changes in a row
xlr_count = 0; %counts how many times xl_r changes in a row
fl_r = f(xl_r);
fu_r = f(xu_r);

while ea_r>es && count_r <= 1000
    
    xr_r = xu_r - (fu_r*(xl_r - xu_r))/(fl_r - fu_r);
    
    if fl_r*f(xr_r)<0
        xu_r = xr_r;
        fu_r = f(xr_r);
        xur_count = xur_count + 1;
            
    elseif fl_r*f(xr_r)>0
        xl_r = xr_r;
        fl_r = f(xr_r);
        xlr_count = xlr_count + 1;
    end

    ea_r = abs((xu_r - xl_r)/(xu_r + xl_r))*100;
    ea_values_r = [ea_values_r, ea_r];
    count_r = count_r + 1;
    
    if xur_count == 2
        fl_r = fl_r/2;
        xur_count = 0;
    end
    
    if xlr_count == 2 
        fu_r = fu_r/2;
        xlr_count = 0;
    end
    
end

root_r = chop(xr_r, n);

if count_r == 1001
    disp('Not able to find a root within 1000 iterations.')
else    
    figure(3);
    plot(1:count_r, ea_values_r, '-o');
    xlabel('Number of repetitions')
    ylabel('Approximate error \epsilon_{a}%')
    title('Convergence of False Position Method')
    grid on;

    fprintf('The root using the False Position Method for %d significant figures is: %.6f\n', n, root_r);
    disp(['Number of repetitions using the False Position Method: ' num2str(count_r)]);
end

%Fixed-Point Iteration
g = @(x) log(2-x); %so that near the starting point |g'(x)|<1
x_f = 0; %starting point
ea_f = 100; %starting value to enter the loop
ea_values_f = []; %approximation error table
count_f = 0;

while ea_f>es && count_f <= 1000
    
    ea_f = abs((g(x_f)-x_f)/g(x_f))*100;
    ea_values_f = [ea_values_f, ea_f];
    x_f = g(x_f);
    count_f = count_f + 1;
    
end

root_f = chop(x_f, n);

if count_f == 1001
    disp('Not able to find a root within 1000 iterations.')
else
    figure(4);
    plot(1:count_f, ea_values_f, '-o');
    xlabel('Number of repetitions')
    ylabel('Approximate error \epsilon_{a}%')
    title('Convergence of Fixed-Point Iteration')
    grid on;
    
    fprintf('The root using the Fixed-Point Iteration for %d significant figures is: %.6f\n', n, root_f);
    disp(['Number of repetitions using the Fixed-Point Iteration: ' num2str(count_f)]);
end

hold off;

result = (root_b==root_r) && (root_r==root_f);

if result == 1
    disp('All methods converge to a consensus value of the root.')
    disp(['The value of the f(x) = exp(x)+x-2 for x = ' num2str(root_b) ' is f(x) = ' num2str(f(root_b))])
end
