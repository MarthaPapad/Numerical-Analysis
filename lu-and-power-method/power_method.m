function [lambda, v] = power_method(A)
    %input: A (target matrix)
    %output: lambda(largest eigenvalue), v (corresponding eigenvector)

    D = det(A); %determinant of A
    [m, n] = size(A); %check if A is a square matrix
    
    if (D ~= 0) && (m == n) %proceeds if A is invertible and square
        disp('Matrix A is amenable to the power method.')
        fprintf('\n');
        
        k = 10; %significant figures
        es = 0.5*10^(2-k); %prespecified percent tolerance
        
        x = ones(n, 1); %initialize random vector
        x(1) = 2;
        x_norm = x/norm(x);
        
        ea = 100; %starting value to enter loop, not added to ea_values
        i = 0; %number of iterations 
        ea_values = []; %empty list to enter ea values in loop
        x_prev_i = A*x_norm; %first x^(i)
        
        %power method iteration
        %runs until ea<es or i>300
        while ea >= es && i <= 300
            
            x_i = A * x_prev_i;
            [max_xi, max_index] = max(x_i); %find max element in x_i and index
            
            lambda = max_xi/x_prev_i(max_index); %eigenvalue
            
            %calculation of approximate percent error
            if i > 0 
                ea = abs((lambda - lambda_prev)/lambda)*100;
                ea_values = [ea_values, ea];
            end
            
            %update values for next iteration
            lambda_prev = lambda;
            x_prev_i = x_i;
            i = i+1;
        end
        
        v = x_i./max_xi; %corresponding eigenvector
        
        if i > 300
            disp('Maximum number of iterations reached.');
        else
            fprintf('Iterations needed for accuracy of %d significant figures is: %d\n', k,i);
            fprintf('\n');
        end
        
        %plot of convergence of the power method
        figure(1);
        plot(2:(length(ea_values)+1),ea_values, '-o');
        xlabel('Number of repetitions')
        ylabel('Approximate error \epsilon_{a}%')
        title('Convergence of Power Method')
        grid on;
        
        
    elseif D == 0
        disp('Matrix A is not invertible.')
    elseif m ~= n
        disp('Matrix A cannot be inverted using the LU decomposition.')
    end

end