function A_inv = LU_inversion(A)
    %input: A (target matrix)
    %output: A_inv (inverse of target matrix)

    D = det(A); %determinant of A
    [m, n] = size(A); %check if A is a square matrix
    
    if (D ~= 0) && (m == n) %proceeds if A is invertible and square
        disp('Matrix A is invertible and amenable to the LU decomposition.')
        
        fprintf('\n');
        
        %check conditions for pivoting 
        condition_1a = (A == A');
        rand_x = ones(n,1);
        ran_index = randi(n);
        rand_x(ran_index) = 2;
        condition_1b = (rand_x' * A * rand_x >0);
        
        condition_2 = 0;
        for i = 1:n
            sum_row = sum(A(i,:)) - A(i,i);
            if (abs(A(i,i)) >= sum_row)
                condition_2 = condition_2 +1;
            end
        end
        
        if ((all(condition_1a(:)) && condition_1b) || (condition_2 == n)) %proceeds if A doesn't need pivoting
            disp('Matrix A does not need pivoting.')
            
            fprintf('\n');
        
            %calculation of L and U matrices
            L = zeros(n); %initialize dimensions of L and U
            U = eye(n); %the diagonal of U is ones, so we use eye() instead of zero()
            
            %populate L and U matrices
            for i = 1:n
                L(i,1) = A(i,1);
                U(1,i) = A(1,i)./L(1,1);
            end
        
            for i = 2:n
                for j = 2:i
                    sum_1 = 0;
                    for k = 1:j-1
                        sum_1 = sum_1 + L(i,k) * U(k,j);
                    end
                    L(i,j) = A(i,j) - sum_1;
                end
                for j = i:n
                    sum_2 = 0;
                    for k = 1:j-1
                        sum_2 = sum_2 + (L(i,k) * U(k,j));
                    end
                    U(i, j) = (A(i, j) - sum_2)/L(i,i);
                end
            end
            
            %display L and U matrices
            disp('L matrix is:');
            disp(L);
            disp('U matrix is:');
            disp(U);
        
            A_LU = L*U; %confirmation matrix to check LU calculations
            confirmation_1 = isequal(A, A_LU); %if true returns 1, if faulse 0
            if confirmation_1
                disp('LU calculations are correct.')
            else
                disp('LU calculations are incorrect.')
            end
         
            fprintf('\n');
            
            A_inv = zeros(n,n); %initialize inverse matrix A_inv
            
            %solving each column of the inverse matrix
            for j = 1:n
                b = zeros(n,1); %initialize b
                b(j) = 1; %each iteration value 1 is in the next row
            
                y = zeros(n,1); %initialize y
                
                %populate y vector
                y(1) = b(1) / L(1,1);
                for i = 2:n
                    sum_3 = 0;
                    for k = 1:i-1
                        sum_3 = sum_3 + L(i,k) * y(k);
                    end
                    y(i) = (b(i) - sum_3) / L(i,i);
                end
            
                x = zeros(n,1); %initialize x 
                
                %populate x vector
                x(n) = y(n) / U(n,n);
                for i = n-1:-1:1
                    sum_4 = 0;
                    for k = i+1:n
                        sum_4 = sum_4 + U(i,k) * x(k);
                    end
                
                    x(i) = (y(i) - sum_4) / U(i,i);
                end
                A_inv(:,j) = x;
            end
        
            I = eye(n); %identity matrix
            AA = A*A_inv; %confirmation matrix to check A_inv calculation 
            confirmation_2 = isequal(round(AA, 13), I); %checks if AA == I
            if confirmation_2
                disp('A^(-1) calculations are correct.');
            else
                disp('A^(-1) calculations are incorrect.');
            end
        
            fprintf('\n');
            
        else
            disp('Matrix A needs pivoting.');
        end
     
    elseif D == 0
        disp('Matrix A is not invertible.')
    elseif m ~= n
        disp('Matrix A cannot be inverted using the LU decomposition.')
    end
end