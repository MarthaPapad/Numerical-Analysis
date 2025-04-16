% Ypologistika Mathhmatika I
% Project 2
% Martha Papadopoulou
% AEM 4438

clear; clc;

A = [4 1 2 3 5; 1 3 1 4 2; 2 1 5 2 3; 3 4 2 4 1; 5 2 3 1 5]; %target matrix

fprintf('Task 1\n -------------------------------------------------------\n');

A_inv = LU_inversion(A); %returns the inverted A
disp('The inverse A matrix using the LU decomposition inversion is:');

%if the calculations were successful displays the inverted A
if A_inv
    disp(A_inv);
end

fprintf('Task 2\n -------------------------------------------------------\n');

[l, v] = power_method(A); %returns the dominant eigenvalue l and corresponding eigenvector v

%if the calculations were successful displays l and v
if l
    fprintf('The largest eigenvalue of matrix A using the power method is: %.8f\n', l);
    fprintf('\n');
    disp('The corresponding eigenvector is:');
    disp(['[', num2str(v'), ']^T']);
end
