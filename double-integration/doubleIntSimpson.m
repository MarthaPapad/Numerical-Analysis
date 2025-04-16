function I = doubleIntSimpson(f, x_l, x_u, y_l, y_u)

    syms x y; % define symbolically
    
    % set number of segments
    n = 100;
    m = 100;
    
    h_x = (x_u - x_l)/n; %set step of x
    I_fl = 0; % sum of first and last term for outer integral
    I_even = 0; % sum of even terms for outer integral
    I_odd = 0; % sum of odd terms for outer integral
    
    % calculate outer integral
    for i = 0:n
        x_val = x_l + i*h_x; % value of x with each iteration
        h_y = (y_u - y_l)/m; % set step of y
        K_fl = subs(subs(f, y, y_l), x, x_val) + subs(subs(f, y, y_u), x, x_val); % sum of first and last term for inner integral
        K_even = 0; % sum of even terms for inner integral
        K_odd = 0; % sum of odd terms for inner integral
        
        % calculate inner integral
        for j = 1:(m-1)
            y_val = y_l + j*h_y; % value of y with each iteration
            Q = subs(subs(f, y, y_val), x, x_val); % f(x_val, y_val(x_val))
            
            if mod(j, 2) == 0 % check if term is even or odd
                K_even = K_even + Q;
            else
                K_odd = K_odd + Q;
            end
        end
        
        L = subs(h_y, x, x_val)/3 * (K_fl + 2*K_even + 4*K_odd); % value of inner integral for a value of x
        
        if i == 0 || i == n % check if term is first or last, even or odd
            I_fl = I_fl + L;
        elseif mod(i,2) == 0
            I_even = I_even + L;
        else
            I_odd = I_odd + L;
        end
    end
    
    I = double(h_x/3 * (I_fl + 2*I_even + 4*I_odd)); % value of the double integral
    
end