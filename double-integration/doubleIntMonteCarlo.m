function I = doubleIntMonteCarlo(f, x_l, x_u, y_l, y_u)
    
    syms x y; % define symbolically
    
    samples = 1000; % number of random samples
    count = 0; % sum of samples that are above f(x,y)
    
    for i = 1:samples
        
        % calculate f for random x, y
        if isnumeric(x_u)
            x_rand = x_l + (x_u - x_l) * rand();
            y_rand = subs(y_l, x, x_rand) + (subs(y_u, x, x_rand) - subs(y_l, x, x_rand)) * rand();
        else
            y_rand = y_l + (y_u - y_l) * rand();
            x_rand = subs(x_l, y, y_rand) + (subs(x_u, y, y_rand) - subs(x_l, y, y_rand)) * rand();
        end
        
        % calculate function f value for random x and y
        f_val = double(subs(f, [x, y], [x_rand, y_rand]));
        
        % calculate maximum and minimum value function f can take
        if isnumeric(x_u)
            f_max = double(subs(subs(f, y, y_u), x, x_u));
            f_min = double(subs(subs(f, y, y_l), x, x_l));
        else
            f_max = double(subs(subs(f, x, x_u), y, y_u));
            f_min = double(subs(subs(f, x, x_l), y, y_l));
        end
        
        % calculate random value for comparison
        z_rand = f_min + (f_max - f_min) * rand();
        
        % check if random value is lower than the function value 
        if z_rand <= f_val
            count = count + 1;
        end

    end 

    % calculate the area
    if isnumeric(x_u)
        area = 1/2 * (x_u - x_l) * (subs(y_u, x, x_u) - subs(y_l, x, x_l));
    else
        area = 1/2 * (y_u - y_l) * (subs(x_u, y, y_u) - subs(x_l, y, y_l));
    end
    
    %calculate integral
    I = double(count * area * (f_max-f_min)/samples);

end