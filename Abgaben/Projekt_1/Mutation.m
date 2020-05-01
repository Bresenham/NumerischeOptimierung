% Function Name: Mutation
%
% Description: Performs Mutation-Selection on a function to find its
% minimum
%
% Inputs:
%     Function, Start value
%
% Version: MATLAB R2020a
% Author: Maximilian Gaul
% Date: 01.05.2020
%---------------------------------------------------------
function min_val = Mutation(f, start_x)
    
    k = 1;
    dims = numel(start_x);
    k_max = 1000000;
    alpha = 0.5;
    
    rand_vec = @(dim, from, to) from + to * rand(1, dim);
    hat = @(x) x + alpha * rand_vec(dims, -0.5, 0.5);
    
    x = start_x;
    
    while k < k_max
    
        x_hat = hat(x);
        if f(x_hat) < f(x)
            x = x_hat;
        end

        fprintf("\tx = [ %s] mit f(x) = %0.4f\n", sprintf("%0.4f ", x), f(x) );
        k = k + 1;
    end
    
    min_val = x;
end
