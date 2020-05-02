% Function Name: Mutation
%
% Description: Performs Mutation-Selection on a function to find a
% minimum
%
% Inputs:
%     Function, Start value
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     01.05.2020
%---------------------------------------------------------
function [min_val, k] = Mutation(f, start_x)
    
    k = 1;
    k_max = 15000;
    alpha = 0.125;
    dim = numel(start_x);
    
    rand_vec = @(dim, from, to) from + (to - from) * rand(1, dim);
    hat = @(x) x + alpha * rand_vec(dim, -0.5, 0.5);
    
    x = start_x;
    f_x_val = 1000;
    f_x_old_val = -1000;

    % Abbruchkriterien: Differenzbetrag oder Anzahl Iterationen
    while norm(f_x_val - f_x_old_val) > 1e-6
    % while k <= k_max
    
        x_hat = hat(x);
        
        f_x = f(x);
        f_x_hat = f(x_hat);
        
        if f_x_hat < f_x
            f_x_old_val = f_x;
            f_x_val = f_x_hat;
            x = x_hat;
            % fprintf("\tx = [ %s] mit f(x) = %0.8f\n", sprintf("%0.4f ", x), f(x));
        end
        
        % Verstärkungsfaktor wird bis zu einer gewissen Grenze auf alpha
        % multipliziert
        alpha = min([alpha * 1.05, 0.5]);
        
        k = k + 1;
    end
    
    min_val = x;
end
