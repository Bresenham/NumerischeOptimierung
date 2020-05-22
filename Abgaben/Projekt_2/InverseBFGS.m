% Function Name: InverseBFGS
%
% Description: Inverse BFGS Quasi-Newton method
%
% Inputs:
%   Function, Gradient, Startvalue x0
% Outputs:
%   for every step: x_n, f(x_n), grad(x_n)
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     22.05.2020
%---------------------------------------------------------

function ret = InverseBFGS(f, grad, x0)

    k = 0;
    x = x0;
    x_old = x0;
    dim = numel(x0);
    
    % Start mit der Einheitsmatrix als inverse zur Approximation der
    % Hesse-Matrix
    B = eye(dim);
    
    % Funktionen phi und phi' zur Berechnung der Schrittweite nach Wolfe
    % und Powell
    
    s = @(x, x_old) x - x_old;
    y = @(x, x_old) grad(x) - grad(x_old);
    
    update1 = @(B, s, y) ( (s - B*y) * s' + s * (s - B*y)' ) / ( y' * s);
    update2 = @(B, s, y) ( (s - B*y)' * y * s * (s') ) / ( y' * s).^2;

    update = @(B, s, y) B + update1(B, s, y) - update2(B, s, y);
    
    while norm( grad(x) ) > 1e-8
        
        % Update macht im 1. Durchlauf wenig Sinn
        if k >= 1
            s_val = s(x, x_old);
            y_val = y(x, x_old);
            B = update(B, s_val, y_val);
        end
        
        d = -( B * grad(x) );
        
        phi = @(a) f(x + a * d);
        phi_grad = @(a) grad(x + a * d)' * d;
        alpha = WolfePowell(phi, phi_grad);
        
        x_old = x;
        x = x + alpha * d;
        
        k = k + 1;
    end
    
    ret = x;
end