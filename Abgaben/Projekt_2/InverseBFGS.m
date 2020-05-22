% Function Name: InverseBFGS
%
% Description: Inverse BFGS Quasi-Newton method
%
% Inputs:
%   Function, Gradient, Startvalue x0
% Outputs:
%   for every step: x_n, f(x_n), grad(x_n)
%
% WolfePowell is defined in 'WolfePowell.m'
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
    k_max = 1000000;

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
    
    % Zur Sicherheit hier auch die Anzahl an Iterationen beschränken
    while norm( grad(x) ) > 1e-8 && k < k_max
                
        % Update macht im 1. Durchlauf wenig Sinn
        if k >= 1
            s_val = s(x, x_old);
            y_val = y(x, x_old);
            B = update(B, s_val, y_val);
        end
        
        % Abstiegsrichtung hier nur ein Matrix-Vektor-Produkt da B der
        % inversen von A (Hesse-Matrix Approximation) entspricht
        d = -( B * grad(x) );
        
        % Prüfe ob es sich bei 'd' um eine Abstiegsrichtung handelt, wenn
        % nicht, verwende den negativen Gradienten und setze die
        % Approximation der Hesse-Matrix auf die Einheitsmatrix zurück
        if ( grad(x)' * d > - norm(d) )
            d = - grad(x);
            B = eye(dim);
        end
        
        % Definition der Funktionen phi und phi' die für Wolfe-Powell
        % benötigt werden
        phi = @(a) f(x + a * d);
        phi_grad = @(a) grad(x + a * d)' * d;
        alpha = WolfePowell(phi, phi_grad);
        
        x_old = x;
        x = x + alpha * d;
        
        k = k + 1;
    end
    
    ret = x;
end