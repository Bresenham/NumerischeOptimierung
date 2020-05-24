% Function Name: WolfePowell
%
% Description: Wolfe Powell method for stepsize calculation
%
% Inputs:
%   Function phi and phi' (its derivative)
% Outputs:
%   step size
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     22.05.2020
%---------------------------------------------------------

% phi = f(x + ad)
% phi_grad = ∇f(x + ad)' * d
function ret = WolfePowell(phi, phi_grad)

    rho = 0.9;
    gamma = 1.5;
    alpha = 1.0;
    sigma = 1e-4;
    phi_zero = phi(0);
    phi_grad_zero = phi_grad(0);

    % Armijo-Bedingung
    A = @(a) phi(a) <= phi_zero + sigma * a * phi_grad_zero;
    
    % Wolfe-Powell-Bedingung
    W = @(a) phi_grad(a) >= rho * phi_grad_zero;
    
    while(1)
        if( A(alpha) && W(alpha) )
            ret = alpha;
            return
        elseif( A(alpha) && ~W(alpha) )
            alpha = gamma * alpha;
        elseif( ~A(alpha) )
            break;
        end
    end
    
    a = 0;
    b = alpha;
    
    while(1)
        alpha = (a + b) / 2.0;
        % Sicherheitscheck damit man nicht in eine Endlosschleife endet
        if( abs(a - b) <= 1e-10 || alpha <= 1e-12 )
            ret = alpha;
            return;
        end
        if( A(alpha) && W(alpha) )
            ret = alpha;
            return
        elseif( A(alpha) && ~W(alpha) )
            a = alpha;
        elseif( ~A(alpha) )
            b = alpha;
        end
    end
    
end