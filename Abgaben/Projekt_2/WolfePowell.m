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

% phi = f(x + td)
% phi_grad = ∇f(x + td)' * d
function ret = WolfePowell(phi, phi_grad)

    rho = 0.3;
    gamma = 1.5;
    alpha = 1.0;
    sigma = 1e-2;
    
    % Armijo-Bedingung
    A = @(a) phi(a) <= phi(0) + sigma * a * phi_grad(0);
    
    % Wolfe-Powell-Bedingung
    W = @(a) phi_grad(a) >= rho * phi_grad(0);
    
    while(1)
        if( A(alpha) && W(alpha) )
            ret = alpha;
            return
        end
        if( A(alpha) && ~W(alpha) )
            alpha = gamma * alpha;
        elseif( ~A(alpha) )
            break
        end
    end
    
    a = 0;
    b = alpha;
    
    while(1)
        alpha = (a + b) / 2.0;
        if( A(alpha) && W(alpha) )
            ret = alpha;
            return
        end
        if( A(alpha) && ~W(alpha) )
            a = alpha;
        elseif( ~A(alpha) )
            b = alpha;
        end
    end
    
end