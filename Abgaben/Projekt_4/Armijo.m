% Function Name: Armijo
%
% Description: Armijo method for stepsize calculation
%
% Inputs:
%   Function phi and phi_grad (derivative of phi)
% Outputs:
%   step size
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     02.07.2020
%---------------------------------------------------------

% phi = f(x + ad)
% phi_grad = ∇f(x + ad)' * d
function ret = Armijo(phi, phi_grad)

    beta = 0.5;
    sigma = 1e-4;
    
    alpha = 1.0;
    
    phi_zero = phi(0);
    phi_grad_zero = phi_grad(0);

    % Armijo-Bedingung
    A = @(a) phi(a) <= phi_zero + sigma * a * phi_grad_zero;
    
    while(1)
        
        if A(alpha)
            ret = alpha;
            return;
        end
        
        alpha = beta * alpha;
        
    end
    
end