% Fach: Numerische Optimierung
% Projekt 4
%
% Autor: Maximilian Gaul
% Date: 02.07.2020
%---------------------------------------------------------------

function ret = GlobNewton(f, grad, hessian, x0, tol)
    
    rho = 1e-8;
    p = 2.1;
    
    x = x0;
    
    k = 1;
    kmax = 1000;
    
    ret = struct("x", x, "f", f(x), "gradient", grad(x), "hessian", hessian(x));

    while norm(grad(x)) > tol && k < kmax
       
        d = hessian(x) \ -grad(x);
        
        if ( grad(x)' * d ) > ( -rho * norm(d).^p )
            d = -grad(x);
        end
        
        phi = @(a) f(x + a * d);
        phi_grad = @(a) grad(x + a * d)' * d;
        alpha = Armijo(phi, phi_grad);
        
        x = x + alpha * d;
        
        k = k + 1;
        
        ret = [ ret; struct("x", x, "f", f(x), "gradient", grad(x), "hessian", hessian(x)) ];
        
    end

end