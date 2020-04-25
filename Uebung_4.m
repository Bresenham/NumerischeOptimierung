f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2)];

ret = gradient_armijo(f_rosen, f_rosen_grad);
fprintf("Result: [%0.4f,%0.4f]\n", ret(1), ret(2));

function ret_val = gradient_armijo(fnct, grad)
    beta = 0.5;
    sig = 1e-4;
    
    x = [-1.2;1];
    k = 0;
    while ( norm(grad(x)) > 1e-8 )
        d = -grad(x);
        
        phi = @(alpha) fnct(x + alpha*d);
        phi_grad = @(alpha) (grad(x + alpha*d).') * d;
        alpha = armijo(phi, phi_grad, beta, sig);
        
        x = x + alpha * d;
        
        fprintf("\tStep %d: [%0.4f,%0.4f] @ %0.4f\n", k, x(1), x(2), fnct(x));
        k = k + 1;
    end
    
    ret_val = x;
end

function ret_val = armijo(phi, phi_grad, beta, sig)
    alpha = 1;
    
    while ~( phi(alpha) <=  phi(0) + sig * alpha * phi_grad(0) )
        alpha = beta * alpha;
    end
    
    ret_val = alpha;
end