% Aufgabe 1
% f(x) = 2x1^2 - 4x1x2 + 4x2^2
% Q1 = [4, -4; -4, 8];
% q1 = [0; 0];

% f(x) = x1^2 + 20x2^2 - 4x1 - 2x2
% Q2 = [0.5, 0; 0, 10];
% q2 = [-4; -2];
% x0 = [-5; 1];

% ret = ConjugateGradient(Q2, q2, x0);
% fprintf("RESULT: [%0.6f, %0.6f]\n", ret(1), ret(2));
% 
% disp( pcg(Q2, -q2) );

% Aufgabe 2
f = @(x) 2*x(1).^2 - 4*x(1)*x(2) + 4*x(2).^2;
grad = @(x) [4*x(1) - 4*x(2); -4*x(1) + 8*x(2)];
ret = ConjugateGradientWolfe(f, grad, [1;1]);

% Funktion Aufgabe 1
function ret = ConjugateGradient(Q, q, x0)
    
    k = 0;
    x = x0;
    grad_f = Q*x + q;
    d = -grad_f;

    while norm(grad_f) > 1e-8 && k < 100
        alpha = - ( grad_f' * d ) / ( d' * (Q*d) );
        x = x + alpha * d;
        grad_f_new = grad_f + alpha*Q*d;
        beta = ( norm(grad_f_new)^2 ) / ( norm(grad_f)^2 );
        d = -grad_f_new + beta * d;
        grad_f = grad_f_new;

        k = k + 1;
    end
    
    ret = x;
end

% phi = f(x + td)
% phi_grad = ∇f(x + td)^T * d
function ret = WolfePowell(phi, phi_grad)

    rho = 0.8;
    gamma = 1.5;
    alpha = 1.0;
    sigma = 1e-3;
    
    A = @(a) phi(a) <= phi(0) + sigma * a * phi_grad(0);
    W = @(a) phi_grad(a) >= rho * phi_grad(0);
    
    while(1)
        if(A(alpha) && W(alpha))
            ret = alpha;
            return
        end
        if(A(alpha) && ~W(alpha))
            alpha = gamma * alpha;
        elseif(~A(alpha))
            break
        end
    end
    
    a = 0;
    b = alpha;
    
    while(1)
        alpha = (a + b) / 2.0;
        if(A(alpha) && W(alpha))
            ret = alpha;
            return
        end
        if(A(alpha) && ~W(alpha))
            a = alpha;
        elseif(~A(alpha))
            b = alpha;
        end
    end
    
end

% Funktion Aufgabe 2
function ret = ConjugateGradientWolfe(f, grad, x0)
    
    k = 0;
    x = x0;
    d = -grad(x0);
    
    while norm( grad(x0) ) > 1e-8 && k < 100
        
        phi = @(a) f(x + a * d);
        phi_grad = @(a) grad(x + a * d)' * d;

        alpha = WolfePowell(phi, phi_grad);
        
        x_new = x + alpha * d;
        beta = norm( grad(x_new) )^2 / norm( grad(x) )^2;
        d = -grad(x_new) + beta * d;
        x = x_new;
        disp(x);
        k = k + 1;
    end
    
    ret = x;
end