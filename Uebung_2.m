% Aufgabe 4
syms x1 x2
f(x1,x2) = 100*(x2 - x1.^2).^2 + (1 - x1).^2;
g(x1,x2) = gradient(f(x1,x2));
h(x1,x2) = hessian(f(x1,x2));

% a)
x_at=[-1,2];
f1 = f(x_at(1), x_at(2));
g1 = g(x_at(1), x_at(2));
h1 = h(x_at(1), x_at(2));
lin1(x1,x2) = expand(f1 + g1'*[x1 - x_at(1); x2 - x_at(2)]);

% b)
quad1(x1,x2) = expand(lin1(x1,x2) + 0.5*[x1 - x_at(1), x2 - x_at(2)]*h1*[x1 - x_at(1); x2 - x_at(2)]);

% c)
x0=[1,1];
f0 = f(x0(1), x0(2));
g0 = g(x0(1), x0(2));
h0 = h(x0(1), x0(2));
lin(x1,x2) = expand(f0 + g0'*[x1 - x0(1); x2 - x0(2)]);

% d)
quad(x1,x2) = expand(lin(x1,x2) + 0.5*[x1 - x0(1), x2 - x0(2)]*h0*[x1 - x0(1); x2 - x0(2)]);

fprintf("Rosenbrock-Funktion\n\tApproximation in [-1,2]:\n\t\t1D: %s\n\t\t2D: %s\n", lin1(x1,x2), quad1(x1,x2));
fprintf("\tApproximation in [1,1]:\n\t\t1D: %s\n\t\t2D: %s\n", lin(x1,x2), quad(x1,x2));

% Aufgabe 5
f = @(x) exp(-x) + 0.5*x.^2;
f_d = @(x) -exp(-x) + x;

syms f(x)
f(x) = exp(-x) + 0.5 * x.^2;
fdx = diff(f);
fddx = diff(fdx);

fprintf("Nullstelle mit Newtonverfahren: %0.4f\n", Newton1D(0,fdx));
fprintf("Nullstelle mit Goldener Schnitt: %0.4f\n", Golden(f,0,1,10));
fprintf("Nullstelle mit Matlab `fzero`: %0.4f\n", fzero(f_d,0));

function out_val = Newton1D(x_start, frst_deriv, scnd_deriv)
    steps = 10;
    h = 0.01;
    if ~exist('scnd_deriv', 'var')
        % second derivative exists
        scnd_deriv = @(x) ( double(frst_deriv(x + h)) - double(frst_deriv(x)) ) / h;
    end
    
    x = x_start;
    while steps > 0
        x = x - ( double(frst_deriv(x)) / double(scnd_deriv(x)) );
        steps = steps - 1;
    end
    
    out_val = x;
end

function x = Golden(f,a,b,iter)

    % Beispielaufruf: x=Golden(@(x)(x−1)^2,−1,3,10)
    g = (sqrt(5)-1)/2 ; % Verhaeltnis des Goldenen Schnitts 0.618
    lambda = a+(1-g)*(b-a);
    mu = a+g*(b-a);
    f_lambda = f(lambda);
    f_mu = f(mu); % Funktionswerte fuer Fallunterscheidung

    for i=1:iter
        fprintf("\t[%0.2f,%0.2f], f(lambda): %0.2f, f(mu): %0.2f\n", a, b, f_lambda, f_mu);
        if f_lambda > f_mu
            a=lambda;
            lambda=mu;
            mu=a+g*(b-a);
            f_lambda=f_mu;
            f_mu=f(mu);
        else
            b=mu;
            mu=lambda;
            lambda=a+(1-g)*(b-a);
            f_mu=f_lambda;
            f_lambda=f(lambda);
        end
    end
    x = (a+b)/2; %Ausgabewert
end