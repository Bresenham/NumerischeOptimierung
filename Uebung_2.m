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

fprintf("Rosenbrock-Funktion Approximation in [-1,2]:\n\t1D: %s\n\t2D: %s\n", lin1(x1,x2), quad1(x1,x2));
fprintf("Rosenbrock-Funktion Approximation in [1,1]:\n\t1D: %s\n\t2D: %s\n", lin(x1,x2), quad(x1,x2));