% Fach: Numerische Optimierung
% Projekt 3 - Penalty-Verfahren & SQP-Verfahren
%
% Autor: Maximilian Gaul
% Date: 12.06.2020
%---------------------------------------------------------------

% Aufgabe 1
% Siehe PDF Dokument

% Aufgabe 2
% Siehe BFGS_Pen.m und ArmijoPen.m

% Aufgabe 3

fprintf("--------------------AUFGABE 3--------------------\n");

r = 5000;
max_iters = 1e5;

f = @(x) 1000 - x(1).^2 - 2 * x(2).^2 - x(3).^2 - x(1) * x(2) - x(1) * x(3);
grad_f = @(x) [ -2 * x(1) - x(2) - x(3), -4 * x(2) - x(1), -2 * x(3) - x(1) ];

f1 = @(x) x(1).^2 + x(2).^2 + x(3).^2 - 25;
grad_h1 = @(x) [ 2 * x(1), 2 * x(2), 2 * x(3) ];

f2 = @(x) 8 * x(1) + 14 * x(2) + 7 * x(3) - 56;
grad_h2 = @(x) [8, 14, 7];

pen_f = @(x, r) f(x) + (0.5 * r) * ( f1(x).^2 + f2(x).^2 );
grad_pen_f = @(x, r) grad_f(x) + r * ( f1(x) * grad_h1(x) + f2(x) * grad_h2(x) );

m_pf = @(x, r) min_pen_f(pen_f, grad_pen_f, x, r);

x0 = [25; 35; 45];
ret = BFGS_Pen(m_pf, r, x0, max_iters);
res = ret.x(end,:);
fprintf("BFGS_PEN returned x=%s with f(x)=%0.4f and f1(x)=%0.4f, f2(x)=%0.4f\n", vec2str(res), f(res'), f1(res'), f2(res'));

% Lösung mit fmincon
confunF = @(x) confuneqF(f1, f2, x);
ret = fmincon(f,x0,[],[],[],[],[],[],confunF);
fprintf("fmincon returned x=%s with f(x)=%0.4f and f1(x)=%0.4f, f2(x)=%0.4f\n", vec2str(ret), f(ret'), f1(ret'), f2(ret'));

% Lösung mit fminunc bei festem r
options = optimoptions("fminunc", "OptimalityTolerance", 1e-6, "MaxFunctionEvaluations", 1e+8, "MaxIterations", max_iters);
pen_f_f = @(x) pen_f(x, r);
ret = fminunc(pen_f_f, x0, options);
fprintf("fminunc returned x=%s with f(x)=%0.4f and f1(x)=%0.4f, f2(x)=%0.4f\n", vec2str(ret), f(ret'), f1(ret'), f2(ret'));

% Aufgabe 4
% Siehe PDF

% Aufgabe 5
% Siehe PDF

% Aufgabe 6

fprintf("--------------------AUFGABE 6--------------------\n");

g = @(x) x(1);

g1 = @(x) x(1).^2 + x(2).^2 - 1;
g2 = @(x) x(1) + x(2) - 1;

x0 = [12; 12];
conNeqG = @(x) confunNeqG(g1, g2, x);
ret = fmincon(g, x0, [], [], [], [], [], [], conNeqG);
fprintf("fmincon returned x=%s with g(x)=%0.4f and g1(x)=%0.4f, g2(x)=%0.4f\n", vec2str(ret), g(ret'), g1(ret'), g2(ret'));

confunG = @(x) confuneqG(g1, g2, x);
ret = fmincon(g, x0, [], [], [], [], [], [], confunG);
fprintf("fmincon returned x=%s with g(x)=%0.4f and g1(x)=%0.4f, g2(x)=%0.4f\n", vec2str(ret), g(ret'), g1(ret'), g2(ret'));

% Aufgabe 7
% Siehe PDF

% Aufgabe 8
% Siehe PDF

t = 0.75;
m = @(x) ( x(1) - 1.5 ).^2 + ( x(2) - t ).^4;
m1 = @(x) -1 + x(1) + x(2);
m2 = @(x) -1 + x(1) - x(2);
m3 = @(x) -1 - x(1) + x(2);
m4 = @(x) -1 - x(1) - x(2);

conNeqM = @(x) confunNeqM(m1, m2, m3, m4, x);
ret = fmincon(m, [2; 3], [], [], [], [], [], [], conNeqM);
fprintf("fmincon returned x=%s with m(x)=%0.4f and m1(x)=%0.4f, m2(x)=%0.4f, m3(x)=%0.4f, m4(x)=%0.4f\n", vec2str(ret), m(ret), m1(ret), m2(ret), m3(ret), m4(ret));

% Aufgabe 9
% Siehe PDF

m_deriv = @(x) [2 * ( x(1) - 1.5 ); 4 * ( x(2) - 0.75).^3];
m_hessian = @(x) [ 2, 0; 0, 12 * ( x(2) - 0.75).^2 ];

x0 = [0; 0];
n = x0;

g = [m1(n), m2(n), m3(n), m4(n)];
grad_g = [1, 1; 1, -1; -1, 1; -1, -1];

% m_f = @(x) m_deriv(n)' * x + 0.5 * (x') * m_hessian(n) * x;
% ret = fmincon(m_f, zeros(2,1), grad_g, -g);
ret = quadprog(m_hessian(n), m_deriv(n), grad_g, -g);

n = n + ret;
fprintf("quadprog returned x=%s with g(x)=%0.4f and g1(x) = %0.4f, g2(x) = %0.4f, g3(x) = %0.4f, g4(x) = %0.4f\n", vec2str(n), m(n), m1(n), m2(n), m3(n), m4(n));

% Aufgabe 10

t = 0.75;
m = @(x) ( x(1) - 1.5 ).^2 + ( x(2) - t ).^4;

options = optimoptions("fmincon", "Algorithm", "sqp");
ret = fmincon(m, [2; 3], [], [], [], [], [], [], conNeqM, options);
fprintf("fmincon returned x=%s with m(x)=%0.4f and m1(x)=%0.4f, m2(x)=%0.4f, m3=%0.4f, m4=%0.4f\n", vec2str(ret), m(ret'), m1(ret'), m2(ret'), m3(ret'), m4(ret'));

% Ungleichheitsbedingungen aus Aufgabe 8 für fmincon
function [c,ceq] = confunNeqM(m1, m2, m3, m4, x)
    % Nonlinear inequality constraints
    c = [m1(x), m2(x), m3(x), m4(x)];
    % Nonlinear equality constraints
    ceq = [];
end

% Ungleichheitsbedingungen aus Aufgabe 4 für fmincon
function [c,ceq] = confunNeqG(g1, g2, x)
    % Nonlinear inequality constraints
    c = [g1(x), g2(x)];
    % Nonlinear equality constraints
    ceq = [];
end

% Gleichheitsbedingungen aus Aufgabe 5 für fmincon
function [c,ceq] = confuneqG(g1, g2, x)
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [g1(x), g2(x)];
end

% Gleichheitsbedingungen aus Aufgabe 2 für fmincon
function [c,ceq] = confuneqF(f1, f2, x)
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [f1(x), f2(x)];
end

% Minimierungsfunktion für BFGS_Pen
% \input: Funktion f, Gradient von f, x, r
function [f, g] = min_pen_f(f_pen, g_pen, x, r)

    f = f_pen(x, r);
    g = g_pen(x, r);

end

% Formatiert einen Vector für die Ausgabe auf der Kommandozeile
% \input: 1D Vector
function str = vec2str(vec)
    
    dim = numel(vec);
    ret_str = "[";
    for i = 1:dim
        if i == dim
            ret_str = append( ret_str, sprintf("%0.4f]", vec(i)) );
        else
            ret_str = append( ret_str, sprintf("%0.4f, ", vec(i)) );
        end
    end
    
    str = ret_str;
end