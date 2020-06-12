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

f = @(x) 1000 - x(1).^2 - 2 * x(2).^2 - x(3).^2 - x(1) * x(2) - x(1) * x(3);
grad_f = @(x) [ -2 * x(1) - x(2) - x(3), -4 * x(2) - x(1), -2 * x(3) - x(1) ];

h1 = @(x) x(1).^2 + x(2).^2 + x(3).^2 - 25;
grad_h1 = @(x) [ 2 * x(1), 2 * x(2), 2 * x(3) ];

h2 = @(x) 8 * x(1) + 14 * x(2) + 7 * x(3) - 56;
grad_h2 = @(x) [8, 14, 7];

pen_f = @(x, r) f(x) + (0.5 * r) * ( h1(x).^2 + h2(x).^2 );
grad_pen_f = @(x, r) grad_f(x) + r * ( grad_h1(x) + grad_h2(x) );

m_pf = @(x, r) min_pen_f(pen_f, grad_pen_f, x, r);

x0 = [25; 35; 45];
ret = BFGS_Pen(m_pf, 5000, x0, 10000);
disp(ret.x);
res = ret.x(end,:);
fprintf("BFGS_PEN returned x=%s with f(x)=%0.4f and h1(x)=%0.4f, h2(x)=%0.4f\n", vec2str(res), f(res'), h1(res'), h2(res'));

% Lösung mit fmincon
confun = @(x) confuneq(h1, h2, x);
ret = fmincon(f,x0,[],[],[],[],[],[],confun);
fprintf("fmincon returned x=%s with f(x)=%0.4f and h1(x)=%0.4f, h2(x)=%0.4f\n", vec2str(ret), f(ret'), h1(ret'), h2(ret'));

% Lösung mit fminunc bei festem r
pen_f_f = @(x) pen_f(x, 5000);
ret = fminunc(pen_f_f, x0);
fprintf("fminunc returned x=%s with f(x)=%0.4f and h1(x)=%0.4f, h2(x)=%0.4f\n", vec2str(ret), f(ret'), h1(ret'), h2(ret'));

% Gleichheitsbedingungen für fmincon
function [c,ceq] = confuneq(h1, h2, x)
    % Nonlinear inequality constraints
    c = [];
    % Nonlinear equality constraints
    ceq = [h1(x), h2(x)];
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