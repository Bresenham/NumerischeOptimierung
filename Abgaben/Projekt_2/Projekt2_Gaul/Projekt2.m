% Fach: Numerische Optimierung
% Projekt 2 - Quasi-Newton-Verfahren & Gauß-Newton-Verfahren
%
% Autor: Maximilian Gaul
% Date: 22.05.2020
%---------------------------------------------------------------

%--------------------Allgemeine Definitionen--------------------

% Funktion von Himmelblau
f_himmel = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
f_himmel_grad = @(x) [ 2 * (x(1).^2 + x(2) - 11) * 2 * x(1) + 2 * (x(1) + x(2).^2 - 7);
                    2 * (x(1).^2 + x(2) - 11) + 2 * (x(1) + x(2).^2 - 7) * 2 * x(2) ];

% Bazaraa-Shetty-Funktion
f_bazaraa = @(x) 100 * (x(1) - 2).^4 + (x(1) - 2 * x(2)).^2;
f_bazaraa_grad = @(x) [ 2 * ( 200 * (x(1) - 2).^3 + x(1) - 2 * x(2) ); 8 * x(2) - 4 * x(1) ];

% 2D Rosenbrock-Funktion
f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [ 400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2) ];

% N-dimensionale Rosenbrock-Funktion
f_rosen_mult = @(x) sum( 100 * ( x(2:length(x))-x(1:(length(x)-1)).^2 ).^2 + ( 1 - x(1:( length(x)-1) ) ).^2 );
f_rosen_mult_grad = @(x) f_rosen_mult_deriv_func(x);

%---------------------------------------------------------------

% Aufgabe 1
% Siehe Datei 'InverseBFGS.m'

% Aufgabe 2
% Siehe Erläuterung im PDF

% Aufgabe 3
% Siehe Erläuterung im PDF

% Aufgabe 4

fprintf("--------------------AUFGABE 4--------------------\n");

x0 = [0; -1];

ret = InverseBFGS(f_himmel, f_himmel_grad, x0);
fprintf("InverseBFGS returned x=%s with f_himmel(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f); 

ret = InverseBFGS(f_bazaraa, f_bazaraa_grad, x0);
fprintf("InverseBFGS returned x=%s with f_bazaraa(x)=%0.6f\n", vec2str(ret(end).x), ret(end).f);

ret = InverseBFGS(f_rosen, f_rosen_grad, x0);
fprintf("InverseBFGS returned x=%s with f_rosen(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f);

rosenbrock_dim = 150;
x0_rosen = zeros(rosenbrock_dim, 1) - ones(rosenbrock_dim, 1);
ret = InverseBFGS(f_rosen_mult, f_rosen_mult_grad, x0_rosen);
fprintf("InverseBFGS returned x=%s with f_rosen_mult(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f);

% Setze die selben Toleranzen und Grenzen wie in 'InverseBFGS' und
% zusätzlich das BFGS-Verfahren zum Updaten der Hesse-Matrix
options = optimoptions("fminunc", "OptimalityTolerance", 1e-8, "MaxFunctionEvaluations", 1e+8, "MaxIterations", 1e+6, "Algorithm", "quasi-newton", "HessUpdate", "bfgs");

ret = fminunc(f_himmel, x0, options);
fprintf("fminunc returned x=%s with f_himmel(x)=%0.6f\n", vec2str(ret), f_himmel(ret));

ret = fminunc(f_rosen_mult, x0_rosen, options);
fprintf("fminunc returned x=%s with f_rosen_mult(x)=%0.6f\n", vec2str(ret), f_rosen_mult(ret));

% Aufgabe 5
% Siehe Erläuterung im PDF

% Aufgabe 6
% Siehe GaussNewton.m

% Aufgabe 7
fprintf("--------------------AUFGABE 7--------------------\n");

f = @(t, x) x(1) * exp( x(2) * t );
f_partial_x1 = @(t, x) exp( x(2) * t);
f_partial_x2 = @(t, x) t * x(1) * exp( x(2) * t );

f_x0 = [0.5; 1];
f_xdata = [0; 1; 2; 3];
f_ydata = [2.0; 0.7; 0.3; 0.1];

ret = GaussNewton(f, {f_partial_x1, f_partial_x2}, f_x0, f_xdata, f_ydata);
fprintf("GaussNewton returned x=%s for function f after %d steps\n", vec2str(ret(end).x), length(ret));

g = @(t, x) x(1) * exp( -(x(2).^2 + x(3).^2)*t ) * ( sinh( x(3).^2 * t ) / ( x(3).^2 ) );
g_partial_x1 = @(t, x) ( exp( -t * (x(2).^2 + x(3).^2) ) * sinh( x(3).^2 * t) ) / ( x(3).^2 );
g_partial_x2 = @(t, x) - ( 2 * t * x(1) * x(2) * exp( -t * (x(2).^2 + x(3).^2) ) * sinh( x(3).^2 * t ) ) / ( x(3).^2 );
g_partial_x3 = @(t, x) ( 2 * x(1) * exp( -t * (x(2).^2 + x(3).^2) ) * ( t * x(3).^2 * cosh(t * x(3).^2) - (t * x(3).^2 + 1) * sinh(t * x(3).^2) ) ) / ( x(3).^2 );

g_x0 = [10; 0.05; 0.1];
g_x0_2 = [5; 0.145; 0.125];
g_x0_3 = [3; 0.1; 0.05];
g_xdata = [6;12;18;24;30;36;42;48;54;60;66;72;78;84;90;96;102;108;114;120;126;132;138;144;150;156;162;168;174;180];
g_ydata = [24.19;35.34;43.43;42.63;49.92;51.53;57.39;59.56;55.60;51.91;58.27;62.99;52.99;53.83;59.37;62.35;61.84;61.62;49.64;57.81;54.79;50.38;43.85;45.16;46.72;40.68;35.14;45.47;42.40;55.21];

ret = GaussNewton(g, {g_partial_x1, g_partial_x2, g_partial_x3}, g_x0, g_xdata, g_ydata);
fprintf("GaussNewton returned x=%s for function g and start value x0=%s after %d steps\n", vec2str(ret(end).x), vec2str(g_x0), length(ret));

ret = GaussNewton(g, {g_partial_x1, g_partial_x2, g_partial_x3}, g_x0_2, g_xdata, g_ydata);
fprintf("GaussNewton returned x=%s for function g and start value x0=%s after %d steps\n", vec2str(ret(end).x), vec2str(g_x0_2), length(ret));

% Aufgabe 8
% Siehe Erläuterung im PDF

% Aufgabe 9
fprintf("--------------------AUFGABE 9--------------------\n");

f_lq_sum = @(x) sum( residuum(f, x, f_xdata, f_ydata).^2 );
f_lq_grad = @(x) 2 * jacobi({f_partial_x1, f_partial_x2}, x, f_xdata)' * residuum(f, x, f_xdata, f_ydata);
ret = InverseBFGS(f_lq_sum, f_lq_grad, [0.5; 1]);
fprintf("InverseBFGS returned x=%s for function f after %d steps\n", vec2str(ret(end).x), length(ret));

g_lq_sum = @(x) sum( residuum(g, x, g_xdata, g_ydata).^2 );
% g_lq_grad = @(x) 2 * ( jacobi({g_partial_x1, g_partial_x2, g_partial_x3}, x, g_xdata)' * residuum(g, x, g_xdata, g_ydata) );
g_lq_grad = @(x) [  ( g_lq_sum([x(1) + 1e-10; x(2); x(3)]) - g_lq_sum(x) ) * (1e10);
                    ( g_lq_sum([x(1); x(2) + 1e-10; x(3)]) - g_lq_sum(x) ) * (1e10);
                    ( g_lq_sum([x(1); x(2); x(3)+ 1e-10]) - g_lq_sum(x) ) * (1e10);
                ];
ret = InverseBFGS(g_lq_sum, g_lq_grad, g_x0_2);
fprintf("InverseBFGS returned x=%s for function g and start value x0=%s after %d steps\n", vec2str(ret(end).x), vec2str(g_x0_2), length(ret));

ret = InverseBFGS(g_lq_sum, g_lq_grad, g_x0_3);
fprintf("InverseBFGS returned x=%s for function g and start value x0=%s after %d steps\n", vec2str(ret(end).x), vec2str(g_x0_3), length(ret));

% Berechnet den Gradienten der N-dimensionalen Rosenbrock-Funktion
% \input: Vector 'x'
% \output: Gradient an der Stelle 'x'
function ret = f_rosen_mult_deriv_func(x)

    start = @(x1, x2) -2 * (1 - x1) + 200 * (x2 - x1.^2) * (-2 * x1);
    middle = @(x1, x2, x3) 200 * (x2 - x1.^2) - 2 * (1 - x2) + 200 * (x3 -x2.^2) * (-2 * x2);
    ending = @(x2, x3) 200 * (x3 - x2.^2);
    
    dim = numel(x);
    
    sum = zeros(dim, 1);
    
    for i = 1:dim
        if i == 1
            sum(i) = sum(i) + start( x(i), x(i + 1) );
        elseif i == dim
            sum(i) = sum(i) + ending( x(i - 1), x(i) );
        else
            sum(i) = sum(i) + middle( x(i - 1), x(i), x(i + 1) );
        end
    end
    
    ret = sum;

end

% Berechnet die Jacobi-Matrix von f für den gegebenen Datensatz
% \input: Funktion, partielle Ableitungen, Wert an dem die Funktion
% ausgewertet werden soll, Datensatz
% \output: Jacobi-Matrix
function ret = jacobi(f_partials, x, xdata)
    
    ydim = numel(xdata);
    xdim = numel(f_partials);

    r = zeros(ydim, xdim);
    for d = 1:ydim
        for i = 1:xdim
            r(d,i) = f_partials{i}(xdata(d), x);
        end
    end
    
    ret = r;
end

% Berechnet das Residuum von f für den gegebenen Datensatz
% \input: Funktion, Wert an dem die Funktion ausgewertet werden soll,
% Datensatz
% \output: Residuum
function ret = residuum(f, x, xdata, ydata)

    xdim = numel(xdata);

    r = zeros(xdim, 1);
    for i = 1:xdim
        r(i) = f(xdata(i), x) - ydata(i);
    end
    
    ret = r;

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