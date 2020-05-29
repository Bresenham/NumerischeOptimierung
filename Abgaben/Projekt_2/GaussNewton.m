% Function Name: GaussNewton
%
% Description: GaussNewton-Method for least-squares calcuation
%
% Inputs:
%   function: The function for which the distance should be minimized
%   Startvalue: x0,
%   Data to fit: xdata, ydata
% Outputs:
%   a vector of structs [ x, residuum(x), jacobi(x) ]
%
% WolfePowell is defined in 'WolfePowell.m'
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     23.05.2020
%---------------------------------------------------------
function ret = GaussNewton(f, f_partials, x0, xdata, ydata)

    k = 0;
    x = x0;
    kmax = 1e+6;
    
    f_jacobi = @(x, xdata, ydata) jacobi(f_partials, x, xdata);
    f_resid = @(x, xdata, ydata) residuum(f, x, xdata, ydata);
    f_lq_sum = @(x) sum( f_resid(x, xdata, ydata).^2 );
    
    gradient = @(x) 2 * f_jacobi(x, xdata, ydata)' * f_resid(x, xdata, ydata);
    
    ret = struct( "x", x, "residuum", f_resid(x, xdata, ydata), "jacobi", f_jacobi(x, xdata, ydata) );
    
    while norm( gradient(x) ) > 1e-8 && k < kmax
        
        j = f_jacobi(x, xdata, ydata);
        resid = f_resid(x, xdata, ydata);
        
        right_side = -j' * resid;
        left_side = j' * j;
        
        d = left_side \ right_side;
        
        phi = @(a) f_lq_sum(x + a * d);
        phi_grad = @(a) gradient(x + a * d)' * d;
        alpha = WolfePowell(phi, phi_grad);
        
        x = x + 0.125 * d;
        
        k = k + 1;
        
        ret = [ ret; struct( "x", x, "residuum", f_resid(x, xdata, ydata), "jacobi", f_jacobi(x, xdata, ydata) ) ];
    end
    
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