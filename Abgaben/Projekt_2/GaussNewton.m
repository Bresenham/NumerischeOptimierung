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
function ret = GaussNewton(f, f_resid, f_jacobi, x0, xdata, ydata)

    k = 0;
    x = x0;
    kmax = 1e+6;
    
    ret = struct( "x", x, "residuum", f_resid(x, xdata, ydata), "jacobi", f_jacobi(x, xdata, ydata) );
    
    gradient = @(x) 2 * f_jacobi(x, xdata, ydata)' * f_resid(x, xdata, ydata);
    
    while norm( gradient(x) ) > 1e-8 && k < kmax
        
        jacobi = f_jacobi(x, xdata, ydata);
        resid = f_resid(x, xdata, ydata);
        
        right_side = jacobi' * resid;
        left_side = jacobi' * jacobi;
        
        d = left_side \ right_side;
        
        phi = @(a) f_resid(x + a * d, xdata, ydata)' * f_resid(x + a * d, xdata, ydata);
        phi_grad = @(a) gradient(x + a * d)' * d;
        alpha = WolfePowell(phi, phi_grad);
        
        x = x - alpha * d;
        fprintf("NORM: %0.4f\n", norm( gradient(x) ) );
        
        k = k + 1;
        
        ret = [ ret; struct( "x", x, "residuum", resid, "jacobi", jacobi ) ];
    end
    
end