% Function Name: Bisektion
%
% Description: Performs bisektion on intervals to search for a minimum
%
% Inputs:
%     Function, Start interval
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     01.05.2020
%---------------------------------------------------------
function min_val = Bisektion(f, start_intrvl)

    k = 1;
    max_k = 10;

    a = start_intrvl(1);
    b = start_intrvl(2);
    
    x_eval = @(a, b) (a + b) / 2.0;
    l_eval = @(a, b) ( a + x_eval(a, b) ) / 2.0;
    r_eval = @(a, b) ( x_eval(a, b) + b ) / 2.0;
    
    while k <= max_k
        
        x = x_eval(a, b);
        % fprintf("\tInterval [%0.4f, %0.4f] mit f(%0.4f) = %0.4f\n", a, b, x, f(x));
        if abs( f(x) - 0.56714) <= 1e-6
            break;
        end
        if f( l_eval(a, b) ) < f( r_eval(a, b) )
            a = a;
            b = x;
        else
            a = x;
            b = b;
        end
        
        k = k + 1;
    end
    
    min_val = x_eval(a, b);
end
