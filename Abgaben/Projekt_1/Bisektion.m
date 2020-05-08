% Function Name: Bisektion
%
% Description: Performs bisektion on intervals to search for a minimum
%
% Inputs:
%     Function, Start interval
% Outputs:
%   Minimum x
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     08.05.2020
%---------------------------------------------------------
function min_val = Bisektion(f, start_intrvl)

    k = 1;
    max_k = 100;

    a = start_intrvl(1);
    b = start_intrvl(2);
    
    % Grenze des neuen Intervalls
    x_mid = @(a, b) (a + b) / 2.0;
    
    % Abbruchkriterium: Anzahl an Iterationen
    while k <= max_k
        
        x = x_mid(a, b);
        fprintf("\tInterval [%0.4f, %0.4f] mit f(%0.4f) = %0.4f\n", a, b, x, f(x));

        if f(a) > f(b)
            % Neues Intervall [(a+b)/2, b]
            a = x;
        else
            % Neues Intervall [a, (a+b)/2]
            b = x;
        end
        
        k = k + 1;
    end
    
    % Ergebnis
    min_val = x_mid(a, b);
end
