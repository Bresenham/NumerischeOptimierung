% Function Name: Mutation
%
% Description: Performs Mutation-Selection on a function to find a
% minimum
%
% Inputs:
%   Function, Start value
% Outputs:
%   Minimum x, Iterations
%
% Version:	MATLAB R2020a
% Author:	Maximilian Gaul
% Date:     01.05.2020
%---------------------------------------------------------
function [min_val, k] = Mutation(f, start_x)
    
    k = 1;
    
    % Anzahl an Iterationen nachdem der Algorithmus automatisch abbricht
    k_max = 15000;
    
    % Schrittweite
    alpha = 0.25;
    
    dim = numel(start_x);
    
    % Erzeugt einen Zufallsvektor der Dimension 'dim' mit Werten
    % von [from, to]
    rand_vec = @(dim, from, to) from + (to - from) * rand(1, dim);
    
    % Berechnet einen neuen x-Wert auf Basis des aktuellen Wertes
    % und einer Schrittweite 'alpha' sowie eines Zufallsvektors
    hat = @(x, alpha) x + alpha * rand_vec(dim, -0.5, 0.5);
    
    x = start_x;
    
    % Startwerte für die Berechnung des Differenzbetrags damit die
    % while-Schleife mindestens 1x ausgeführt wird
    f_x_val = 1000;
    f_x_old_val = -1000;

    % Abbruchkriterien: Differenzbetrag oder Anzahl Iterationen
    while norm(f_x_val - f_x_old_val) > 1e-6
    % while k <= k_max
    
        x_hat = hat(x, alpha);
        
        f_x = f(x);
        f_x_hat = f(x_hat);
        
        if f_x_hat < f_x
            % Es wurde ein besserer x-Wert gefunden
            %   - Variablen für den Differenzbetrag anpassen
            %   - x anpassen
            f_x_old_val = f_x;
            f_x_val = f_x_hat;
            x = x_hat;
            
            % Es kann sinnvoller sein, die Ausgabe nur bei einem besseren
            % x-Wert zu tätigen
            fprintf("\tx = [ %s] mit f(x) = %0.8f\n", sprintf("%0.4f ", x), f(x));
        end
        
        % fprintf("\tx = [ %s] mit f(x) = %0.8f\n", sprintf("%0.4f ", x), f(x));
        
        % Verstärkungsfaktor wie in der Ausarbeitung beschrieben,
        % wird bis zu einer gewissen Grenze auf alpha multipliziert
        % alpha = min([alpha * 1.05, 0.5]);
        
        k = k + 1;
    end
    
    % Rückgabewert Minimum
    min_val = x;
end
