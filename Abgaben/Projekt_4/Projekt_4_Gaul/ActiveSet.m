% Fach: Numerische Optimierung
% Projekt 4
%
% Autor: Maximilian Gaul
% Date: 04.07.2020
%---------------------------------------------------------------

function ret = ActiveSet(Q, q, G, U, b, r, x0)

    f = @(x) 0.5 * (x') * Q * x + (q') * x;
    gradient = @(x) Q * x + q;

    x = x0;
    k = 0;
    k_max = 1000;
    
    ret = struct("x", x, "f", f(x));
    
    N = getActiveIneqRestrictions(U, x0, r);
    
    while k < k_max

        A = buildLeftMatrix( Q, G, U(N, 1:end) );

        right_side = zeros(size(A,1), 1);
        right_side(1:size(Q,1)) = -gradient(x);

        res = A \ right_side;

        [d, mu, lambda] = extractFromResult(res, size(Q, 1), size(G, 1), N, size(U, 1));
        
        if isEnd( d, lambda(N) )
            return
        else
            [is_true, min_lambda_idx] = isPoint4(d, lambda, N);
            if is_true
                N(min_lambda_idx) = [];
            elseif isPoint5(d, x, U, r, G, b)
                x = x + d;
            elseif isPoint6(d, x, U, r, G, b)
                [t_min, t_min_idx] = calculateStepLengthPoint6(r, U, N, x, d);
                if t_min_idx > -1
                    x = x + t_min * d;
                    N = [N, t_min_idx];
                else
                    disp("ERROR: Invalid descent step length!");
                    return;
                end
            end
            
            k = k + 1;
        end

        ret = [ ret; struct("x", x, "f", f(x)) ];
        
    end

end

% Berechnet minimale Schrittweite falls Punkt 6 im Active-Set Algorithmus erfüllt ist
function [t_min, t_min_idx] = calculateStepLengthPoint6(r, U, N, x, d)

    t_min = 1e45;
    t_min_idx = -1;
    
    for i = 1:size(U,1)
        if ~ismember(i, N) && U(i, :) * d > 0
            t = ( r(i) - U(i, :) * x ) / ( U(i, :) * d );
            if t < t_min
                t_min = t;
                t_min_idx = i;
            end
        end
    end

end

% Prüft, ob es sich im aktuellen Durchlauf um Fall 6 im Active-Set Algorithmus handelt
function ret = isPoint6(d, x, U, r, G, b)

    d_not_zero = norm( d ) > 1e-8;
    
    d_x = x + d;
    
    has_wrong_ineq_constr = false;
    for i = 1:size(U, 1)
        if U(i, :) * d_x > r(i)
            has_wrong_ineq_constr = true;
        end
    end
    
    has_wrong_eq_constr = false;
    for k = 1:size(G, 1)
        % Floating-Point Werte auf Gleichheit überprüfen
        if ~( abs( G(k, :) * d_x - b(k) ) < 1e-8 )
            has_wrong_eq_constr = true;
        end
    end
    
    ret = d_not_zero && ( has_wrong_ineq_constr || has_wrong_eq_constr );
    
end

% Prüft, ob es sich im aktuellen Durchlauf um Fall 5 im Active-Set Algorithmus handelt
function ret = isPoint5(d, x, U, r, G, b)
    
    d_not_zero = norm( d ) > 1e-8;
    
    d_x = x + d;
    
    has_wrong_ineq_constr = false;
    for i = 1:size(U, 1)
        if U(i, :) * d_x > r(i)
            has_wrong_ineq_constr = true;
        end
    end
    
    has_wrong_eq_constr = false;
    for k = 1:size(G, 1)
        % Floating-Point Werte auf Gleichheit überprüfen
        if ~( abs( G(k, :) * d_x - b(k) ) < 1e-8 )
            has_wrong_eq_constr = true;
        end
    end
    
    ret = d_not_zero && ~( has_wrong_ineq_constr || has_wrong_eq_constr );
end

% Prüft, ob es sich im aktuellen Durchlauf um Fall 4 im Active-Set Algorithmus handelt
function [ret, min_lambda_idx] = isPoint4(d, lambda, N)

    d_is_zero = norm( d ) <= 1e-8;
    
    min_lambda_idx = -1;
    min_lambda = 1e45;
    
    for k = 1:length(N)
        if lambda( N(k) ) < min_lambda
            min_lambda = lambda( N(k) );
            min_lambda_idx = k;
        end
    end
    
    ret = d_is_zero && ( min_lambda < 0 );

end
    
% Prüft, ob der Algorithmus fertig ist
function ret = isEnd(d, lambda)

    d_is_zero = norm( d ) <= 1e-8;
    
    has_lambda_below_zero = false;
    for k = 1:length(lambda)
        % Floating-Point Rundungen auf 0 überprüfen
        if lambda(k) < (-1e-8)
            has_lambda_below_zero = true;
        end
    end
    
    ret = d_is_zero && ~( has_lambda_below_zero );
    
end

% Liest aus dem gelösten Gleichungssystem Werte für d, mu und lambda
function [d, mu, lambda_full] = extractFromResult(res, Q_row_len, G_row_len, N, U_row_len)

    d = res(1:Q_row_len);
    mu = res(Q_row_len+1:Q_row_len+G_row_len);
    lambda_full = zeros(U_row_len, 1) - 10 * ones(U_row_len, 1) ;
    lambdas = res(Q_row_len+G_row_len+1:end);
    
    for i = 1:length(lambdas)
        lambda_full( N(i) ) = lambdas(i);
    end

end

% Große, linke Matrix im zu lösenden Gleichungssystem
function A = buildLeftMatrix(Q, G, U)

    [Q_row_len, Q_col_len] = size(Q);
    [G_row_len, G_col_len] = size(G);
    [G_t_row_len, G_t_col_len] = size(G');
    [U_row_len, U_col_len] = size(U);
    [U_t_row_len, U_t_col_len] = size(U');
    
    edge_len = Q_col_len + G_t_col_len + U_t_col_len;
    A = zeros(edge_len, edge_len);
    
    A( 1:Q_col_len, 1:Q_col_len ) = Q;
    A( 1:G_t_row_len, Q_col_len+1:(Q_col_len+G_t_col_len) ) = G';
    A(1:U_t_row_len, Q_col_len+G_t_col_len+1:Q_col_len+G_t_col_len+U_t_col_len) = U';
    
    A(Q_row_len+1:Q_row_len+G_row_len, 1:G_col_len) = G;
    A(Q_row_len+G_row_len+1:Q_row_len+G_row_len+U_row_len, 1:U_col_len) = U;

end

% Liefert die Menge an aktiven Ungleichungsnebenbedingungen
function N_ret = getActiveIneqRestrictions(U, x0, r)

    N_ret = [];
    
    for i = 1:size(U, 1)
		% Floating-Point Rundungen auf Gleichheit überprüfen
        if abs( U(i, :) * x0 - r(i) ) < 1e-8
            N_ret = [N_ret, i];
        end
    end

end