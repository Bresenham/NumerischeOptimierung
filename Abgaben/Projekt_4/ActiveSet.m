% Fach: Numerische Optimierung
% Projekt 4
%
% Autor: Maximilian Gaul
% Date: 04.07.2020
%---------------------------------------------------------------

function ret = ActiveSet(Q, q, G, U, b, r, x0)

    x = x0;
    
    % Aktive Ungleichungsnebenbedingungen
    N = getActiveIneqRestrictions(U, x0, r);
    A = buildLeftMatrix(Q, G, U);
    disp(A);

end

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

function N = getActiveIneqRestrictions(U, x0, r)

    N = [];
    for i = 1:size(U, 1)
        if U(i, :) * x0 == r(i)
            N = [N, i];
        end
    end

end