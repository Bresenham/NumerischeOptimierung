% Fach: Numerische Optimierung
% Projekt 4
%
% Autor: Maximilian Gaul
% Date: 04.07.2020
%---------------------------------------------------------------

function ret = ActiveSet(Q, q, G, U_ext, b, r, x0)

    f = @(x) 0.5 * (x') * Q * x + (q') * x;
    gradient = @(x) Q * x + q;

    x = x0;
    k = 0;
    k_max = 10000;
    
    ret = struct("x", x, "f", f(x));
    
    [N, U] = getActiveIneqRestrictions(U_ext, x0, r);
    
    while k < k_max

        A = buildLeftMatrix(Q, G, U);

        right_side = zeros(size(A,1), 1);
        right_side(1:size(Q,1)) = -gradient(x);

        res = A \ right_side;

        [d, mu, lambda] = extractFromResult(res, size(Q, 1), size(G, 1), size(U, 1));
        
        if isEnd(d, lambda)
            return
        else
            [is_true, min_lambda_idx] = isPoint4(d, lambda, N);
            if is_true
                N(min_lambda_idx) = [];
                U(min_lambda_idx,:) = [];
            elseif isPoint5(d, x)
                x = x + d;
            elseif isPoint6(d, x)
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

function [t_min, t_min_idx] = calculateStepLengthPoint6(r, U, N, x, d)

    t_min = 1e45;
    t_min_idx = -1;
    
    for i = 1:size(U,1)
        if ~ismember(i, N) && norm( U(i, :) * d ) > 0
            t = ( r(i) - U(i, :) * x ) / ( U(i, :) * d );
            if t < t_min
                t_min = t;
                t_min_idx = i;
            end
        end
    end

end

function ret = isPoint6(d, x)

    ret = true;
    
    if norm( d ) < 1e-8
        ret = false;
    end
    
    d_x = x + d;
    
    for i = 1:length(d)
       if d_x(i) >= 0
           ret = false;
       end
    end
    
end

function ret = isPoint5(d, x)

    ret = true;
    
    if norm( d ) < 1e-8
        ret = false;
    end
    
    d_x = x + d;
    
    for i = 1:length(d)
       if d_x(i) < 0
           ret = false;
       end
    end

end

function [ret, min_lambda_idx] = isPoint4(d, lambda, N)

    ret = true;
    
    if norm( d ) > 1e-8
        ret = false;
    end
    
    min_lambda_idx = -1;
    min_lambda = 1e45;
    
    for k = 1:length(N)
        if lambda( N(k) ) < min_lambda
            min_lambda = lambda( N(k) );
            min_lambda_idx = k;
        end
    end
    
    if min_lambda >= 0
        ret = false;
    end

end
    

function ret = isEnd(d, lambda)

    ret = true;
    
    if norm( d ) > 1e-8
        ret = false;
    end
    
    for k = 1:length(lambda)
        if lambda(k) < 0
            ret = false;
        end
    end
    
end

function [d, mu, lambda] = extractFromResult(res, Q_row_len, G_row_len, U_row_len)

    d = res(1:Q_row_len);
    mu = res(Q_row_len+1:G_row_len);
    lambda = res(Q_row_len+G_row_len+1:Q_row_len+G_row_len+U_row_len);

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

function [N_ret, U_ret] = getActiveIneqRestrictions(U, x0, r)

    N_ret = [];
    U_ret = [];
    
    for i = 1:size(U, 1)
        if U(i, :) * x0 == r(i)
            N_ret = [N_ret, i];
            U_ret = [U_ret; U(i, :)];
        end
    end

end