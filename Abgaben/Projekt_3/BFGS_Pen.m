% Fach: Numerische Optimierung
% Projekt 3 - Penalty-Verfahren & SQP-Verfahren
%
% Autor: Maximilian Gaul
% Date: 12.06.2020
%---------------------------------------------------------------

function res = BFGS_Pen(fun,r,x0,maxit)

    % BFGS-Verfahren
    % Parameter:
    % fun: function handle
    % r: Penalty-Parameter
    % x0: Startwert
    % maxit: Maximale Anzahl der Iterationen

    glob = 1;
    x(1,:)=x0;x(2,:)=inf*ones(1,length(x0));
    k=1;sig=1;rho=0.01;
    B=eye(length(x));
    [f,g]=fun(x(k,:)',r);
    grad(k,:)=g';
    while (k<maxit && (k==1||norm(g)>1e-6))
        d = B*g';  % Suchrichtung
        % Falls Abstiegsrichtung die Winkelbedingung nicht erfuellt, verwende
        % negativen Gradienten als Abstiegsrichtung.
        if glob % Falls globale Variante verwendet, bestimme Schrittweite
            if ( grad(k,:) * d ) < ( rho * norm( grad(k,:) ) * norm(d) )
                d=grad(k,:)';
                disp(['Negativer Gradient in Schritt ',num2str(k),' verwendet.']);
            end
            % Aufruf der Armijo-Schrittweitensteuerung (mit Extra-Parameter
            % fuer Penalty-Funktion
            sig = ArmijoPen(fun,r,x(k,:)',-d,sig);
            sigma(k)=sig;
        end
        x(k+1,:)=x(k,:)-sig*d';  % Neue Naeherung
        fv(k)=f;
        [f,g]=fun(x(k+1,:)',r); % Neue Werte von Funktion und Gradient
        grad(k+1,:)=g';
        % Update der Hesse-Approximation
        B = local_update_matrix(B,x(k,:),x(k+1,:),grad(k,:),grad(k+1,:));
        k=k+1;
    end
    res.x = x;
    res.fv = fv';
    res.grad = grad;
    if glob
        res.sigma = sigma;
    end
    
end


function res = local_update_matrix(mat,x_old,x_new,grad_old,grad_new)
    
    y = (grad_new - grad_old)';
    s = (x_new - x_old)';
    bracket = s - mat*y;
    res = mat + (bracket*s'+s*bracket')/(y'*s) - (((bracket'*y)*s)*s')/((y'*s)^2);

end