%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mehrgitterverfahren in 1 Dimension
%% A die Systemmatrix, b die rechte Seite
%% u_0 Anfangswert tol die Abbruchskriterium 
%% maxIter die maximale Anzahl von Iterationen
%% V_Cycle: mu = 1; W_Cycle: mu = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_0, resmax] = MGM1d(l, A4lvl, u_0, b, tol, maxIte, mu)
    resmax = zeros(maxIte,1);
    normb = norm(b);
    resmax(1) = norm(b - A4lvl{l}*u_0)/normb;
    for ite = 1:maxIte
        u_0 = V_Cycle(l, A4lvl, u_0, b, mu);
        resmax(ite+1) = norm(b - A4lvl{l}*u_0)/normb; 
        if resmax(ite+1) < tol
            break;
        end
    end
    if ite == maxIte
        fprintf("Keine Konvergenz");
    end
end

function u0 = V_Cycle(l, A4lvl, u0, b, mu)
    % 2-mal Vorglättung mit GaussSeidel
    u0 = gaussSeidel(A4lvl{l},b,u0,2);

    %Grob-Gitter_korrektur
    d1 = b - A4lvl{l}*u0;  % Berechnung der Defekt.
    d0 = Restriction(d1);  % Restriktion des Defekts.
    if l == 2
        v = A4lvl{1}\d0;
    else
        for i = 1:mu
            v = V_Cycle(l-1,A4lvl,zeros(length(d0),1),d0, 1);
        end
    end
    u0 = u0 + Prolongation(v);  % Prolongation der Fehler.

    % 2-mal A-Posteriori-Glättung mit GaussSeidel
    u0 = gaussSeidel(A4lvl{l},b,u0,2);
end

function b = Restriction(d)
    N = (length(d) - 1)/2;
    b = zeros(N,1);
    for i = 1:N
        b(i) = (d(2*i-1) + 2*d(2*i) + d(2*i+1))/2;
    end
end

function v = Prolongation(d)
    N = length(d)*2 + 1;
    v = zeros(N,1);
    v(2:2:N) = d(1:1:length(d));
    for i = 3 : 2 : N-1
        v(i) = (d((i+1)/2) + d((i - 1)/2))/2;
    end
    v(1) = d(1)/2;
    v(end) = d(end)/2;
end
function x1 = gaussSeidel(A,b,x0,Itmax)
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    G = (D - L)\U;
    f = (D-L)\b;
    for i = 1:Itmax
        x1 = G*x0 + f;
        x0 = x1;
    end
end
