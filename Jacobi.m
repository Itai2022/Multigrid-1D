%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gedämpfte Jacobiverfahren zur Lösung von Ax = b
%% A die Systemmatrix, b die rechte Seite
%% x_old Anfangswert tol die Abbruchskriterium 
%% maxIter die maximale Anzahl von Iterationen
%% w ist Dämpfungsparameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function [x_new, res0, res10, res100, iter] = Jacobi(A, b, x_old, tol, maxIter, w)
function [x_new] = Jacobi(A, b, x_old,maxIter, w)
    len = length(b);
    x_new = zeros(len,1);
    for iter = 1:maxIter
        for k  = 1:len
            sigma = A(k, 1:k-1)*x_old(1:k-1) + A(k, k+1:len)*x_old(k+1:len);       
            x_new(k) = (1 - w) * x_old(k) + (w/A(k,k))*(b(k) - sigma);
        end
        x_old = x_new;
    end
end