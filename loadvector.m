%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aufstellen des Loadvektors
%% Benutze die Quadratur Formel:
%% int_a^b f(x) dx = b-a/2(f((a+b)/2 - (b-a)/2sqrt(3)) + f((a+b)/2 + (b-a)/2sqrt(3)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [b] = loadvector(x, f)
    N = length(x);
    b = zeros(N-2, 1);
    h = x(2) - x(1);
    for i = 2:N-1
        w1 = (x(i-1)+x(i))/2;
        w2 = h/(2*sqrt(3));
        sum1 = f(w1-w2)*(w1-w2-x(i-1))/h;
        sum2 = f(w1+w2)*(w1+w2-x(i-1))/h;
        w3 = (x(i)+x(i+1))/2;
        sum3 = f(w3-w2)*(x(i+1)-w3+w2)/h;
        sum4 = f(w3+w2)*(x(i+1)-w3-w2)/h;
        b(i-1) = (sum1+sum2+sum3+sum4)*h/2;
    end

end

