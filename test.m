N = 64;
x_cor = linspace(0,1,N+1);
A = N*(diag(-1*ones(N,1), -1) + diag(2*ones(N+1,1))+ diag(-1*ones(N,1),1)); 
A(N+1, N+1) = 1;
A(N+1,N) = 0;
A(1,1) = 1;
A(1,2) = 0;
%b = -1/N*ones(N-1,1);
b = zeros(N+1,1);
b(2:N) = loadvector(x_cor,@f);
% b(1) = 0;
% b(end) = 0;
x0 = sin(x_cor*10*pi);
%x = zeros(1, N+1);
% x0(1) = 0;
% x0(end) = 0;
[x,r0,r10,r100] = Jacobi(A,b,x0',1e-6,100, 1);
%[x,r0,r10,r100] = GS(A,b,x0',100);
%x = A\b;



plot(x_cor', r0, "x-", x_cor', r10, "*-", x_cor', r100,"o-");
legend('0 Iteration','10 Iterationen','100 iterationen')
xlabel("Gitterpunkte")
ylabel("Residum")
title("Residuumskurve gegen die Gitterpunkte")
%plot(x_cor',x);

function val = f(x)
    %val = -1;
    %val = ((16*pi^2 -4)*sin(x*4*pi) - 16*pi*cos(x*4*pi))*exp(2*x);
    val = (sin(pi*x) + sin(16*pi*x))/2;
end
