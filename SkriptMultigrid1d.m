function SkriptMultigrid1d
    N = 64;
    l = 4;
    x = linspace(0,1,N+1);
    A = LEVELM(l,N);
    b = loadvector(x,@f);

    [u,resmax] = MGM1d(l,A,zeros(N-1,1),b,1e-6, 100, 1);
    [ug,resmaxg] = GaussSeidel(A{l},b,zeros(N-1,1), 100);

    echt = sin(pi*x(1:N-1))/(2*pi^2) + sin(16*pi*x(1:N-1))/(512*pi^2);
    draw(x, u, ug, echt, resmax, resmaxg, N);

end

function val = f(x)
    %val = -1;
    %val = ((16*pi^2 - 4)*sin(x*4*pi) - 16*pi*cos(x*4*pi)).*exp(2*x);
    val = (sin(pi*x) + sin(16*pi*x))/2;
end