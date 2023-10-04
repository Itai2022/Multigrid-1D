function draw(x,u,ug,echt,resmax,resmaxg,N)    
    figure;
    plot(x(1:N-1), u, x(1:N-1), ug,x(1:N-1), echt);
    legend("Multigrid", "GaussSeidel", "Echte Lösung");
    xlabel("Giterpunkte");
    ylabel("Nährungslösng");
    
    figure;
    plot((0:length(find(resmax)) - 1)', resmax(resmax>0), "x-",(0:length(find(resmaxg)) - 1)', resmaxg, "+-");
    legend("Multigrid", "Gauss-Seidel")
    xlabel("Anzahl den Iterationen");
    ylabel("Maximaler Wert im Residuum");
end