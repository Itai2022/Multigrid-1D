function [A4lvl] = LEVELM(l, N)
    A4lvl = cell(l,1);
    for i = l:-1:1
        A4lvl{i} = N*(diag(-1*ones(N-2,1), -1) + ...
                   diag(2*ones(N-1,1)) + ...
                   diag(-1*ones(N - 2,1),1)); 
        N = N / 2;
    end
end