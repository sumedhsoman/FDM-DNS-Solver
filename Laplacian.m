function A = Laplacian(n,h,a1,a2)
% Status - COMPLETE
% a1/a2: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a1 0;ones(n-2,1)*[-1 2 -1];0 a2 -1],-1:1,n,n)'/h^2;    
end