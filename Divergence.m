%% Function to compute divergence of the velocities, for the Poisson's equation
function div = Divergence(a1,a2,a3,a4,ua,va,nx,ny,dx,dy)
Adivx = Laplacian((nx+2)*(ny+2),dx,a1,a2);
divx = Adivx*AveragingOperator2(ua,1,nx+1,ny+2);
Adivy = Laplacian((ny+2)*(nx+2),dy,a3,a4);
divy = Adivy*AveragingOperator2(va,2,nx+2,ny+1);
div = divx+divy;
end