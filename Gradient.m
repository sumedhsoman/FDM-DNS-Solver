%% Function to compute Gradient
function out = Gradient(phi,k,nx,ny,dx,dy,a1,a2,a3,a4)
if k == 1
    Agradx = Laplacian(((ny+2)*(nx+1)),dx,a1,a2);
    out = Agradx*reshape(Averaging_phi(phi,1,nx+2,ny+2),[],1);
elseif k == 2
    Agrady = Laplacian(((ny+1)*(nx+2)),dy,a3,a4);
    out = Agrady*reshape(Averaging_phi(phi,2,nx+2,ny+2),[],1); 
end
end