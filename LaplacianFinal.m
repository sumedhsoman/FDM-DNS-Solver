function Coeff = LaplacianFinal(a1,a2,a3,a4,nx,ny,dx,dy)
% Status-COMPLETE   
%Function to generate Laplacian matrix, with 5 point stencil, and
%modifications for the central term, probably adding a diagonal matrix?
% a1 a2 a3 a4 are for denoting the boundary conditions.
  Coeff = kron(speye(ny),Laplacian(nx,dx,a1,a2))+kron(Laplacian(ny,dy,a3,a4),speye(nx));
end