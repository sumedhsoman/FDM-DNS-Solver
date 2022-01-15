%% Defining Coefficients
alpha1 = 4/15;
alpha2 = 1/15;
alpha3 = 1/6;
gamma1 = 8/15;
gamma2 = 5/12;
gamma3 = 3/4;
zeta1 = -17/60;
zeta2 = -5/12;
Re = 30;
nx = 30;
ny = 30;
dx = 0.01;
dy = 0.01;
dt = 1e-05;
t = 0.001;
u = reshape(ones(ny+2,nx+1),[],1);
v = reshape(ones(ny+1,nx+2),[],1);
U = 1;
CFL = zeros(t/dt,1);
%% Modified Runge-Kutta Algorithm
for i = 0:dt:t
    %i % Display time step
%% Predictor Step 1 
%% Velocity x

LHSx = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))));
RHSx = dt*(gamma1*ConvectionOperator(u,v,1,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))*u);
ua = RHSx\(LHSx);
ua = ua';

%% Predictor Step 1 Velocity y

LHSy = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))));
RHSy = dt*(gamma1*ConvectionOperator(u,v,2,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy))*v);
va = RHSy\LHSy;
va = va';

%% Poisson's Equation

phi = (Divergence(1,1,1,1,ua,va,nx,ny,dx,dy))\(gamma1*(LaplacianFinal(1,1,1,1,nx+2,ny+2,dx,dy)));
phi = phi';

%% Apply phi boundary conditions

phi([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phi(2:nx+1) = phi(2+nx+2:2*(nx+2)-1);
phi(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1) = phi(2+(nx+2)*(ny):(nx+2)*(ny+1)-1);
phi(nx+3:nx+2:1+(nx+2)*(ny)) = phi(nx+4:nx+2:2+(nx+2)*(ny));
phi(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phi(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);

%% Update Velocities
utilde = ua-dt*(gamma1*Gradient(phi,1,nx,ny,dx,dy,1,1,1,1));
vtilde = va-dt*(gamma1*Gradient(phi,2,nx,ny,dx,dy,1,1,1,1));
%% Update intermediate BCs

% Firstly enforce invalid cell velocities as zero.
utilde([1 nx+1 1+(ny+1)*(nx+1) (nx+1)*(ny+2)],:) = 0;
vtilde([1 nx+2 1+(ny)*(nx+2) (nx+2)*(ny+1)],:) = 0;
% BC for u
utilde(2:nx,:) = 2*U-utilde(2+nx+1:nx+nx+1,:);
utilde((1+(nx+1)):nx+1:1+(ny)*(nx+1),:) = 0;
utilde((2*(nx+1)):nx+1:(ny+1)*(nx+1),:) = 0;
utilde(2+(nx+1)*(ny+1):(ny+2)*(nx+1)-1,:) = 0-utilde(2+(nx+1)*(ny):(ny+1)*(nx+1)-1,:);
% BC for v
vtilde(2:nx,:) = 0;
vtilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:) = 0-vtilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:);
vtilde((2*(nx+2)):nx+2:(ny)*(nx+1),:) = 0-vtilde((2*(nx+2)-1):nx+2:(ny)*(nx+1)-1,:);
vtilde(2+(nx+2)*(ny):(ny+1)*(nx+2)-1,:) = 0;
% Reinforce phi boundary conditions 
phi([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phi(2:nx+1) = phi(2+nx+2:nx+1+nx+2);
phi(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1) = phi(2+(nx+2)*(ny):(nx+2)*(ny+1)-1);
phi(nx+3:nx+2:1+(nx+2)*(ny)) = phi(nx+4:nx+2:2+(nx+2)*(ny));
phi(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phi(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);

%% Predictor Step 2
%% Velocity x
LHSx = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))));
RHSx = dt*(gamma2*ConvectionOperator(utilde,vtilde,1,dx,dy,nx,ny)+zeta1*ConvectionOperator(u,v,1,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))*utilde);
ub = RHSx\(LHSx);
ub = ub';
%% Predictor Step 1 Velocity y
LHSy = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))));
RHSy = dt*(gamma2*ConvectionOperator(utilde,vtilde,2,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy))*vtilde);
vb = RHSy\LHSy;
vb = vb';
%% Poisson's Equation
Int = -(zeta1*(LaplacianFinal(1,1,1,1,nx+2,ny+2,dx,dy))*phi);
phitilde = (Divergence(1,1,1,1,ub,vb,nx,ny,dx,dy)-Int)\(gamma2*LaplacianFinal(1,1,1,1,nx+2,ny+2,dx,dy));
phitilde = phitilde';
%% Apply phi boundary conditions
phitilde([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phitilde(2:nx+1,:) = phitilde(2+nx+2:nx+1+nx+2);
phitilde(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1,:) = phitilde(2+(nx+2)*(ny):(nx+2)*(ny+1)-1,:);
phitilde(nx+3:nx+2:1+(nx+2)*(ny)) = phitilde(nx+4:nx+2:2+(nx+2)*(ny));
phitilde(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phitilde(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);
%% Update Velocities
udoubletilde = ub-dt*(Gradient((zeta1*phi+gamma2*phitilde),1,nx,ny,dx,dy,1,1,1,1));
vdoubletilde = vb-dt*(Gradient((zeta1*phi+gamma2*phitilde),2,nx,ny,dx,dy,1,1,1,1));
%% Update intermediate BCs
% Firstly enforce invalid cell velocities as zero.
udoubletilde([1 nx+1 1+(ny+1)*(nx+1) (nx+1)*(ny+2)],:) = 0;
vdoubletilde([1 nx+2 1+(ny)*(nx+2) (nx+2)*(ny+1)],:) = 0;
% BC for u
udoubletilde(2:nx,:) = 2*U-udoubletilde(2+nx+1:nx+nx+1,:);
udoubletilde((1+(nx+1)):nx+1:1+(ny)*(nx+1),:) = 0;
udoubletilde((2*(nx+1)):nx+1:(ny+1)*(nx+1),:) = 0;
udoubletilde(2+(nx+1)*(ny+1):(ny+2)*(nx+1)-1,:) = 0-udoubletilde(2+(nx+1)*(ny):(ny+1)*(nx+1)-1,:);
% BC for v
vdoubletilde(2:nx,:) = 0;
vdoubletilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:) = 0-vdoubletilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:);
vdoubletilde((2*(nx+2)):nx+2:(ny)*(nx+1),:) = 0-vdoubletilde((2*(nx+2)-1):nx+2:(ny)*(nx+1)-1,:);
vdoubletilde(2+(nx+2)*(ny):(ny+1)*(nx+2)-1,:) = 0;
% Reinforce phi boundary conditions 
phitilde([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phitilde(2:nx+1) = phitilde(2+nx+2:nx+1+nx+2);
phitilde(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1) = phitilde(2+(nx+2)*(ny):(nx+2)*(ny+1)-1);
phitilde(nx+3:nx+2:1+(nx+2)*(ny)) = phitilde(nx+4:nx+2:2+(nx+2)*(ny));
phitilde(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phitilde(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);
%% Predictor Step 3 
%% Velocity x
LHSx = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))));
RHSx = dt*(gamma3*ConvectionOperator(udoubletilde,vdoubletilde,1,dx,dy,nx,ny)+zeta2*ConvectionOperator(utilde,vtilde,1,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+1,ny+2,dx,dy))*udoubletilde);
uc = RHSx\(LHSx);
uc = uc';
%% Predictor Step 3 Velocity y
LHSy = dt*(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy)+0.5*(dx+dy)*ones(size(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy))));
RHSy = dt*(gamma3*ConvectionOperator(udoubletilde,vdoubletilde,2,dx,dy,nx,ny)+(alpha1/Re)*(LaplacianFinal(2,2,2,2,nx+2,ny+1,dx,dy))*vdoubletilde);
vc = RHSy\LHSy;
vc = vc';
%% Poisson's Equation
Int2 = -(zeta2*(LaplacianFinal(1,1,1,1,nx+2,ny+2,dx,dy))*phitilde);
phidoubletilde = (Divergence(1,1,1,1,uc,vc,nx,ny,dx,dy)+Int2)\(gamma3*(LaplacianFinal(1,1,1,1,nx+2,ny+2,dx,dy)));
phidoubletilde = phidoubletilde';
%% Apply phi boundary conditions
phidoubletilde([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phidoubletilde(2:nx+1) = phidoubletilde(2+nx+2:2*(nx+2)-1);
phidoubletilde(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1) = phidoubletilde(2+(nx+2)*(ny):(nx+2)*(ny+1)-1);
phidoubletilde(nx+3:nx+2:1+(nx+2)*(ny)) = phidoubletilde(nx+4:nx+2:2+(nx+2)*(ny));
phidoubletilde(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phidoubletilde(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);
%% Update Velocities
utritilde = uc-dt*(Gradient((zeta2*phitilde+gamma3*phidoubletilde),1,nx,ny,dx,dy,1,1,1,1));
vtritilde = vc-dt*(Gradient((zeta2*phitilde+gamma3*phidoubletilde),2,nx,ny,dx,dy,1,1,1,1));
%% Update intermediate BCs
% Firstly enforce invalid cell velocities as zero.
utritilde([1 nx+1 1+(ny+1)*(nx+1) (nx+1)*(ny+2)],:) = 0;
vtritilde([1 nx+2 1+(ny)*(nx+2) (nx+2)*(ny+1)],:) = 0;
% BC for u
utritilde(2:nx,:) = 2*U-utritilde(2+nx+1:nx+nx+1,:);
utritilde((1+(nx+1)):nx+1:1+(ny)*(nx+1),:) = 0;
utritilde((2*(nx+1)):nx+1:(ny+1)*(nx+1),:) = 0;
utritilde(2+(nx+1)*(ny+1):(ny+2)*(nx+1)-1,:) = 0-utritilde(2+(nx+1)*(ny):(ny+1)*(nx+1)-1,:);
% BC for v
vtritilde(2:nx,:) = 0;
vtritilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:) = 0-vtritilde((1+(nx+2)):nx+2:1+(ny-1)*(nx+2),:);
vtritilde((2*(nx+2)):nx+2:(ny)*(nx+1),:) = 0-vtritilde((2*(nx+2)-1):nx+2:(ny)*(nx+1)-1,:);
vtritilde(2+(nx+2)*(ny):(ny+1)*(nx+2)-1,:) = 0;
% Reinforce phi boundary conditions 
phidoubletilde([1 nx+2 (ny+1)*(nx+2)+1 (ny+2)*(nx+2)]) = 0;
phidoubletilde(2:nx+1) = phidoubletilde(2+nx+2:nx+1+nx+2);
phidoubletilde(2+(nx+2)*(ny+1):(nx+2)*(ny+2)-1) = phidoubletilde(2+(nx+2)*(ny):(nx+2)*(ny+1)-1);
phidoubletilde(nx+3:nx+2:1+(nx+2)*(ny)) = phidoubletilde(nx+4:nx+2:2+(nx+2)*(ny));
phidoubletilde(2*(nx+2):nx+2:(nx+2)*(ny+1)) = phidoubletilde(2*(nx+2)-1:nx+2:(nx+2)*(ny+1)-1);
% Update final velocity
u = utritilde;
v = vtritilde;
CFL = mean((dt/dx).*u,'all');
disp("Timestep")
disp(i)
disp("Courant No",CFL)
disp(CFL)
uvis = reshape(u',[],nx+1);
vvis = reshape(v',[],nx+2);
end





