%% Defining Coefficients
alpha1 = 4/15;
alpha2 = 1/15;
alpha3 = 1/6;
gamma1 = 8/15;
gamma2 = 5/12;
gamma3 = 3/4;
zeta1 = -17/60;
zeta2 = -5/12;
Re = 100;
nx = 5;
ny = 5;
dx = 1;
dy = 1;
dt = 1e-03;
%% Modified Runge-Kutta Algorithm
% Predictor Step 1, Velocity x
LHSx = dt*(alpha1/Re)*LaplacianFinal();
RHSx = dt*(gamma1*H(u,1,dx,dy)+(alpha1/Re)*(LaplacianFinal(2,2,1)));
ua = RHSx\LHSx;
% Predictor Step 1 Velocity y
LHSy = dt*(alpha1/Re)*LaplacianM();
RHSy = dt*(gamma1*H(v)+(alpha1/Re)*(Laplace()));
va = RHSy\LHSy;







