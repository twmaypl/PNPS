function [Conc]=ModelPnP_PN_ODE15s 

global kB Temp e_unit
global DiffIon ValIon StericG HOT
global AlphaBetaIon AlphaBetaPhi
global tstart
global N CFL c_tau  
global IonBc 
global x TauIon M tend DiffMat_cgl DiffMat_lgl
global Valence bp bm
global kBxT TotNumIon lc_p lc_m S T
global SolOperPoisson 

% % values from PNP_Steric
global phi Current Euler dt  

%%Determine the total number of Ion species
TotNumIon = length(DiffIon);

% Valence of Ion: ValIon=[zp zn]
Valence= ValIon * e_unit;

% Code Begin
% ===================================================================

kBxT = kB*Temp;                     
                    
% Allocate memory for Boundary Condition Variables
  TauIon = zeros([2 1]);
 
% set up Legendre pseusospectral computation elements 
CGL_x = -cos(pi*(0:N)'/N); % Chebyshev-Gauss-Labatto Points
[LGL_x,w,P]=lglnodes(N);  % Legendre-Gauss_Lobatto
LGL_x=flipud(LGL_x);  % LGL_x: Legendre-Gauss-Lobatta grid points
w=flipud(w);          %     w: LGL quadrature weights
M = diag(w); Minv = inv(M);  % M: Mass matrix; Minv: M^{-1}

[dpy, ~] = lepoly(N,LGL_x);    % First derivative of legendre polynomials dy
[dpx, ~] = lepoly(N,CGL_x);    % First derivative of Chebyshev polynomials dx
lc_p = ((1+CGL_x) .* dpx(:))/(2*dpy(N+1)); % Legendre intepolation polynomial using Chebyshev points
lc_m = ((1-CGL_x) .* dpx(:))/(2*dpy( 1 )); % Legendre intepolation polynomial using Chebyshev points

S = LegIntPol(CGL_x,LGL_x,N); % Legendre Interpolation Polynomial
S = S';         % S is the lgl to cgl transform matrix 
T = inv(S);     % T is the cgl to lgl transform matrix


% set up differentiation grid points 
DiffMat_lgl = collocD(LGL_x);     % D: Differentiation matrix
DiffMat_cgl = collocD(CGL_x);   % CGL Differentiation matrix

% setup grid points in physical space and compute transformation metrics
x = CGL_x;


% setup tau parameters based on the imposed BCs
TauIon(1) = 1 / (4 * AlphaBetaIon(1,1) * w( 1 )...
                       + AlphaBetaIon(1,2))/ w( 1 ); % general BC at x=-1
TauIon(2) = 1 / (4 * AlphaBetaIon(2,1) * w(N+1)...
                       + AlphaBetaIon(2,2))/ w(N+1); % general BC at x=-1
 
    if AlphaBetaIon(1,1) == 1
       TauIon(1) = TauIon(1) * c_tau;
    end
    if AlphaBetaIon(2,1) == 1
       TauIon(2) = TauIon(2) * c_tau;
    end


% set up time step and parameters
dx_min=abs(x(2)-x(1));

if Euler == 1
dt = CFL * 0.5/(max(DiffIon)) *dx_min^2;
end


% set up poission solution operator     
% % Robin B.C.
    TauPhi_m = 1/AlphaBetaPhi(1,2); % constant TauPhi_m for Robin
    TauPhi_p = 1/AlphaBetaPhi(2,2); % constant TauPhi_m for Robin
    e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
    eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';

    B_p = TauPhi_p * Minv; 
    bp = S*B_p*eN;
    B_m = TauPhi_m * Minv; 
    bm = S*B_m*e0;
  
    L = - (DiffMat_lgl * DiffMat_lgl ...
        - B_p * (AlphaBetaPhi(2,1) * I_p + AlphaBetaPhi(2,2) ...
          * I_p * DiffMat_lgl) ...
        - B_m * (AlphaBetaPhi(1,1) * I_m - AlphaBetaPhi(1,2) ...
          * I_m * DiffMat_lgl));
    G = M*L ;
%     disp(max(max(G-G')));
    SolOperPoisson= inv(S*L*T);

% set up initial condition
time = tstart; [Conc0] = ConcPhiInit(x,IonBc);
u=Conc0;

tcount = 0;


% Solution figure handle
phi=zeros(N+1,1); Current_x=zeros(size(x));

% VisualProfile(Conc0, phi, Current_x, x, StericG, HOT, dt, tcount, 0);
% pause
fprintf ('N = %i \n', N);
dtime = dt;
while (time  < tend) %
   if time+dtime >= tend
       dtime = tend-time;
   end

fprintf ('time = %d \n', time);

if Euler==0
maxstp=dtime/2;   
tspan=(time:maxstp:time+dtime)';
options=odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',maxstp);%,'Mass',M);
[ttemp,new_utemp]=ode15s(@pnp1d,tspan,u,options);
new_u = new_utemp(end,:)';
else
[dudt] = pnp1d(time,u); 
new_u = u + dtime * dudt;
end
u = new_u;

tcount = tcount+1;
time=time+dtime;

Conc = u(1:N+1);
Current_x = ValIon * Current; 

h=figure(5);
VisualProfile(Conc, phi, Current_x, ...
    x, StericG, HOT, dt, tcount, time);

end

    fclose('all');
    fprintf (' \n');
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function dudt=pnp1d(t,u)
global x N M T
global Valence bp bm DiffIon
global PhiBc HOT dCdxBc lc_p lc_m
global kBxT StericG 
global SolOperPoisson TauIon
global phi Current %phi_ex

Conc = u;


e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';

Qs = Qsource(x,t,Conc);
rhs_b = Qs + Conc*Valence;
rhs_b = rhs_b + bp * PhiBc(2) + bm * PhiBc(1);

phi = SolOperPoisson * rhs_b;

e0 = zeros(N+1,1); e0(1)=1; eN = zeros(N+1,1); eN(N+1)=1;

% Compute dE/dc
dEdc = kBxT * (log(Conc) + 1) + Valence * phi;

Mob = diag(DiffIon/kBxT * Conc);

% add steric effect term
dEdc_ext = StericG * Conc;
    
%add High order term
dEdc_ext(1:N+1,1) = dEdc_ext(1:N+1,1) ...
    + HOT(1,1) * (- chebfft1(chebfft1(Conc,1),1) ...
    + lc_p ./ M(end,end) * (eN' * chebfft1(Conc,1) - dCdxBc(2,1)) ...
    - lc_m ./ M( 1 , 1 ) * (e0' * chebfft1(Conc,1) - dCdxBc(1,1)));

dEdc = dEdc+dEdc_ext;

J = chebfft1(dEdc,1);
Current = - Mob * J;

DivCurrent = - chebfft1 (Current,1);

%Compute divergence current DivCurrent(1:N+1,1:TotNumIon)

g_m = 0; g_p = 0; 

   r_m = (e0' * Mob * e0) * (e0' * dEdc) - g_m;

   r_p = (eN' * Mob * eN) * (eN' * dEdc) - g_p;

    
    % add penalty BC to flux
  
   DivCurrent = DivCurrent - TauIon(1)/M( 1 , 1 ) * (lc_m * r_m) ...
                           - TauIon(2)/M(N+1,N+1) * (lc_p * r_p);


dudt = DivCurrent;

end
