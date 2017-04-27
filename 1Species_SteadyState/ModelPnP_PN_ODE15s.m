function [Conc] = ModelPnP_PN_ODE15s 

global kB Temp e_unit dielec 
global DiffIon ValIon StericG HOT
global AlphaBetaIon AreaStr AlphaBetaPhi
global xmin xmax tstart savesteps
global N CFL c_tau dt Euler
global IonBc PhiBc dCdxBc 
global  x TauIon M tend
global Valence bp bm Area
global kBxT TotNumIon DiffMat
global Jacobian SolOperPoisson dxi_dx phi_ex C_ex w

% % values from PNP_Steric
global phi Current

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
[LGL_x,w,P]=lglnodes(N);  
LGL_x=flipud(LGL_x);  % LGL_x: Legendre-Gauss-Lobatta grid points
w=flipud(w);          %     w: LGL quadrature weights
M = diag(w); Minv = inv(M);  % M: Mass matrix; Minv: M^{-1}
                      
% set up differentiation grid   points 
DiffMat=collocD(LGL_x);         % D: Differentiation matrix

% setup grid points in physical space and compute transformation metrics
x = xmin + (xmax-xmin)/(LGL_x(N+1)-LGL_x(1)) * (LGL_x(1:N+1)+1);
dx_dxi = DiffMat*x; dxi_dx = 1./dx_dxi; Jacobian= dx_dxi; 

% compute cross section area
Area     =  eval(AreaStr);

% test case
phi_ex = @(x) - 1/ValIon * (log(x) + 1 + StericG * x ...
    - HOT * (DiffMat * DiffMat * x));
C_ex = @(x) 1+exp(-1) + cos(pi*x);

% setup tau parameters based on the imposed BCs
TauIon(1) = 1 / (4 * AlphaBetaIon(1,1) * w( 1 )...
                       + AlphaBetaIon(1,2)); % general BC at x=-1
TauIon(2) = 1 / (4 * AlphaBetaIon(2,1) * w(N+1)...
                       + AlphaBetaIon(2,2)); % general BC at x=-1

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

    
% % Robin B.C.
    TauPhi_m = 1/AlphaBetaPhi(1,2); % constant TauPhi_m for Robin
    TauPhi_p = 1/AlphaBetaPhi(2,2); % constant TauPhi_m for Robin
    e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
    eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';
    A = diag(Jacobian .* dxi_dx .* dielec .* Area .* dxi_dx); 
    
    B_p = TauPhi_p * Minv;% * A;
    bp = B_p*eN;
    B_m = TauPhi_m * Minv;% *  A;
    bm = B_m*e0;
    L = - (DiffMat * A * DiffMat ...
        - B_p * (AlphaBetaPhi(2,1) * I_p + AlphaBetaPhi(2,2) ...
          * I_p * diag(dxi_dx) * DiffMat) ...
        - B_m * (AlphaBetaPhi(1,1) * I_m - AlphaBetaPhi(1,2) ...
          * I_m * diag(dxi_dx) * DiffMat));
    G = M*L ;
% %     disp(max(max(G-G')));
% %     disp(eig(G));
% %     pause 
SolOperPoisson= inv(G);

% set up initial condition
time = tstart; [Conc0] = ConcPhiInit(x,IonBc);
u=Conc0;

tcount = 0;



% Solution figure handle
phi=zeros(N+1,1); Current_x=zeros(size(x)); 
% VisualProfile(Conc0(:,1), phi, Current_x, x, StericG, HOT, dt, tcount, 0);
% pause

fprintf ('N = %i \n', N);
dtime = dt;
while (time  < tend) %
   if time+dtime >= tend
       dtime = tend-time;
   end

fprintf('time = %d \n', time);

if Euler == 0
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

Conc = u;


% % % % % % % % %Electrio current
Current_x = ValIon * Current;  

h=figure(5);
VisualProfile(Conc, phi, Current_x, x, StericG, HOT, dt, tcount, time);

end

fprintf ('\n');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function dudt=pnp1d(t,Conc)
global x N  M phi_ex
global Valence bp bm Area DiffIon AlphaBetaIon
global PhiBc HOT dCdxBc 
global kBxT StericG DiffMat 
global Jacobian SolOperPoisson dxi_dx TauIon
% % values from PNP_Steric
global phi Current 

e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';

Qs = Qsource(x,t, Conc);
rhs_b = Qs + Conc*Valence;
% rhs_b = Jacobian .* rhs_b .* Area + bp * PhiBc(2) + bm * PhiBc(1);
rhs_b = rhs_b + bp * (eN' * phi_ex(Conc)) ...
              + bm * (e0' * phi_ex(Conc));
phi = SolOperPoisson * M * rhs_b; 

A = Jacobian .* Area .* dxi_dx;
e0 = zeros(N+1,1); e0(1)=1; eN = zeros(N+1,1); eN(N+1)=1;

% Compute dE/dc
dEdc = kBxT * (log(Conc) + 1) ...
    + Valence * phi;
Mob = DiffIon/kBxT * Conc;

% add steric effect term
dEdc_ext = StericG * Conc;
    
%add High order term

dEdc_ext(1:N+1,1) = dEdc_ext(1:N+1,1) ...
    + HOT(1,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx)...
                * DiffMat * Conc - dCdxBc(2,1)) ...
                - e0 / M( 1 , 1 ) * dxi_dx( 1 ) * (e0' * diag(dxi_dx)...
                * DiffMat * Conc - dCdxBc(1,1)));

% compute current J = -  D/kBT * Conc * D_x(dEdc) )

dEdc = dEdc+dEdc_ext;

Current = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);

%Compute divergence current DivCurrent(1:N+1,1:TotNumIon)

DivCurrent = - diag(1./(Jacobian.*Area)) ...
    * (DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current));


g_m = 0; g_p = 0;
 
    r_m = (AlphaBetaIon(1,1) ...
        * Mob( 1 ) * dEdc( 1 ) ...
        - AlphaBetaIon(1,2) * (-Current( 1 ))  - g_m );

    r_p = (AlphaBetaIon(2,1) ...
        * Mob(N+1) * dEdc(N+1)  ...
        + AlphaBetaIon(2,2) * (-Current(N+1)) - g_p );
    
% add penalty BC to flux 
    DivCurrent( 1 ) = DivCurrent( 1 ) - TauIon(1) / M( 1 , 1 ) * r_m;
    DivCurrent(N+1) = DivCurrent(N+1) - TauIon(2) / M(N+1,N+1) * r_p;

dudt = DivCurrent;

end
