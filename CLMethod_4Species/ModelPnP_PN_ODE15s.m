% The solver of Poisson-Nerst-Planck equation couple wiht steric effect
% update log:
%        date: 2015/10/08
%              Output resuts in Tecplot360 ASCII form
%        date: 2015/10/05
%              Compute energy and total current
%        date: 2015/09/28
%              add high order term
%              plot I_ex, I_pnp, I_totle
%        date: 2015/08/26
%              add to 4 species
%              add one more ion species to 3 species
%        date: 2015/08/15
%              bug fixed : dEdc + gij * Cj 
% Equations:
%     - 1/A d/dx (A e d(phi)/dx) = rho
%     dE/dc = A (KbT(log u +1) + z e phi + sum{g_{ij} u}
%     Mob = Diffusion u / KbT
%     Flux = - Mob d/dx(1/A dE/dc)
%     A du/dt = -d/dx(A Flux)
% Using Legendre-Gauss-Lobatto points and Lagrange interpolation
% polynomials to semi-descritise the equations
%  
% Poisson solver 
%
%     L phi = -[DAD - B_{+}(a_{+}I_{+} + b_{+}I_{+}D) 
%                   - B_{-}(a_{-}I_{-} + b_{-}I_{-}D)] phi = rho
%     Dirichlet : a_{+-} = 1, b_{+-} = 0 
%            B_{+-} = +-M^{-1} (D^{T} _+ tau_{+-} I) A
%     Then phi = L^{-1} rho
%
% Nerst-Planck with Steric effect
%
%     du/dt = - 1/A [D A (-Mob D dE/dc) 
%             - 1/A_{-} tau_{-} A_{-} (Mob dE/dc_{-} - g_{-})
%             - 1/A_{+} tau_{+} A_{+} (Mob dE/dc_{+} - g_{+})]
%     Solved by ODE15s
%   

function ModelPnP_PN_ODE15s(dt) 

global kB Temp e_unit dielec %td epslion strong
global DiffIon ValIon StericG HOT
global AlphaBetaIon AreaStr AlphaBetaPhi
global xmin xmax tstart savesteps
global N CFL delta c_tau IonBc PhiBc ddxAdCdxBc dCdxBc
global err x TauIon M tend
global Valence bp bm Area
global kBxT TotNumIon DiffMat
global Jacobian SolOperPoisson dxi_dx

%% values from PNP_Steric
global phi I_PNP I_ext Current dEdt En_total dCdT
                                
%% Determine the total number of Ion species
TotNumIon = length(DiffIon);

% Valence of Ion: ValIon=[zp zn]
Valence= ValIon * e_unit;

kBxT = kB*Temp;                     
% c_tau = 2.0;     % scaling parameter for adjusting 
                 % penalty parameters tau
                    
% Allocate memory for Boundary Condition Variables
  TauIon = zeros([2 TotNumIon]);
 
% Code Begin
% ===================================================================

% set up Legendre pseusospectral computation elements 
[LGL_x,w,P]=lglnodes(N);  
LGL_x=flipud(LGL_x);  % LGL_x: Legendre-Gauss-Lobatta grid points
w=flipud(w);          %     w: LGL quadrature weights
M = diag(w); Minv = inv(M);
                      
% set up differentiation grid   points 
DiffMat=collocD(LGL_x);         % D: Differentiation matrix

% setup grid points in physical space and compute transformation metrics
x = xmin + (xmax-xmin)/(LGL_x(N+1)-LGL_x(1)) * (LGL_x(1:N+1)+1);
dx_dxi = DiffMat*x; dxi_dx = 1./dx_dxi; Jacobian= dx_dxi; 

% compute cross section area
Area     =  eval(AreaStr);

% setup tau parameters based on the imposed BCs

TauIon(1:2,:) = 1 / (4 * AlphaBetaIon(1:2,1,:) * w( 1 )...
                       + AlphaBetaIon(1:2,2,:))/ w( 1 ); % general BC at x=-1
for IonNum = 1:TotNumIon
 
    if AlphaBetaIon(1,1,IonNum) == 1
       TauIon(1,IonNum) = TauIon(1,IonNum) * c_tau;
    end
    if AlphaBetaIon(2,1,IonNum) == 1
       TauIon(2,IonNum) = TauIon(2,IonNum) * c_tau;
    end

end

% set up time step and parameters
dx_min=abs(x(2)-x(1));
stand_dt = CFL * 0.5/(max(DiffIon)) *dx_min^delta;
Energy_OLD = 0.0;
%%%%%nc = ceil((tend-tstart)/dt);

% set up poission solution operator 
% % % if strong,
% % %     disp('strongly solution');
% % %     L = - diag((1./Jacobian)) * DiffMat ...
% % %         * (diag(Area .* dielec .* Jacobian .* dxi_dx .* dxi_dx) * DiffMat);
% % %     L( 1 ,1:N+1)=0.0; L( 1 , 1 ) = 1;
% % %     L(N+1,1:N+1)=0.0; L(N+1,N+1) = 1;
% % % else
% % %     disp('weakly solution');
% % Dirichlet B.C.
% %     TauPhi_m = c_tau * 1/w(1);   % constant TauPhi_m for Dirichlet
% %     TauPhi_p = c_tau * 1/w(N+1); % constant TauPhi_p for Dirichlet
% %     e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
% %     eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';
% %     A = diag(Jacobian .* dxi_dx .* dielec .* Area .* dxi_dx); 
% % 
% %     B_p  = - Minv * (transpose(DiffMat) - TauPhi_p * eye(size(DiffMat))) * A;
% %     bp = B_p*eN;
% %     B_m = + Minv * (transpose(DiffMat) + TauPhi_m * eye(size(DiffMat))) * A;
% %     bm = B_m*e0;
% %     L = - (DiffMat * A * DiffMat - B_p * I_p - B_m * I_m);
% %     G = M*L;

% % Robin B.C.
    TauPhi_m = 1/AlphaBetaPhi(1,2); % constant TauPhi_m for Robin
    TauPhi_p = 1/AlphaBetaPhi(2,2); % constant TauPhi_m for Robin
    e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
    eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';
    A = diag(Jacobian .* dxi_dx .* dielec .* Area .* dxi_dx); 
    
    B_p = TauPhi_p * M \ A;
    bp = B_p*eN;
    B_m = TauPhi_m * M \ A;
    bm = B_m*e0;
    L = - (DiffMat * A * DiffMat ...
        - B_p * (AlphaBetaPhi(2,1) * I_p + AlphaBetaPhi(2,2) ...
          * I_p * diag(dxi_dx) * DiffMat) ...
        - B_m * (AlphaBetaPhi(1,1) * I_m - AlphaBetaPhi(1,2) ...
          * I_m * diag(dxi_dx) * DiffMat));
    G = M*L ;

% % % end
SolOperPoisson= inv(G); 

% set up initial condition
time = tstart; [Conc0] = ConcPhiInit(x,IonBc);
u=[Conc0(:,1);Conc0(:,2);Conc0(:,3);Conc0(:,4)];

tcount = 0;
% utemp = u;
Energy = 0;
Energy_34 = 0;
Current_t = 0;
Current_OLD = 0;

% Solution figure handle
phi=zeros(N+1,1); Current_x=zeros(size(x)); I = zeros(N+1,1);I_ext=I;I_PNP=I;
 Energy_t=0; Energy_err=0; Current_total=0;
VisualProfile(Conc0(:,1),Conc0(:,2),Conc0(:,3),Conc0(:,4),phi,Current_x, ...
    I_ext, I_PNP, Current_t, Energy, x, StericG, HOT, dt, tcount, time,1, ...
    Energy_t, Energy_err, Current_total);
%pause

A = Jacobian .* Area .* dxi_dx;

% Set up output folder and file name
StrSteric = [num2str(StericG(1,1)) ' ' num2str(StericG(1,2)) ' ' ...
        num2str(StericG(3,3)) ' '  num2str(StericG(3,4))];
StrHOT = [num2str(HOT(1,1)) ' ' num2str(HOT(1,2)) ' ' ...
        num2str(HOT(3,3)) ' ' num2str(HOT(3,4))];
txtfolder = ['./Results/N' num2str(N) '/tend=' ...
        num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT ...
        '/Conc1=' num2str(Conc0(1,1)) ' Conc3=' num2str(Conc0(1,3)) '/'];
    
if ~exist(txtfolder, 'dir')
   mkdir(txtfolder);
%    warning (['creating new folder' txtfolder]);
end

filename = ['C1-=' num2str(Conc0(1,1)) ' C1+=' num2str(Conc0(end,1)) ...
    ' C3-=' num2str(Conc0(1,2)) ' C3+=' num2str(Conc0(end,2)) ...
    ' Phi-=' num2str(PhiBc(1)) ' Phi+=' num2str(PhiBc(2)) ' Qs=' num2str(1e+0)];
outputInfo = [txtfolder '/' filename '.txt'];
InfotxtID = fopen(outputInfo, 'w');
outputName = [txtfolder '/' filename '.dat'];
fileID = fopen(outputName, 'w');
% ENERGY = repmat(Energy_t, length(x),1);
fprintf(fileID, 'TITLE = "2 species of PNP-Steric with high order term" \n');
fprintf(fileID, ['VARIABLES = "X" "Conc1" "Conc2" "Conc3" "Conc4" "phi"'...
    ' "Current_x" "I_ext" "I_PNP" "CURRENT" "ENERGY" \n']);


while (time  < tend) %
   if time+dt >= tend
       dt = tend-time;
   end
  
maxstp=dt/2;
fprintf ('time = %d \n', time);
fprintf (InfotxtID, 'time = %d \n', time);

tspan=(time:maxstp:time+dt)';
options=odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',maxstp);%,'Mass',M);


[ttemp,new_utemp]=ode15s(@pnp1d,tspan,u,options);
new_u=new_utemp(end,:)';
dCdt_err = max(max(abs(new_u-u)));
fprintf('relative err of Conc = %d \n', dCdt_err);
fprintf(InfotxtID, 'relative err of Conc = %d \n', dCdt_err);
u = new_u;
% u=utemp;

tcount=tcount+1;
time=time+dt;

Conc1 = u(1:N+1); Conc2 = u(N+2:2*N+2);
Conc3 = u(2*N+3:3*N+3); Conc4 = u(3*N+4:4*N+4);
Conc(1:N+1,1) = Conc1(:,end);
Conc(1:N+1,2) = Conc2(:,end);
Conc(1:N+1,3) = Conc3(:,end);
Conc(1:N+1,4) = Conc4(:,end);
% testC1 = [Conc1<0];
% testC2 = [Conc2<0];
% testC3 = [Conc3<0];
% testC4 = [Conc4<0];
% disp([nnz(testC1), nnz(testC2)]);

fprintf('dEdt = %d \n', dEdt);
fprintf(InfotxtID, 'dEdt = %d \n', dEdt);
%% check the steady-state %%+++++++++++++++++++++++++++++++++++++++++++++++
dConcdt = max(max(abs(Current_OLD - Current)));
fprintf('relative error of Current = %d \n', dConcdt);
fprintf(InfotxtID, 'relative error of Current = %d \n', dConcdt);
fprintf(InfotxtID, 'dCdt = %d \n', max(max(abs(dCdT))));
fprintf('dCdt = %d \n', max(max(abs(dCdT))));

%% Electrio current %%+++++++++++++++++++++++++++++++++++++++++++++++++++++
Current_x = sum(repmat(ValIon,N+1,1) .* Current,2);
Current_total = sum(w.*Current_x);
Current_t = [Current_t; Current_total];
fprintf('Current_total = %d \n', Current_total);
fprintf(InfotxtID, 'Current_total = %d \n', Current_total);

A = Jacobian .* Area .* dxi_dx;
Energy_t = sum(w .* A .* En_total) ;
Energy_err = Energy_t - Energy_OLD;
Energy = [Energy; Energy_t];
Energy_OLD = Energy_t;
fprintf('Energy_err = %d \n', Energy_err);
fprintf('Energy_t = %d \n', Energy_t); 
fprintf(InfotxtID,'Energy_err = %d \n', Energy_err);
fprintf(InfotxtID,'Energy_t = %d \n', Energy_t); 

fprintf(' \n');
fprintf(InfotxtID, ' \n');
%%% =======================================================================

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Qs = Qsource(x,time);
if mod(tcount,savesteps) == 0
   matfolder = ['./Results/N' num2str(N) '/tend=' ...
        num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT ...
        '/Conc1=' num2str(Conc0(1,1)) ' Conc3=' num2str(Conc0(1,3)) '/matFile/' ...
        '/Q=' num2str(Qs(1)) ' z1=' num2str(ValIon(1)) ' z3=' num2str(ValIon(3)) ...
        ' phi_l=' num2str(PhiBc(1)) ' phi_r=' num2str(PhiBc(2)) ]; 
   
   if ~exist(matfolder, 'dir')
       mkdir(matfolder);
   %    print("creating new folder", matfolder);
   end
   
   save([matfolder, '/', num2str(tcount/savesteps),'.mat'], ...
       'x', 'Conc','phi', 'Current_x', 'I_ext', 'I_PNP', ...
       'Energy', 'Current_t', 'time', 'Energy_34', 'Energy_err');
end   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Output results to .txt file
    CURRENT = repmat(Current_total, length(x),1);
    ENERGY = repmat(Energy_t, length(x),1);
    fprintf(fileID, ['ZONE T="' num2str(time) '", I=' num2str(N+1) '\n']);
    A = [x Conc(:,1) Conc(:,2) Conc(:,3) Conc(:,4) phi Current_x I_ext I_PNP CURRENT ENERGY];
    fprintf(fileID, '%d %d %d %d %d %d %d %d %d %d %d \n', A');
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% if mod(tcount,savesteps) == 0
% Plotout = ['./Results/N' num2str(N) '/tend=' ...
%         num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT ...
%         '/Conc1=' num2str(Conc0(1,1)) ' Conc3=' num2str(Conc0(1,3)) '/PlotsProfile/' ...
%         '/Q=' num2str(Qs(1)) ' z1=' num2str(ValIon(1)) ' z3=' num2str(ValIon(3)) ...
%         ' phi_l=' num2str(PhiBc(1)) ' phi_r=' num2str(PhiBc(2))];
% if ~exist(Plotout, 'dir')
%         mkdir(Plotout);
% end
% h=figure(1);
% VisualProfile(Conc(:,1),Conc(:,2),Conc(:,3),Conc(:,4),phi,Current_x, ...
%     I_ext, I_PNP, Current_t, Energy, x, StericG, HOT, dt, tcount, time, 1, ...
%     Energy_t, Energy_err, Current_total);
% 
% saveas(h, [Plotout '/' num2str(tcount/savesteps) '.eps'], 'psc2');
% 
% end
end

fclose('all');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function dudt=pnp1d(t,u)
global x N err M
global Valence bp bm IonBc PhiBc Area DiffIon AlphaBetaIon
global kBxT TotNumIon StericG DiffMat ValIon 
global  ddxAdCdxBc HOT ChemBc dCdxBc e_unit
global Jacobian SolOperPoisson dxi_dx TauIon
%% values from PNP_Steric
global phi I_PNP I_ext Current dEdt En_total dCdT

%%
Conc = zeros([N+1 TotNumIon]);
Conc(1:N+1,1) = u(1:N+1);
Conc(1:N+1,2) = u(N+2:2*N+2);
Conc(1:N+1,3) = u(2*N+3:3*N+3);
Conc(1:N+1,4) = u(3*N+4:4*N+4);

rhs_b = zeros(size(x)); phi = zeros(size(x)); Qs = zeros(size(x));
dEdc = zeros(size(Conc)); Mob  = zeros(size(Conc));  
DivCurrent = zeros(size(Conc));

Qs = Qsource(x,t);
% rhs_b = Qs + Conc(:,1)*Valence(1) + Conc(:,2)*Valence(2);
rhs_b = Qs + sum(Conc .* repmat(Valence,length(Conc(:,1)),1),2); 
rhs_b = Jacobian .* rhs_b .* Area + bp * PhiBc(2) + bm * PhiBc(1);
% phi = SolOperPoisson * rhs_b;
phi = SolOperPoisson * M * rhs_b;

A = Jacobian .* Area .* dxi_dx;
e0 = zeros(N+1,1); e0(1)=1; eN = zeros(N+1,1); eN(N+1)=1;
    
% Compute dE/dc
dEdc = kBxT * (log(Conc) + 1) ...
    + repmat(Valence,N+1,1) .* repmat(phi,1,TotNumIon);

Mob = repmat(DiffIon,N+1,1)/kBxT .* Conc;
Current_PNP = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);
% % % %Electrio current
I_PNP = sum(repmat(ValIon,N+1,1) .* Current_PNP,2);

% add steric effect term
dEdc_ext(1:N+1,1) = sum(repmat(StericG(1,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
dEdc_ext(1:N+1,2) = sum(repmat(StericG(2,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
dEdc_ext(1:N+1,3) = sum(repmat(StericG(3,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
dEdc_ext(1:N+1,4) = sum(repmat(StericG(4,1:4),N+1,1) .* Conc (1:N+1,1:4),2);

% % % % % Couculate dCdxBC
% % % % % dCdx = dEdc;
% % % % % dCdxHOT(1:N+1,1) = ...
% % % % %     + HOT(1,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1)) ...
% % % % %     + HOT(1,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2)) ...
% % % % %     + HOT(1,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3)) ...
% % % % %     + HOT(1,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4));
% % % % % 
% % % % % dCdxHOT(1:N+1,2) = ...
% % % % %     + HOT(2,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1)) ...
% % % % %     + HOT(2,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2)) ...
% % % % %     + HOT(2,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3)) ...
% % % % %     + HOT(2,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4));
% % % % % 
% % % % % dCdxHOT(1:N+1,3) = ...
% % % % %     + HOT(3,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1)) ...
% % % % %     + HOT(3,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2)) ...
% % % % %     + HOT(3,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3)) ...
% % % % %     + HOT(3,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4));
% % % % % 
% % % % % dCdxHOT(1:N+1,4) = ...
% % % % %     + HOT(4,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1)) ...
% % % % %     + HOT(4,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2)) ...
% % % % %     + HOT(4,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3)) ...
% % % % %     + HOT(4,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4));
% % % % % % dCdxHOT = [0 0 0 0; dCdxHOT(2:N,1:4); 0 0 0 0];
% % % % % dCdx = dCdx + dCdxHOT;
% % % % % 
% % % % % Mob = repmat(DiffIon,N+1,1)/kBxT .* Conc;
% % % % % dCdxBc(1,1:4) = -Mob( 1 ,1:end) .* ((e0' * DiffMat*DiffMat*dCdx)./(e0' *DiffMat*dCdx));
% % % % % dCdxBc(2,1:4) = -Mob(N+1,1:end) .* ((eN' * DiffMat*DiffMat*dCdx)./(eN' *DiffMat*dCdx));
% % disp(dCdxBc)
% % % % pause
%add High order term

dEdc_ext(1:N+1,1) = dEdc_ext(1:N+1,1) ...
    + HOT(1,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(2,1)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(1,1))) ...
    + HOT(1,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(2,2)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(1,2)))...
    + HOT(1,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(2,3)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(1,3)))...
    + HOT(1,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(2,4)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(1,4))); 
            
dEdc_ext(1:N+1,2) = dEdc_ext(1:N+1,2) ...
    + HOT(2,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(2,1)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(1,1))) ...
    + HOT(2,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(2,2)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(1,2))) ...
    + HOT(2,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(2,3)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(1,3))) ...
    + HOT(2,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(2,4)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(1,4)));
            
dEdc_ext(1:N+1,3) = dEdc_ext(1:N+1,3) ...
    + HOT(3,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(2,1)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(1,1))) ...
    + HOT(3,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(2,2)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(1,2))) ...
    + HOT(3,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(2,3)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(1,3))) ...
    + HOT(3,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(2,4)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(1,4)));
            
dEdc_ext(1:N+1,4) = dEdc_ext(1:N+1,4) ...
    + HOT(4,1) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,1) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(2,1)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,1) - dCdxBc(1,1))) ...
    + HOT(4,2) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,2) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(2,2)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,2) - dCdxBc(1,2))) ...
    + HOT(4,3) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,3) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(2,3)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,3) - dCdxBc(1,3))) ...
    + HOT(4,4) * (- diag(1./A) * DiffMat * diag(A) * DiffMat * Conc(1:N+1,4) ...
                + eN / M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(2,4)) ...
                - e0 / M(1,1) * dxi_dx(1) * (e0' * diag(dxi_dx) * DiffMat * Conc(1:N+1,4) - dCdxBc(1,4)));

% compute current J = -  D/kBT * Conc * D_x(dEdc) )

Current_ext = - Mob .* (diag(dxi_dx)*DiffMat*dEdc_ext);
I_ext = sum(repmat(ValIon,N+1,1) .* Current_ext,2);

dEdc = dEdc+dEdc_ext;
Current = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);

% % % %Electrio current
% % % I_total = repmat(ValIon(1),N+1,1) .* Current(:,1) + repmat(ValIon(2),N+1,1) .* Current(:,2);  

%Compute divergence current DivCurrent(1:N+1,1:TotNumIon)
DivCurrent = - diag(1./(Jacobian.*Area)) ...
    * (DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current));

for IonNum=1:TotNumIon
    
    % modify flux function with BCs  
    % evaluate exact u for imposing Dirichlet bc for u
    if ChemBc == 0 
    if IonNum == 1
        g_m_Cpp = HOT(1,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
               + HOT(1,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2)) ...
               + HOT(1,3) * (- 1/(e0'*A) * ddxAdCdxBc(1,3)) ...
               + HOT(1,4) * (- 1/(e0'*A) * ddxAdCdxBc(1,4));
           
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,1) * (kBxT * ...
           (log(IonBc(1,1)) + 1) + Valence(1) * PhiBc(1) + ...
           sum(StericG(1,1:4) .* IonBc(1,1:4))+g_m_Cpp);
       
       g_p_Cpp = HOT(1,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
               + HOT(1,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2)) ...
               + HOT(1,3) * (- 1/(eN'*A) * ddxAdCdxBc(2,3)) ...
               + HOT(1,4) * (- 1/(eN'*A) * ddxAdCdxBc(2,4)); 
           
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,1) * (kBxT * ...
           (log(IonBc(2,1)) + 1) + Valence(1) * PhiBc(2) + ...
           sum(StericG(1,1:4) .* IonBc(2,1:4)) + g_p_Cpp);            
    end 

    if IonNum == 2
        g_m_Cpp = HOT(2,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
               + HOT(2,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2)) ...
               + HOT(2,3) * (- 1/(e0'*A) * ddxAdCdxBc(1,3)) ...
               + HOT(2,4) * (- 1/(e0'*A) * ddxAdCdxBc(1,4));
           
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,2) * (kBxT * ...
           (log(IonBc(1,2)) + 1) + Valence(2) * PhiBc(1) + ...
           sum(StericG(2,1:4) .* IonBc(1,1:4)) + g_m_Cpp);
       
       g_p_Cpp = HOT(2,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
               + HOT(2,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2)) ...
               + HOT(2,3) * (- 1/(eN'*A) * ddxAdCdxBc(2,3)) ...
               + HOT(2,4) * (- 1/(eN'*A) * ddxAdCdxBc(2,4));
       
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,2) * (kBxT * ...
           (log(IonBc(2,2)) + 1) + Valence(2) * PhiBc(2) + ...
           sum(StericG(2,1:4) .* IonBc(2,1:4)) + g_p_Cpp);
    end
    
    if IonNum == 3
        g_m_Cpp = HOT(3,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
               + HOT(3,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2)) ...
               + HOT(3,3) * (- 1/(e0'*A) * ddxAdCdxBc(1,3)) ...
               + HOT(3,4) * (- 1/(e0'*A) * ddxAdCdxBc(1,4));
           
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,3) * (kBxT * ...
           (log(IonBc(1,3)) + 1) + Valence(3) * PhiBc(1) + ...
           sum(StericG(3,1:4) .* IonBc(1,1:4)) + g_m_Cpp);
           
       g_p_Cpp = HOT(3,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
               + HOT(3,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2)) ...
               + HOT(3,3) * (- 1/(eN'*A) * ddxAdCdxBc(2,3)) ...
               + HOT(3,4) * (- 1/(eN'*A) * ddxAdCdxBc(2,4));
           
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,3) * (kBxT * ...
           (log(IonBc(2,3)) + 1) + Valence(3) * PhiBc(2) + ...
           sum(StericG(3,1:4) .* IonBc(2,1:4)) + g_p_Cpp);            
    end 
    
    if IonNum == 4
        g_m_Cpp = HOT(4,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
               + HOT(4,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2)) ...
               + HOT(4,3) * (- 1/(e0'*A) * ddxAdCdxBc(1,3)) ...
               + HOT(4,4) * (- 1/(e0'*A) * ddxAdCdxBc(1,4));
           
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,4) * (kBxT * ...
           (log(IonBc(1,4)) + 1) + Valence(4) * PhiBc(1) + ...
           sum(StericG(4,1:4) .* IonBc(1,1:4)) + g_m_Cpp);
           
       g_p_Cpp = HOT(4,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
               + HOT(4,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2)) ...
               + HOT(4,3) * (- 1/(eN'*A) * ddxAdCdxBc(2,3)) ...
               + HOT(4,4) * (- 1/(eN'*A) * ddxAdCdxBc(2,4));
           
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,4) * (kBxT * ...
           (log(IonBc(2,4)) + 1) + Valence(4) * PhiBc(2) + ...
           sum(StericG(4,1:4) .* IonBc(2,1:4)) + g_p_Cpp);            
    end 

    else
        g_m = 0; g_p = 0;
    end
    
    r_m = Area( 1 ) * (AlphaBetaIon(1,1,IonNum) ...
        * (Mob( 1 ,IonNum) * dEdc( 1 ,IonNum) - g_m) ...
        - AlphaBetaIon(1,2,IonNum) * (-Current( 1 ,IonNum) - 0) );

    r_p = Area(N+1) * (AlphaBetaIon(2,1,IonNum) ...
        * (Mob(N+1,IonNum) * dEdc(N+1,IonNum) - g_p) ...
        + AlphaBetaIon(2,2,IonNum) * (-Current(N+1,IonNum) - 0) );
    
    % add penalty BC to flux
    DivCurrent( 1 ,IonNum) = DivCurrent( 1 ,IonNum) - TauIon(1,IonNum) * r_m;
    DivCurrent(N+1,IonNum) = DivCurrent(N+1,IonNum) - TauIon(2,IonNum) * r_p;
    
end

% dEdt = sum( w .* sum(((- DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current)) .* dEdc),2));
dEdt = sum(M * sum(((- DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current)) .* dEdc),2));
%%%% compute energy ====================================================
E_NP = kBxT * sum(Conc.*log(Conc),2); 
E_Poisson = 1/2 * e_unit * (Qs + sum(repmat(ValIon,N+1,1).*Conc,2)) .* phi;
% % % E_Steric = 1/2 * (StericG(1,1) * Conc(:,1) .* Conc(:,1) + ...
% % %                 StericG(1,2) * Conc(:,1) .* Conc(:,2) + ...
% % %                 StericG(1,3) * Conc(:,1) .* Conc(:,3) + ...
% % %                 StericG(1,4) * Conc(:,1) .* Conc(:,4) + ...
% % %                 StericG(2,1) * Conc(:,2) .* Conc(:,1) + ...
% % %                 StericG(2,2) * Conc(:,2) .* Conc(:,2) + ...
% % %                 StericG(2,3) * Conc(:,2) .* Conc(:,3) + ...
% % %                 StericG(2,4) * Conc(:,2) .* Conc(:,4) + ...
% % %                 StericG(3,1) * Conc(:,3) .* Conc(:,1) + ...
% % %                 StericG(3,2) * Conc(:,3) .* Conc(:,2) + ...
% % %                 StericG(3,3) * Conc(:,3) .* Conc(:,3) + ...
% % %                 StericG(3,4) * Conc(:,3) .* Conc(:,4) + ...
% % %                 StericG(4,1) * Conc(:,4) .* Conc(:,1) + ...
% % %                 StericG(4,2) * Conc(:,4) .* Conc(:,2) + ...
% % %                 StericG(4,3) * Conc(:,4) .* Conc(:,3) + ...
% % %                 StericG(4,4) * Conc(:,4) .* Conc(:,4));
% % %            
% % % E_HOT = 1/2 * (HOT(1,1) * (DiffMat * Conc(:,1)) .* (DiffMat * Conc(:,1)) + ...
% % %                HOT(1,2) * (DiffMat * Conc(:,1)) .* (DiffMat * Conc(:,2)) + ...
% % %                HOT(1,3) * (DiffMat * Conc(:,1)) .* (DiffMat * Conc(:,3)) + ...
% % %                HOT(1,4) * (DiffMat * Conc(:,1)) .* (DiffMat * Conc(:,4)) + ...
% % %                HOT(2,1) * (DiffMat * Conc(:,2)) .* (DiffMat * Conc(:,1)) + ...
% % %                HOT(2,2) * (DiffMat * Conc(:,2)) .* (DiffMat * Conc(:,2)) + ...
% % %                HOT(2,3) * (DiffMat * Conc(:,2)) .* (DiffMat * Conc(:,3)) + ...
% % %                HOT(2,4) * (DiffMat * Conc(:,2)) .* (DiffMat * Conc(:,4)) + ...
% % %                HOT(3,1) * (DiffMat * Conc(:,3)) .* (DiffMat * Conc(:,1)) + ...
% % %                HOT(3,2) * (DiffMat * Conc(:,3)) .* (DiffMat * Conc(:,2)) + ...
% % %                HOT(3,3) * (DiffMat * Conc(:,3)) .* (DiffMat * Conc(:,3)) + ...
% % %                HOT(3,4) * (DiffMat * Conc(:,3)) .* (DiffMat * Conc(:,4)) + ...
% % %                HOT(4,1) * (DiffMat * Conc(:,4)) .* (DiffMat * Conc(:,1)) + ...
% % %                HOT(4,2) * (DiffMat * Conc(:,4)) .* (DiffMat * Conc(:,2)) + ...
% % %                HOT(4,3) * (DiffMat * Conc(:,4)) .* (DiffMat * Conc(:,3)) + ...
% % %                HOT(4,4) * (DiffMat * Conc(:,4)) .* (DiffMat * Conc(:,4)));
En_total = E_NP + E_Poisson;% + E_Steric + E_HOT;
%% 

dudt = [DivCurrent(:,1);DivCurrent(:,2);DivCurrent(:,3);DivCurrent(:,4)];
dCdT = dudt;
% %  %% check the steady-state %%+++++++++++++++++++++++++++++++++++++++++++++++
% % 
% % fprintf('dCdt = %f \n', max(max(abs(dudt))));
% % if max(max(abs(dudt))) <= 1e-8
% %     fprintf('steady-state occur at t = %f \n', t);
% % end
    
% err = max(abs(dudt));
% format short e
%   disp([t err])
% pause;
end

