function ModelPnP_PN_ODE15s(dt) 
% function [x Conc phi] = ModelPnP_PN_ODE15s %(Bc)
% Update on 2015/10/14
% Change log 
%        Change output folder and format of plots
% Update on 2015/10/08
% Change log
%        Change output file .txt format to tecplot360 format
% Update on 2015/10/05
% Change log
%        Compute energy and total current
% Update on 2015/09/03
% Change log
%          Output results every savesteps steps


global kB Temp e_unit dielec %td epslion strong
global DiffIon ValIon StericG HOT
global AlphaBetaIon AreaStr AlphaBetaPhi
global xmin xmax tstart savesteps
global N CFL delta c_tau  
global IonBc PhiBc ddxAdCdxBc dCdxBc DiffMat err 
global x TauIon M tend
global Valence bp bm Area
global kBxT TotNumIon lc_p lc_m S T
global Jacobian SolOperPoisson dxi_dx

% % values from PNP_Steric
global phi I_PNP I_ext Current dEdt En_total dCdT

%%Determine the total number of Ion species
TotNumIon = length(DiffIon);

% Valence of Ion: ValIon=[zp zn]
Valence= ValIon * e_unit;

% Code Begin
% ===================================================================

kBxT = kB*Temp;                     
% c_tau = 2.0;     % scaling parameter for adjusting 
                 % penalty parameters tau
                    
% Allocate memory for Boundary Condition Variables
  TauIon = zeros([2 TotNumIon]);
 
  

% set up Legendre pseusospectral computation elements 
CGL_x = -cos(pi*(0:N)'/N); % Chebyshev-Gauss-Labatto Points
[LGL_x,w,P]=lglnodes(N);  % Legendre-Gauss_Lobatto
LGL_x=flipud(LGL_x);  % LGL_x: Legendre-Gauss-Lobatta grid points
w=flipud(w);          %     w: LGL quadrature weights
M = diag(w); Minv = inv(M);  % M: Mass matrix; Minv: M^{-1}

[dpy, ~] = lepoly(N,LGL_x);    % First derivative of legendre polynomials dy
[dpx, ~] = lepoly(N,CGL_x);    % First derivative of Chebyshev polynomials dx
lc_p = ((1+LGL_x) .* dpx(:))/(2*dpy(N+1)); % Legendre intepolation polynomial using Chebyshev points
lc_m = ((1-LGL_x) .* dpx(:))/(2*dpy(N+1)); % Legendre intepolation polynomial using Chebyshev points

S = LegIntPol(CGL_x,LGL_x,N); % Legendre Interpolation Polynomial
S = S';          % S is the lgl to cgl transform matrix 
T = inv(S);     % T is the cgl to lgl transform matrix
                      
% set up differentiation grid   points 
DiffMat_lgl=collocD(LGL_x);         % D: Differentiation matrix
[DiffMat_cgl,~] = cheb(N); 
DiffMat_cgl = -DiffMat_cgl;% CGL Differentiation matrix

% setup grid points in physical space and compute transformation metrics
x = xmin + (xmax-xmin)/(CGL_x(N+1)-CGL_x(1)) * (CGL_x(1:N+1)+1);
% y = xmin + (xmax-xmin)/(LGL_x(N+1)-LGL_x(1)) * (LGL_x(1:N+1)+1);
dx_dxi = chebfft1(x,1); dxi_dx = 1./dx_dxi; Jacobian= dx_dxi; 
% dx_dxi_lgl = DiffMat_lgl*y; dxi_dx_lgl = 1./dx_dxi; Jacobian_lgl= dx_dxi; 

% compute cross section area
Area     =  eval(AreaStr);


% setup tau parameters based on the imposed BCs
% TauPhi_m = c_tau * 1/w(1);
% TauPhi_p = c_tau * 1/w(N+1); 
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
% stand_dt = CFL * 0.5/(max(DiffIon)) *dx_min^delta;
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
% %     TauPhi_m = c_tau * 1/w(1);   % constant tau at xmin
% %     TauPhi_p = c_tau * 1/w(N+1); % constant tau at xmax
% %     e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
% %     eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';
% %     A = diag(Jacobian .* dxi_dx .* dielec .* Area .* dxi_dx); 
% % 
% %     B_p = - Minv * (transpose(DiffMat) - TauPhi_p * eye(size(DiffMat))) * A;
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
% disp(A)
    B_p = TauPhi_p * Minv; %M \ A;
    bp = S*B_p*eN;
    B_m = TauPhi_m * Minv; %M \ A;
    bm = S*B_m*e0;

%     disp(S*A*T - S*T)
    L = - (DiffMat_lgl * A * DiffMat_lgl ...
        - B_p * (AlphaBetaPhi(2,1) * I_p + AlphaBetaPhi(2,2) ...
          * I_p * diag(dxi_dx) * DiffMat_lgl) ...
        - B_m * (AlphaBetaPhi(1,1) * I_m - AlphaBetaPhi(1,2) ...
          * I_m * diag(dxi_dx) * DiffMat_lgl));
% %     L = - (DiffMat_cgl * S * A * T * DiffMat_cgl ...
% %         - S * B_p * (AlphaBetaPhi(2,1) * I_p + AlphaBetaPhi(2,2) ...
% %           * I_p * diag(dxi_dx) * T * DiffMat_cgl * S) * T ...
% %         - S * B_m * (AlphaBetaPhi(1,1) * I_m - AlphaBetaPhi(1,2) ...
% %           * I_m * diag(dxi_dx) * T * DiffMat_cgl * S) * T);
    G = M*L ;
%     disp (max(max(G-G')))
%     disp(eig(G))

% % % end 
SolOperPoisson= inv(S*L*T);

% set up initial condition
time = tstart; [Conc0] = ConcPhiInit(x,IonBc);
u=[Conc0(:,1);Conc0(:,2)];

tcount = 0;
% utemp = u0;
Energy = 0;
Current_t = 0;
Energy_OLD = 0;
Current_OLD = 0;
% % E_NP_t = 0;
% % E_Poisson_t = 0;
% % E_Steric_t =0;
% % E_HOT_t = 0;
% phi_t = 0;


% Solution figure handle
phi=zeros(N+1,1); Current_x=zeros(size(x)); I = zeros(N+1,1);I_ext=I;I_PNP=I;
dEdt=zeros(size(x)); Energy_t=0; Energy_err=0; Current_total=0;
VisualProfile(Conc0(:,1), Conc0(:,2), phi, Current_x, I_ext, I_PNP, ...
    Current_t, Energy, x, StericG, HOT, dt, tcount, 0, 1, ...
    Energy_t, Energy_err, Current_total);
pause

A = Jacobian .* Area .* dxi_dx;

% Set up output folder and file name
StrSteric = [num2str(StericG(1,1)) ' ' num2str(StericG(1,2)) ' ' ...
        num2str(StericG(2,1)) ' '  num2str(StericG(2,2))];
StrHOT = [num2str(HOT(1,1)) ' ' num2str(HOT(1,2)) ' ' ...
        num2str(HOT(2,1)) ' ' num2str(HOT(2,2))];
txtfolder = ['./Results/N' num2str(N) '/tend=' ...
        num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT '/textFile/'];
    
if ~exist(txtfolder, 'dir')
   mkdir(txtfolder);
%    warning (['creating new folder' txtfolder]);
end

filename = ['C1-=' num2str(Conc0(1,1)) ' ' ' C1+=' num2str(Conc0(end,1)) ' '...
    ' C2-=' num2str(Conc0(1,2)) ' ' ' C2+=' num2str(Conc0(end,2)) ' ' ...
    ' Phi-=' num2str(PhiBc(1)) ' ' ' Phi+=' num2str(PhiBc(2))];
outputInfo = [txtfolder '/' filename '.txt'];
InfotxtID = fopen(outputInfo, 'w');
outputName = [txtfolder filename '.dat'];
fileID = fopen(outputName, 'w');
% ENERGY = repmat(Energy_t, length(x),1);
fprintf(fileID, 'TITLE = "2 species of PNP-Steric with high order term" \n');
fprintf(fileID, ['VARIABLES = "X" "Conc1" "Conc2" "phi" "Current_x"' ...
    '"I_ext" "I_PNP" "CURRENT" "ENERGY" "\n']);
% fprintf(fileID, ['ZONE T="' num2str(time) '", I=' num2str(N+1) '\n']);
%     OUTPUT = [x Conc0(:,1) Conc0(:,2) phi I_total ENERGY];
%     fprintf(fileID, '%f %f %f %f %f %f \n', OUTPUT');

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
new_u = new_utemp(end,:)';
dCdt_err = max(max(abs(new_u-u)));
fprintf('relative err of Conc = %d \n', dCdt_err);
fprintf(InfotxtID, 'relative err of Conc = %d \n', dCdt_err);
u = new_u;

tcount = tcount+1;
time=time+dt;

Conc1 = u(1:N+1); Conc2 = u(N+2:2*N+2);
Conc(1:N+1,1) = Conc1(:,end);
Conc(1:N+1,2) = Conc2(:,end);
% testC1 = [Conc1<0];
% testC2 = [Conc2<0];
%disp([nnz(testC1), nnz(testC2)]);

fprintf('dEdt = %d \n', dEdt);
fprintf(InfotxtID, 'dEdt = %d \n', dEdt);

% % check the steady-state %%+++++++++++++++++++++++++++++++++++++++++++++++
dConcdt = max(max(abs(Current_OLD - Current)));
fprintf('relative error of Current = %d \n', dConcdt);
fprintf(InfotxtID, 'relative error of Current = %d \n', dConcdt);
fprintf(InfotxtID, 'dCdt = %d \n', max(max(abs(dCdT))));
fprintf('dCdt = %d \n', max(max(abs(dCdT))));

% % % % % % % % %Electrio current
Current_x = repmat(ValIon(1),N+1,1) .* Current(:,1) + repmat(ValIon(2),N+1,1) .* Current(:,2);  
Current_total = sum(w.*Current_x);
Current_t = [Current_t; Current_total];
fprintf('Current_total = %d \n', Current_total);
fprintf(InfotxtID, 'Current_total = %d \n', Current_total);

A = Jacobian .* Area .* dxi_dx;
Energy_t = sum(w .* A .* En_total) ;   
Energy_err = Energy_t - Energy_OLD;
Energy_OLD = Energy_t;
Energy = [Energy; Energy_t];
fprintf('Energy_err = %d \n', Energy_err);
fprintf('Energy_t = %d \n', Energy_t); 
fprintf(InfotxtID,'Energy_err = %d \n', Energy_err);
fprintf(InfotxtID,'Energy_t = %d \n', Energy_t); 

fprintf(' \n');
fprintf(InfotxtID, ' \n');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
Qs = Qsource(x,time);
if mod(tcount,savesteps) == 0
   matfolder = ['./Results/N' num2str(N) '/tend=' ...
        num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT ...
        '/Conc1=' num2str(Conc0(1,1)) ' Conc2=' num2str(Conc0(1,2)) '/matFile/'...
        '/Q=' num2str(Qs(1)) ' phi_l=' num2str(PhiBc(1)) ' phi_r=' num2str(PhiBc(2)) ];
   
   if ~exist(matfolder, 'dir')
       mkdir(matfolder);
   %    print("creating new folder", matfolder);
   end
   
   save([matfolder num2str(tcount/savesteps),'.mat'], ...
       'Conc','phi','x', 'Current_x', 'I_ext', 'I_PNP', 'Energy', 'Current_t', 'time');
end   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Output results to .txt file
    ENERGY = repmat(Energy_t, length(x),1);
    CURRENT = repmat(Current_total, length(x),1);
    fprintf(fileID, ['ZONE T="' num2str(time) '", I=' num2str(N+1) '\n']);
    A = [x Conc(:,1) Conc(:,2) phi Current_x I_ext I_PNP CURRENT ENERGY];
    fprintf(fileID, '%d %d %d %d %d %d %d %d %d \n', A');
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
if mod(tcount,savesteps)==0
    Plotout = ['./Results/N' num2str(N) '/tend=' ...
        num2str(tend) '/StericG=' StrSteric '/HOT=' StrHOT ...
        '/Conc1=' num2str(Conc0(1,1)) ' Conc2=' num2str(Conc0(1,2)) '/PlotsProfile/' ...
        '/Q=' num2str(Qs(1)) ' phi_l=' num2str(PhiBc(1)) ' phi_r=' num2str(PhiBc(2))];
        
if ~exist(Plotout, 'dir')
        mkdir(Plotout);
end
h=figure(5);
VisualProfile(Conc(:,1), Conc(:,2), phi, Current_x, I_ext, I_PNP, Current_t, Energy, ...
    x, StericG, HOT, dt, tcount, time, 1, Energy_t, Energy_err, Current_total);

saveas(h, [Plotout  '/' num2str(tcount/savesteps) '.eps'], 'psc2');
end
end

    fclose('all');

end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
function dudt=pnp1d(t,u)
global x N M lc_p lc_m S T
global Valence bp bm Area DiffIon AlphaBetaIon
global IonBc PhiBc HOT ChemBc dCdxBc 
global kBxT TotNumIon StericG ValIon e_unit
global Jacobian SolOperPoisson dxi_dx TauIon
global ddxAdCdxBc %DiffMat err
% % values from PNP_Steric
global phi I_PNP I_ext Current dEdt En_total dCdT

Conc = zeros([N+1 TotNumIon]);
Conc(1:N+1,1) = u(1:N+1);
Conc(1:N+1,2) = u(N+2:2*N+2);

rhs_b = zeros(size(x)); phi = zeros(size(x)); Qs = zeros(size(x));
dEdc = zeros(size(Conc)); Mob  = zeros(size(Conc));  
DivCurrent = zeros(size(Conc));

e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';

Qs = Qsource(x,t);
rhs_b = Qs + Conc(:,1)*Valence(1) + Conc(:,2)*Valence(2);
rhs_b = Jacobian .* rhs_b .* Area + bp * PhiBc(2) + bm * PhiBc(1);

phi = SolOperPoisson * rhs_b;


A = Jacobian .* Area .* dxi_dx;
e0 = zeros(N+1,1); e0(1)=1; eN = zeros(N+1,1); eN(N+1)=1;

% Compute dE/dc
dEdc = kBxT * (log(Conc) + 1) ...
    + repmat(Valence,N+1,1) .* repmat(phi,1,TotNumIon);
Mob = repmat(DiffIon,N+1,1)/kBxT .* Conc;

Current_PNP = - Mob .* (diag(dxi_dx)*chebfft1(dEdc,2));

% disp(chebfft1(dEdc,2))
% pause
% Current_PNP = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);
% % % %Electrio current
I_PNP = sum(repmat(ValIon,N+1,1) .* Current_PNP,2);

% add steric effect term
dEdc_ext(1:N+1,1) = sum(repmat(StericG(1,1:2),N+1,1) .* Conc (1:N+1,1:2),2);
dEdc_ext(1:N+1,2) = sum(repmat(StericG(2,1:2),N+1,1) .* Conc (1:N+1,1:2),2);
    
%add High order term
dEdc_ext(1:N+1,1) = dEdc_ext(1:N+1,1) ...
    + HOT(1,1) * (- diag(1./A) * chebfft2(diag(A), Conc(1:N+1,1),1) ...
                + lc_p ./ M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx)...
                * chebfft1(Conc(1:N+1,1),1) - dCdxBc(2,1)) ...
                - lc_m ./ M( 1 , 1 ) * dxi_dx( 1 ) * (e0' * diag(dxi_dx)...
                * chebfft(Conc(1:N+1,1)) - dCdxBc(1,1))) ...
    + HOT(1,2) * (- diag(1./A) * chebfft2(diag(A), Conc(1:N+1,2),1) ...
                + lc_p ./ M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx)...
                * chebfft1(Conc(1:N+1,2),1) - dCdxBc(2,2)) ...
                - lc_m ./ M( 1 , 1 ) * dxi_dx( 1 ) * (e0' * diag(dxi_dx)...
                * chebfft1(Conc(1:N+1,2),1) - dCdxBc(1,2))); 
            
dEdc_ext(1:N+1,2) = dEdc_ext(1:N+1,2) ...
    + HOT(2,1) * (- diag(1./A) * chebfft2(diag(A),Conc(1:N+1,1),1) ...
                + lc_p ./ M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) ...
                * chebfft1(Conc(1:N+1,1),1) - dCdxBc(2,1)) ...
                - lc_m ./ M( 1 , 1 ) * dxi_dx( 1 ) * (e0' * diag(dxi_dx) ...
                * chebfft1(Conc(1:N+1,1),1) - dCdxBc(1,1))) ...
    + HOT(2,2) * (- diag(1./A) * chebfft2(diag(A), Conc(1:N+1,2),1) ...
                + lc_p ./ M(end,end) * dxi_dx(N+1) * (eN' * diag(dxi_dx) ...
                * chebfft1(Conc(1:N+1,2),1) - dCdxBc(2,2)) ...
                - lc_m ./ M( 1 , 1 ) * dxi_dx( 1 ) * (e0' * diag(dxi_dx) ...
                * chebfft1(Conc(1:N+1,2),1) - dCdxBc(1,2)));             

% compute current J = -  D/kBT * Conc * D_x(dEdc) )
Current_ext = - Mob .* (diag(dxi_dx)*chebfft1(dEdc_ext,2));
I_ext = zeros(size(sum(repmat(ValIon,N+1,1) .* Current_ext,2)));

dEdc = dEdc+dEdc_ext;
Current = - Mob .* (diag(dxi_dx)*chebfft1(dEdc,2));
% Current = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);
% % % %Electrio current
% % % I_total = repmat(ValIon(1),N+1,1) .* Current(:,1) + repmat(ValIon(2),N+1,1) .* Current(:,2);  

%Compute divergence current DivCurrent(1:N+1,1:TotNumIon)
DivCurrent = - diag(1./(Jacobian.*Area)) ...
    * (chebfft1(diag(Jacobian .* Area .* dxi_dx) * Current,2));
% DivCurrent = - diag(1./(Jacobian.*Area)) ...
%     * (DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current));

for IonNum=1:TotNumIon
    
    % modify flux function with BCs  
    % evaluate exact u for imposing Dirichlet bc for u
%     if ChemBc == 0
%     if IonNum == 1
% % %         g_m_Cpp = HOT(1,1) * (- 1/Area(1) * ddxAdCdxBc(1,1)) ...
% % %             + HOT(1,2) * (- 1/Area(1) * ddxAdCdxBc(1,2));
%        g_m_Cpp = HOT(1,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
%                + HOT(1,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2));
%         
%         g_m = DiffIon(IonNum)/kBxT * IonBc(1,1) * (kBxT * ...
%            (log(IonBc(1,1)) + 1) + Valence(1) * PhiBc(1) + ...
%            sum(StericG(1,1:2) .* IonBc(1,1:2)) + g_m_Cpp);
%        
% % %        g_p_Cpp = HOT(1,1) * (- 1/Area(N+1) * ddxAdCdxBc(2,1)) ...
% % %            + HOT(1,2) * (- 1/Area(N+1) * ddxAdCdxBc(2,2));
%        g_p_Cpp = HOT(1,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
%                + HOT(1,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2)); 
%        
%        g_p = DiffIon(IonNum)/kBxT * IonBc(2,1) * (kBxT * ...
%            (log(IonBc(2,1)) + 1) + Valence(1) * PhiBc(2) + ...
%            sum(StericG(1,1:2) .* IonBc(2,1:2)) + g_p_Cpp);
% 
%     end 
% 
%     if IonNum == 2
% %         g_m_Cpp = HOT(2,1) * (- 1/Area(1) * ddxAdCdxBc(1,1)) ...
% %             + HOT(2,2) * (- 1/Area(1) * ddxAdCdxBc(1,2));
%         g_m_Cpp = HOT(2,1) * (- 1/(e0'*A) * ddxAdCdxBc(1,1)) ...
%                 + HOT(2,2) * (- 1/(e0'*A) * ddxAdCdxBc(1,2));
%       
%         g_m = DiffIon(IonNum)/kBxT * IonBc(1,2) * (kBxT * ...
%            (log(IonBc(1,2)) + 1) + Valence(2) * PhiBc(1) + ...
%            sum(StericG(2,1:2) .* IonBc(1,1:2)) + g_m_Cpp);
% 
%         
% %        g_p_Cpp = HOT(2,1) * (- 1/Area(N+1) * ddxAdCdxBc(2,1)) ...
% %            + HOT(2,2) * (- 1/Area(N+1) * ddxAdCdxBc(2,2));
%        g_p_Cpp = HOT(2,1) * (- 1/(eN'*A) * ddxAdCdxBc(2,1)) ...
%                + HOT(2,2) * (- 1/(eN'*A) * ddxAdCdxBc(2,2));
%         
%        g_p = DiffIon(IonNum)/kBxT * IonBc(2,2) * (kBxT * ...
%            (log(IonBc(2,2)) + 1) + Valence(2) * PhiBc(2) + ...
%            sum(StericG(2,1:2) .* IonBc(2,1:2)) + g_p_Cpp);
%     end
%     else
        g_m = 0; g_p = 0; 
%     end
    r_m = Area( 1 ) * (AlphaBetaIon(1,1,IonNum) ...
        * (Mob( 1 ,IonNum) * dEdc( 1 ,IonNum) - g_m) ...
        - AlphaBetaIon(1,2,IonNum) * (-Current( 1 ,IonNum) - 0) );

    r_p = Area(N+1) * (AlphaBetaIon(2,1,IonNum) ...
        * (Mob(N+1,IonNum) * dEdc(N+1,IonNum) - g_p) ...
        + AlphaBetaIon(2,2,IonNum) * (-Current(N+1,IonNum) - 0) );
    
    % add penalty BC to flux

    DivCurrent(:,IonNum) = DivCurrent(:,IonNum) + TauIon(1,IonNum) * lc_m .* r_m ...
                            + TauIon(2,IonNum) * lc_p .* r_p;

%     DivCurrent( 1 ,IonNum) = DivCurrent( 1 ,IonNum) - TauIon(1,IonNum) * r_m;
%     DivCurrent(N+1,IonNum) = DivCurrent(N+1,IonNum) - TauIon(2,IonNum) * r_p;
    
end

dEdt = sum(M * sum(((- chebfft1(diag(Jacobian .* Area .* dxi_dx) * Current,2)) .* dEdc),2));
% dEdt = sum(M * sum(((- DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current)) .* dEdc),2));
%%%% compute energy ====================================================
E_NP = kBxT * sum(Conc.*log(Conc),2); 
E_Poisson = 1/2 * e_unit * (Qs + sum(repmat(ValIon,N+1,1).*Conc,2)) .* phi;
En_total = E_NP + E_Poisson;% + E_Steric + E_HOT;

dudt = [DivCurrent(:,1);DivCurrent(:,2)];
dCdT = dudt;
% err = max(abs(dudt));
% format short e
%   disp([t err])
% pause;
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function w = chebfft1(v,n)
% w = zeros(size(v));
% for i = 1:n
%     w(i)=chebfft(-v(i));
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function d2u = chebfft2(k,u,n)
% % DiffMatt * K * DiffMatt * U 
% % d/dx (k du/dx)
% 
% du = chebfft1(u,n);
% d2u = chebfft1(k*du,n);
% end












