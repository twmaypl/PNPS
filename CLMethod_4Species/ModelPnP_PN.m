% The solver of Poisson-Nerst-Planck equation couple wiht steric effect
% update log:
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
%     u = u + dt du/dt
%   

function [err] = ModelPnP_PN (IonBc, PhiBc)


global kB Temp e_unit dielec savesteps
global DiffIon ValIon StericG
global AlphaBetaIon AreaStr 
global xmin xmax tstart tend
global N CFL c_tau

StericG 
%Determine the total number of Ion species
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
[LGL_x,w,~]=lglnodes(N);  
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
TauPhi_m = c_tau * 1/w(1);   % constant tau at xmin
TauPhi_p = c_tau * 1/w(N+1); % constant tau at xmax

% TauIon(1:2,:) = 1 / (4 * AlphaBetaIon(1:2,1,:) * w( 1 )...
%              + AlphaBetaIon(1:2,2,:))/ w( 1 ); % general BC at x=-1
%          
TauIon(1,:) = 1 / (4 * AlphaBetaIon(1,1,:) * w(  1  ) * Jacobian(1)...
             + AlphaBetaIon(1,2,:))/ w(  1  ) ; % general BC 
TauIon(2,:) = 1 / (4 * AlphaBetaIon(2,1,:) * w( N+1 ) * Jacobian(N+1)...
             + AlphaBetaIon(2,2,:))/ w( N+1 ) ; % general BC at x=-1
         
                   
for IonNum = 1:TotNumIon
 
    if AlphaBetaIon(1,1,IonNum) == 1
       TauIon(1,IonNum) = TauIon(1,IonNum) * c_tau;
    end
    if AlphaBetaIon(2,1,IonNum) == 1
       TauIon(2,IonNum) = TauIon(2,IonNum) * c_tau;
    end
end

% set up poission solution operator 
    disp('weakly solution');
    
    e0 = zeros(N+1,1); e0(1)=1; I_m = e0 * e0';  
    eN = zeros(N+1,1); eN(N+1)=1; I_p = eN * eN';
    A = diag(Jacobian .* dxi_dx .* dielec .* Area .* dxi_dx); 

    B_p = - Minv * (transpose(DiffMat) - TauPhi_p * eye(size(DiffMat))) * A;
    bp = B_p*eN;
    B_m = + Minv * (transpose(DiffMat) + TauPhi_m * eye(size(DiffMat))) * A;
    bm = B_m*e0;
    L = - (DiffMat * A * DiffMat - B_p * I_p - B_m * I_m);
    G = M*L;
   SolOperPoisson= inv(G); 
   
testpara = 0;
while testpara ==0 

% set up time step and parameters
dx_min=abs(x(2)-x(1));
disp(CFL)
dt = CFL * 0.5/(max(DiffIon)) *dx_min^2
    
% set up initial condition
time = tstart; [Conc] = ConcPhiInit(x, IonBc);


% Solution figure handle
% % phi=zeros(N+1,1); I = zeros(N+1,1);
% % VisualProfile(Conc(:,1),Conc(:,2),Conc(:,3),Conc(:,4),phi,I,x,StericG,time,1);
% % pause
tcount = 0;

%Check Point
AA = nnz(sum(repmat(ValIon, N+1, 1) .* Conc, 2));
if AA ~=0
    warning('Not neuture everywhere')
end


while (time  < tend) && (testpara==0) %
   if time+dt >= tend
       dt = tend-time;
   end
   
   rhs_b = zeros(size(x)); phi = zeros(size(x)); Qs = zeros(size(x));
   dEdc = zeros(size(Conc));
   Mob  = zeros(size(Conc));  DivCurrent = zeros(size(Conc));
   
   % solve poission equation for phi
    % set up right hand side
    Qs = Qsource(x,time);
%     disp(max(Qs-exp(-x)))
    rhs_b = Qs + sum(Conc .* repmat(Valence,length(Conc(:,1)),1),2);
    rhs_b = Jacobian .* rhs_b .* Area + bm * PhiBc(1) + bp * PhiBc(2);
    phi = SolOperPoisson * M * rhs_b;

   % Compute dE/dc
   dEdc = kBxT * (log(Conc) + 1) ...
        + repmat(Valence,N+1,1) .* repmat(phi,1,TotNumIon);
   % add steric effect term
   dEdc(1:N+1,1) = dEdc(1:N+1,1) + sum(repmat(StericG(1,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
   dEdc(1:N+1,2) = dEdc(1:N+1,2) + sum(repmat(StericG(2,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
   dEdc(1:N+1,3) = dEdc(1:N+1,3) + sum(repmat(StericG(3,1:4),N+1,1) .* Conc (1:N+1,1:4),2);
   dEdc(1:N+1,4) = dEdc(1:N+1,4) + sum(repmat(StericG(4,1:4),N+1,1) .* Conc (1:N+1,1:4),2);

   % compute current J = -  D/kBT * Conc * D_x(dEdc) )
   Mob = repmat(DiffIon,N+1,1)/kBxT .* Conc;
   
   Current = - Mob .* (diag(dxi_dx)*DiffMat*dEdc);  
   
   if mod(tcount,savesteps) == 0 || time == tend
   %Electrio Current
   I_total = sum(repmat(ValIon,N+1,1) .* Current,2);
%    I_total = repmat(ValIon(1),N+1,1) .* Current(:,1) + ...
%        repmat(ValIon(2),N+1,1) .* Current(:,2);
   end
                 
   %Compute divergence current DivCurrent(1:N+1,1:TotNumIon)
   DivCurrent = - diag(1./(Jacobian.*Area)) ...
       * (DiffMat * (diag(Jacobian .* Area .* dxi_dx) * Current));
 
   for IonNum=1:TotNumIon
       
    if IonNum == 1
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,1) * (kBxT * ...
           (log(IonBc(1,1)) + 1) + Valence(1) * PhiBc(1) + ...
           sum(StericG(1,1:4) .* IonBc(1,1:4)));
           
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,1) * (kBxT * ...
           (log(IonBc(2,1)) + 1) + Valence(1) * PhiBc(2) + ...
           sum(StericG(1,1:4) .* IonBc(2,1:4)));            
    end 

    if IonNum == 2
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,2) * (kBxT * ...
           (log(IonBc(1,2)) + 1) + Valence(2) * PhiBc(1) + ...
           sum(StericG(2,1:4) .* IonBc(1,1:4)));
       
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,2) * (kBxT * ...
           (log(IonBc(2,2)) + 1) + Valence(2) * PhiBc(2) + ...
           sum(StericG(2,1:4) .* IonBc(2,1:4)));
    end
    
    if IonNum == 3
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,3) * (kBxT * ...
           (log(IonBc(1,3)) + 1) + Valence(3) * PhiBc(1) + ...
           sum(StericG(3,1:4) .* IonBc(1,1:4)));
       
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,3) * (kBxT * ...
           (log(IonBc(2,3)) + 1) + Valence(3) * PhiBc(2) + ...
           sum(StericG(3,1:4) .* IonBc(2,1:4)));
    end
    
    if IonNum == 4
        g_m = DiffIon(IonNum)/kBxT * IonBc(1,4) * (kBxT * ...
           (log(IonBc(1,4)) + 1) + Valence(4) * PhiBc(1) + ...
           sum(StericG(4,1:4) .* IonBc(1,1:4)));
       
       g_p = DiffIon(IonNum)/kBxT * IonBc(2,4) * (kBxT * ...
           (log(IonBc(2,4)) + 1) + Valence(4) * PhiBc(2) + ...
           sum(StericG(4,1:4) .* IonBc(2,1:4)));
%        disp ([g_m g_p])
    end
     
    r_m = Area( 1 ) * (AlphaBetaIon(1,1,IonNum) ...
        * (Mob( 1 ,IonNum) * dEdc( 1 ,IonNum) - g_m) ...
        - AlphaBetaIon(1,2,IonNum) * (-Current( 1 ,IonNum) - 0) );

    r_p = Area(N+1) * (AlphaBetaIon(2,1,IonNum) ...
        * (Mob(N+1,IonNum) * dEdc(N+1,IonNum) - g_p) ...
        + AlphaBetaIon(2,2,IonNum) * (-Current(N+1,IonNum) - 0) );
   
    % add penalty BC to flux
    DivCurrent( 1 ,IonNum) = DivCurrent( 1 ,IonNum) - ...
        1/Area( 1 ) .* 1/Jacobian( 1 ) * TauIon(1,IonNum) * r_m;
    DivCurrent(N+1,IonNum) = DivCurrent(N+1,IonNum) - ...
        1/Area(N+1) .* 1/Jacobian(N+1) * TauIon(2,IonNum) * r_p;

   end

    

% pause;
   % update solution in time u = u + dt* DivCurrent; 
   Conc = M*Conc + (dt * M * DivCurrent);
   Conc = M\(Conc);
% %    Conc = Conc + dt * DivCurrent;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if mod(tcount,savesteps) == 0
   matfolder = ['./Results/matFile/N=', num2str(N), '/tend=', ...
    num2str(tend), '/StericG=', num2str(reshape(StericG,1,...
    length(StericG(:,1))*length(StericG(1,:))))]; %num2str(StericG(1,1)), ' ', ...
%        num2str(StericG(1,2)), ' ', num2str(StericG(2,1)), ' ', ...
%        num2str(StericG(2,2))];
   
   if ~exist(matfolder, 'dir')
       mkdir(matfolder);
   %    print("creating new folder", matfolder);
   end
   
   save([matfolder, '/', num2str(tcount/savesteps),'.mat'], ...
       'Conc','phi','x', 'I_total', 'time');
   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    txtfolder = ['./Results/textFile/N', num2str(N), '/tend=', ...
        num2str(tend), '/StericG', num2str(reshape(StericG,1,...
    length(StericG(:,1))*length(StericG(1,:))))]; %num2str(StericG(1,1)),',', ...
%         num2str(StericG(1,2)), ',',num2str(StericG(2,1)), ',', ...
%         num2str(StericG(2,2))];
    
    if ~exist(txtfolder, 'dir')
        mkdir(txtfolder);
    %    print("creating new folder");
    end
    if mod(tcount,savesteps)==0
    outputName = [txtfolder,'/',num2str(tcount/savesteps),'.txt'];
    fileID = fopen(outputName, 'w');
    TIME = repmat(time, length(x), 1);
    A = [x phi Conc(:,1) Conc(:,2) TIME I_total];
    fprintf(fileID, '%f %f %f %f %f %f \n', A');
    fclose('all');
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
   if mod(tcount,savesteps) == 0
      VisualProfile(Conc(:,1),Conc(:,2),Conc(:,3),Conc(:,4),phi,I_total,x,StericG,time,1);
   end
  
% %    if (mod(time,savesteps))
% %        disp (time);
% %    end

   time=time + dt;   
   tcount = tcount + 1;
   

   
   testC = Conc <= 0;
   if (nnz(testC)>0) || (time==tend) 
       testpara=1;
%        disp (time)
%        disp (Conc)
%        return;
   end
    if time == tend
       disp ([time testpara])
    end
end
if (testpara==1) && (time < tend)
    testpara=0;
    CFL = CFL*0.95
    if CFL < 0.35/10
        warning('CFL too small')
        pause
    end
end
end

Plotout = ['./Results/Plot/N=' num2str(N) '/tend=' num2str(tend)];
if ~exist(Plotout, 'dir')
        mkdir(Plotout);
end
h=figure(1),
VisualProfile(Conc(:,1),Conc(:,2),Conc(:,3),Conc(:,4),phi,I_total,x,StericG,time,1);

saveas(h, [Plotout '/StericG=' num2str(reshape(StericG,1,...
    length(StericG(:,1))*length(StericG(1,:)))) '.eps']);%num2str(StericG(1,1)) ',' ...
%         num2str(StericG(1,2)) ',' num2str(StericG(2,1)) ','  ...
%         num2str(StericG(2,2)) '.eps']); 
    err = abs(Conc(1,1)-Bc(1,1));
end


