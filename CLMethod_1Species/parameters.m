% Input parameter for PNP-Steric simulation
% Update on 2015/09/03

function [input]=parameters
input.c_tau=2.0;                       %scaling parameter for adjjusting penalty parameters tau
input.NN=80;                           % Number of points  
input.CFL=0.35; 
input.savesteps = 1;                   % Output every # of steps
input.Euler = 0;                       % If 1, using Euler method, otherwise using ode15s
input.dt = 1e-3;                       % dt for ODE15s

input.xmin     = -1;                   % Left bound
input.xmax     = 1;                    % Right bound

input.tstart = 0.0;                    % Start time
input.tend   = 1e-1;                   % End time

% Physical parameters
input.kB=1.0;
input.Temp=1.0;
input.e_unit =  1.0;
input.dielec = 1.0;
input.epslion = 1E-13;

% Number of species
input.Ns = 1;       

% Diffusion constants;
Dp = 1;
input.DiffIon =Dp;%*10^(-3);

% Valence constants
zp = 1.0;
input.ValIon=zp;

% Steric effect parameters
input.StericG= 8 ;

% High order term parameters
input.HOT= 0;%1e-4;             

% Boundary Condition
input.AlphaBetaIon = zeros(2,2);
input.AlphaBetaPhi = zeros(2,2);
input.AlphaBetaIon(1:2,1:2)=[1 0; 1 0];
input.AlphaBetaPhi(1:2,1:2) = [1 1e-1; 1 1e-1];

input.AreaStr='1.0*ones(size(x))';
% % % % Boundary Conditions
ion1bc_l = 3;%exp(input.xmin); 
ion1bc_r = 3;%exp(input.xmax);

dCdx1bc_l = 0;
dCdx1bc_r = 0;

phi_l = -1e-4; %-1e-4;%input.xmin; 
phi_r = -1e-4; %1e-4;%input.xmax;

input.IonBc=[ion1bc_l; ion1bc_r];
input.dCdxBc = [dCdx1bc_l; dCdx1bc_r];
input.PhiBc=[phi_l phi_r]; 



 
