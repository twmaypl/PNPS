% Input parameter for PNP-Steric simulation
% Update on 2016/01/11

function [input]=parameters
input.c_tau=2.0;                       %scaling parameter for adjjusting penalty parameters tau
input.NN=80;                           % Number of points  
input.CFL=0.35/150;
input.ep = 0.001;                      % Parameter for initial condition
input.savesteps = 1;                   % Output every # of steps
input.dt = 1e-4;%0.1;
input.Euler = 0;                       % if 1, using Euler, otherwise using ode15s

input.xmin     = -1;                    % Left bound
input.xmax     = 1; %1.5;              % Right bound
% SpaceDomain =[xmin xmax];

input.tstart = 0.0;                    % Start time
input.tend   = 5e-4;%2.0;             % End time


input.kB=1.0;
input.Temp=1.0;
input.e_unit =  1.0;
input.epslion = 1E-13;


% Number of species
input.Ns = 2;       

% Diffusion constants
input.DiffIon = zeros(1,input.Ns);
Dp = 1; Dn = 1;
input.DiffIon =[Dp Dn];%*10^(-3);

% Valence constants
input.ValIon = zeros(1,input.Ns);
zp = 2.0; zn= 2.0;
input.ValIon=[zp zn];

% Steric effect parameters
input.StericG = [1 2; 2 1];%[1 4; 4 1];


% High order term parameters
input.HOT = [5e-3 0; 0 5e-3];%[0 0; 0 0];%
     

% Boundary Condition
input.AlphaBetaIon = zeros(2,2,input.Ns);
input.AlphaBetaPhi = zeros(2,2);
input.AlphaBetaIon(1:2,1:2,1)=[1 0; 1 0];
input.AlphaBetaIon(1:2,1:2,2)=[1 0; 1 0];
% input.AlphaBetaPhi(1:2,1:2) = [1 1; 1 1];
input.AlphaBetaPhi(1:2,1:2) = [1 1e-1; 1 1e-1];

ion1bc_l = 3;%exp(input.xmin); 
ion1bc_r = 3;%exp(input.xmax);
ion2bc_l = 3;%exp(input.xmin); 
ion2bc_r = 3;%exp(input.xmax);

dCdx1bc_l = 0;
dCdx1bc_r = 0;
dCdx2bc_l = 0;
dCdx2bc_r = 0; 

phi_l = -1e-4; %-1e-4;%input.xmin; 
phi_r = -1e-4; %1e-4;%input.xmax;

input.IonBc=[ion1bc_l ion2bc_l; ion1bc_r ion2bc_r];
input.dCdxBc = [dCdx1bc_l dCdx2bc_l; dCdx1bc_r dCdx2bc_r];
input.PhiBc=[phi_l phi_r]; 

 
