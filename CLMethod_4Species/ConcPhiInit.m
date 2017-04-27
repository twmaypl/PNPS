% Initial Conditions for ions' concentration
% Update on 2015/09/28
% Change log:
%    Version 3
%      4 species, initial concentrations are liner decay from left to right
%
%    Version 2 
%      Due to version 1 lock of the electric neutural inside the
%      computational domain. We change the initial concentrations of both 
%      ions decrease from left boundary to the right. The electric neutural 
%      which satisfy that zp*Conc1 + zn* Conc2 = 0.
%
%    Version 1 
%      Concentration at first and last 5% of the computational domain are
%      function of tanh decay (increase) from boundary condition (nearly to 
%      zero) to nearly to zero (boundary condition). Between these regions
%      are set to be nearly to zero.
%

function [Conc] = ConcPhiInit(z, IonBc)

global e_unit DiffIon ValIon 
global N ep
  
Valence= ValIon * e_unit;

a=z(end) - z(1); %z=0.5*a*(x+1); 
df=0;
R1=IonBc(2,1);   L1=IonBc(1,1);   R2=IonBc(2,2);   L2=IonBc(1,2);
R3=IonBc(2,3);   L3=IonBc(1,3);   R4=IonBc(2,4);   L4=IonBc(1,4);

Conc = zeros([N+1 length(DiffIon)]);

Conc(:,1) = L1+(R1-L1)/a*(z-(z(1)));
Conc(:,2) = L2+(R2-L2)/a*(z-(z(1)));
Conc(:,3) = L3+(R3-L3)/a*(z-(z(1)));
Conc(:,4) = L4+(R4-L4)/a*(z-(z(1)));
% load('./input/7300.mat');%, 'Conc')

%Conc(:,1) = L1*ones(size(z));
%Conc(:,2) = L2*ones(size(z));
%Conc(:,3) = L3*ones(size(z));
%Conc(:,4) = L4*ones(size(z));
% Conc(:,1) = (tanh((z-(0.95-df)*a)*100)+1+ep)/2*(R1-df)+...
%     (1+ep-(tanh((z-(0.05+df)*a)*100)+1)/2)*(L1-df)+df;
% Conc(:,2) = (tanh((z-(0.95-df)*a)*100)+1+ep)/2*(R2-df)+...
%     (1+ep-(tanh((z-(0.05+df)*a)*100)+1)/2)*(L2-df)+df;
% Conc(:,3) = (tanh((z-(0.95-df)*a)*100)+1+ep)/2*(R3-df)+...
%     (1+ep-(tanh((z-(0.05+df)*a)*100)+1)/2)*(L3-df)+df;
% Conc(:,4) = (tanh((z-(0.95-df)*a)*100)+1+ep)/2*(R4-df)+...
%     (1+ep-(tanh((z-(0.05+df)*a)*100)+1)/2)*(L4-df)+df;
      
end
