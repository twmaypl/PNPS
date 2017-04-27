% External source of electrical potencial
% Update on 2015/08/13

function [Qs] = Qsource(x,t,C)

global phi_ex DiffMat DiffIon 
global ValIon e_unit 

zp = ValIon*e_unit;
Qs = - (DiffMat * DiffMat * phi_ex(C) + zp * C);
end

