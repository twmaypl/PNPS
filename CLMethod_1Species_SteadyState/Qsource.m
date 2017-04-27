% External source of electrical potencial
% Update on 2015/08/13

function [Qs] = Qsource(x,t,C)

global e_unit DiffMat_cgl
global DiffIon ValIon phi_ex

zp = ValIon*e_unit;
Qs = -chebfft1(chebfft1(phi_ex(C),1),1) - zp * C;
end

