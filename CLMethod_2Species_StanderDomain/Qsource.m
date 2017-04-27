% External source of electrical potencial
% Update on 2016/01/11

function [Qs] = Qsource(x,t,C)

global e_unit 
global ValIon phi_ex

zp = ValIon(1)*e_unit;
zn = ValIon(2)*e_unit;
% Qs = - chebfft1(chebfft1(phi_ex(C(:,1)),1),1) - zp * C(:,1) - zn * C(:,2);
Qs = zeros(size(x));
end

