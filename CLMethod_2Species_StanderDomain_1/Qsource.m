% External source of electrical potencial
% Update on 2015/08/13

function [Qs] = Qsource(x,t,C)

global e_unit 
global DiffIon ValIon phi_ex

zp = ValIon(1)*e_unit;
zn = ValIon(2)*e_unit;
% Qs = ones(size(x))*-1e-4;
% Qs = - chebfft1(chebfft1(phi_ex(C(:,1)),1),1) - zp * C(:,1) - zn * C(:,2);
Qs = zeros(size(x));%-zp * (exp(x)) - zn * (exp(x));
end

