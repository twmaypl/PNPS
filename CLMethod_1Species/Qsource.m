% External source of electrical potencial
% Update on 2015/08/13

function [Qs] = Qsource(x,t,C)

global e_unit ValIon
% global DiffIon DiffMat_cgl phi_ex

zp = ValIon*e_unit;
% Qs = ones(size(x))*-1e-4;
Qs = zeros(size(x));
% Qs = -chebfft1(chebfft1(phi_ex(C),1),1) - zp * C;
end

