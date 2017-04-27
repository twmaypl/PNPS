% External source of electrical potencial
% Update on 2015/09/28

function [Qs] = Qsource(x,t)

global e_unit
global DiffIon ValIon

zp = ValIon(1)*e_unit;
zn = ValIon(2)*e_unit;
Qs = zeros(size(x));
% Qs = ones(size(x))*1e+0;
% Qs = -zp * (1.1+exp(x)) -zn * (1.1+exp(-x));
% Qs = zeros(size(x));%-zp * (exp(x)) - zn * (exp(x));
end

