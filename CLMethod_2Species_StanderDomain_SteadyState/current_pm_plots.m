clear all; close all;
istart = 5000; iend = 6000; istep=1;
inum = iend-istart+1;
Current_p = zeros(1,inum); Current_m = zeros(1,inum);
for i=1:inum
% %     load(['./1Results/N100/tend=20/StericG=1 2 1 4/HOT=0.001 0 0.01 0' ...
% %     '/333Conc1=0.001 Conc3=3/matFile/Q=0 z1=1 z3=2 phi_l=-0.0001 ' ...
% %     'phi_r=-0.0001/' num2str(istart-1+i) '.mat']);
    load (['./Results/N100/tend=12/StericG=1 2 1 4/HOT=0.001 0 0.01 0/' ...
        'Conc1=0.001 Conc3=3/matFile/Q=0 z1=1 z3=2 phi_l=-0.001 phi_r=-0.001/'...
        num2str(istart-1+i) '.mat']);
Current_p(i) = Current_x(end);
Current_m(i) = Current_x(1);
end
X = [istart/1000:0.001*istep:iend/1000];
figure(1),
plot(X,Current_p(1:istep:end),'LineWidth', 4);
title('x=1'); xlabel('time'); ylabel('Current');
figure(2),
plot(X,-Current_m(1:istep:end), 'LineWidth', 4);
title('x=-1'); xlabel('time'); ylabel('Current');
figure(3),
plot(X,Current_p(1:istep:end),'b', 'LineWidth', 4);
hold on,
plot(X,-Current_m(1:istep:end),'r', 'LineWidth', 4);
hold off
title('blue x=1, red x=-1'); xlabel('time'); ylabel('Current');