function VisualProfile(Ion1, Ion2, Phi, Current_x, I_ex, I_pnp, Current_t, Energy, ...
    x, StericG, HOT, dt, tcount, time, option, Energy_t, Energy_err, Current_total)
persistent ConcP_FigHandle ConcN_FigHandle PhiFigHandle TimeStrHandle

if option == 1 
    StrSteric = [num2str(StericG(1,1)) ' ' num2str(StericG(1,2)) ' ' ...
        num2str(StericG(2,1)) ' '  num2str(StericG(2,2))];
    StrHOT = [num2str(HOT(1,1)) ' ' num2str(HOT(1,2)) ' ' ...
        num2str(HOT(2,1)) ' ' num2str(HOT(2,2))];
figure(5)

subplot(4,2,[1,2])
ConcP_FigHandle = plot(x,Ion1,'r-',x,Ion1,'k.'); hold on,
ConcN_FigHandle = plot(x,Ion2,'blue-',x,Ion2,'k.'); hold off
set (ConcP_FigHandle,'EraseMode','Xor','MarkerSize',6)
set (ConcN_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('p(red) & n(blue)');
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str) 
drawnow

% % subplot(4,1,2)
% % ConcN_FigHandle = plot(x,Ion2,'k-',x,Ion2,'r.');
% % set (ConcN_FigHandle,'EraseMode','Xor','MarkerSize',6)
% % xlabel('x'); ylabel('n(x,t)'); 
% % drawnow

subplot(4,2,[3,4])
PhiFigHandle = plot(x,Phi,'k-',x,Phi,'r.');
set (PhiFigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('phi(x,t)');
drawnow

subplot(4,2,[5,6])
JFigHandle1 = plot(x,I_pnp,'r-',x,I_pnp,'r.');
% hold on,
% JFigHandle2 = plot(x,I_ex, 'blue-',x,I_ex,'blue.');
% hold on,
% JFigHandle3 = plot(x,Current_x, 'green-',x,Current_x,'green.');
% hold off
set (JFigHandle1,'EraseMode','Xor','MarkerSize',6)
% set (JFigHandle2,'EraseMode','Xor','MarkerSize',6)
% set (JFigHandle3,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('Current');
% title('I_{pnp}-red; I_{ext}-blue; I-green')
drawnow

subplot(4,2,7)
EFH=plot([1:tcount]*dt, Energy(2:tcount+1),'k-', [1:tcount]*dt, Energy(2:tcount+1),'r.'); %hold on,
set (EFH,'EraseMode','Xor','MarkerSize',6)
xlabel('t'); ylabel('Energy');
title(['Energy = ' num2str(Energy_t) ' ' 'Energyerr = ' num2str(Energy_err)])
drawnow

subplot(4,2,8)
CFH=plot([1:tcount]*dt, Current_t(2:tcount+1),'k-', [1:tcount]*dt, Current_t(2:tcount+1),'r.'); %hold on,
set (CFH,'EraseMode','Xor','MarkerSize',6)
xlabel('t'); ylabel('Current');
title(['Current =' num2str(Current_total)])
drawnow
% % plotfolder = ['./Results/SepPlot/time=' num2str(time) ...
% %     '/StericG=' StrSteric];
% % if ~exist(plotfolder, 'dir')
% %        mkdir(plotfolder);
% %    %    print("creating new folder", matfolder);
% %    end
% % saveas(h, [plotfolder '/HOT = ' StrHOT '.eps']);

end

if option ==2
    StrSteric = [num2str(StericG(1,1)) ' ' num2str(StericG(1,2)) ' ' ...
        num2str(StericG(2,1)) ' ' num2str(StericG(2,2))];
    StrHOT = [num2str(HOT(1,1)) ' ' num2str(HOT(1,2)) ' ' ...
        num2str(HOT(2,1)) ' ' num2str(HOT(2,2))];
    plotfolder = ['./Results/SepPlot/time=' num2str(time) ...
    '/StericG=' StrSteric];
   if ~exist(plotfolder, 'dir')
       mkdir(plotfolder);
   %    print("creating new folder", matfolder);
   end

    ConcP_h=figure(1);
ConcP_FigHandle = plot(x,Ion1,'k-',x,Ion1,'r.');
set (ConcP_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('C1(x,t)');
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str) 
drawnow
saveas(ConcP_h, [plotfolder '/HOT = ' StrHOT  ' C1.eps']);
    
ConcN_h=figure(2);
ConcN_FigHandle = plot(x,Ion2,'k-',x,Ion2,'r.');
set (ConcN_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('C2(x,t)'); 
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str)
drawnow
saveas(ConcN_h, [plotfolder '/HOT = ' StrHOT  ' C2.eps']);
    

Phi_h=figure(3);
PhiFigHandle = plot(x,Phi,'k-',x,Phi,'r.');
set (PhiFigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('\phi(x,t)');
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str)
drawnow
saveas(Phi_h, [plotfolder '/HOT = ' StrHOT ' Phi.eps']);
    
J_h=figure(4);
JFigHandle = plot(x,I,'k-',x,I,'r.');
set (PhiFigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('J(x,t)');
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str)
drawnow
saveas(J_h, [plotfolder '/HOT = ' StrHOT ' J.eps']);
end
end
