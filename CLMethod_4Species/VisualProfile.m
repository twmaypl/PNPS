% function VisualProfile(Ion1,Ion2,Ion3,Ion4, Phi,I,x,StericG, HOT, time,option)
function VisualProfile(Ion1, Ion2, Ion3, Ion4, Phi, Current_x, I_ex, I_pnp, ...
    Current_t, Energy, x, StericG, HOT, dt, tcount, time, option, ...
    Energy_t, Energy_err, Current_total )

persistent ConcP_FigHandle ConcN_FigHandle PhiFigHandle TimeStrHandle

if option == 1 
    StrSteric = [num2str(StericG(1,1)) ', ' ...
        num2str(StericG(1,2)) ', ' num2str(StericG(3,3)) ', ' ...
        num2str(StericG(3,4))];
    StrHOT = [num2str(HOT(1,1)) ', ' num2str(HOT(1,2)) ', ' ...
        num2str(HOT(3,3)) ', ' num2str(HOT(3,4))];
    
    figure(1)

subplot(3,2,1)
Conc1_FigHandle = plot(x,Ion1,'r-',x,Ion1,'r.'); hold on,
Conc2_FigHandle = plot(x,Ion2,'blue-',x,Ion2,'blue.'); hold off
set (Conc1_FigHandle,'EraseMode','Xor','MarkerSize',6)
set (Conc2_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('Ion1(red) & Ion2(blue)');
time_str = ['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str) 
drawnow

subplot(3,2,2)
Conc3_FigHandle = plot(x,Ion3,'r-',x,Ion3,'r.'); hold on,
Conc4_FigHandle = plot(x,Ion4,'blue-',x,Ion4,'blue.'); hold off
set (Conc3_FigHandle,'EraseMode','Xor','MarkerSize',6)
set (Conc4_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('Ion3(red) & Ion4(blue)'); 
drawnow

subplot(3,2,3)
PhiFigHandle = plot(x,Phi,'k-',x,Phi,'r.');
set (PhiFigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('phi(x,t)');
drawnow

subplot(3,2,4)
JFigHandle1 = plot(x,I_pnp,'r-',x,I_pnp,'r.');
hold on,
JFigHandle2 = plot(x,I_ex, 'blue-',x,I_ex,'blue.');
hold on,
JFigHandle3 = plot(x,Current_x, 'green-',x,Current_x,'green.');
hold off
set (JFigHandle1,'EraseMode','Xor','MarkerSize',6)
set (JFigHandle2,'EraseMode','Xor','MarkerSize',6)
set (JFigHandle3,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('Current');
title('I_{pnp}-red; I_{ext}-blue; I-green')
drawnow

subplot(3,2,5)
EFH=plot([1:tcount]*dt, Energy(2:tcount+1),'k-', [1:tcount]*dt, Energy(2:tcount+1),'r.'); %hold on,
set (EFH,'EraseMode','Xor','MarkerSize',6)
xlabel('t'); ylabel('Energy');
title(['Energy = ' num2str(Energy_t) ' ' 'Energyerr = ' num2str(Energy_err)])
drawnow

subplot(3,2,6)
CFH=plot([1:tcount]*dt, Current_t(2:tcount+1),'k-', [1:tcount]*dt, Current_t(2:tcount+1),'r.'); %hold on,
set (CFH,'EraseMode','Xor','MarkerSize',6)
xlabel('t'); ylabel('Current');
title(['Current =' num2str(Current_total)])
drawnow

end

if option ==2

subplot(3,1,1);
    
    set(ConcP_FigHandle,'YData',Ion1);
    set(gca,'YLim',[0 max(Ion1)*1.1]);
    time_str=['time = ' num2str(time)];
    title(time_str);
    drawnow
subplot(3,1,2);    
    set(ConcN_FigHandle,'YData',Ion2);
    set(gca,'YLim',[0 max(Ion2)*1.1]);
    drawnow
subplot(3,1,3);
    set(PhiFigHandle,'YData',Phi);
    set(gca,'YLim',[min(Phi)-(max(Phi)-min(Phi))*0.1 max(Phi)+(max(Phi)-min(Phi))*0.1]);
    drawnow
end
end
