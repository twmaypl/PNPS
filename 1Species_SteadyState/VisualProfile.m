function VisualProfile(Ion1, Phi, Current_x, x, StericG, HOT, dt, tcount, time)
persistent ConcP_FigHandle ConcN_FigHandle PhiFigHandle TimeStrHandle


    StrSteric = num2str(StericG);
    StrHOT = num2str(HOT);
figure(5)

subplot(3,1,1)
ConcP_FigHandle = plot(x,Ion1,'r-',x,Ion1,'k.');
set (ConcP_FigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('p(red) & n(blue)');
time_str=['time = ' num2str(time) ' StericG = ' StrSteric ' HOT = ' StrHOT];
title(time_str) 
drawnow

subplot(3,1,2)
PhiFigHandle = plot(x,Phi,'k-',x,Phi,'r.');
set (PhiFigHandle,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('phi(x,t)');
drawnow

subplot(3,1,3)
JFigHandle1 = plot(x,Current_x, 'green-',x,Current_x,'green.');
set (JFigHandle1,'EraseMode','Xor','MarkerSize',6)
xlabel('x'); ylabel('Current');
drawnow

