close all;
%load solvers.mat

%  plot(ode113(:,1),ode113(:,2),'r'); hold on; plot(ode113(:,1),ode113(:,3),'r'); hold on;
%  plot(ode23(:,1),ode23(:,2),'g'); hold on; plot(ode23(:,1),ode23(:,3),'g'); hold on;
%  plot(ode23(:,1),ode45(:,2),'b'); hold on; plot(ode23(:,1),ode45(:,3),'b'); hold on;
%  legend('ode113 x1','ode113 x2','ode23 x1','ode23 x2','ode45 x1','ode45 x2');
%  title('different solvers from Simulink for the same problem');

 plot(ode113(:,2),'r'); hold on; plot(ode113(:,3),'r'); hold on;
 plot(ode23(:,2),'g'); hold on; plot(ode23(:,3),'g'); hold on;