close all;
clear all;

load forceDeformation.txt
t = forceDeformation(:,1);
x = forceDeformation(:,2);
v = forceDeformation(:,3);
a = forceDeformation(:,4);
x_e = forceDeformation(:,5);
LF_x = forceDeformation(:,6);
LF_v = forceDeformation(:,7);
LF_a = forceDeformation(:,8);
% v_e = forceDeformation(:,6);
% a_e = forceDeformation(:,7);
clear forceDeformation;

% plot(f,dx,'r-');
% hold on;
% plot(f_e,x_e,'g-');
% plot(f-f_e,x-x_e,'r-');
% title('force deformation');

% figure;
% plot(t,a,'r-');
% hold on;
% plot(t,a_e,'g-');
% relErr_a = 100*(a-a_e)/max(abs(a_e));
% plot(t,relErr_a,'g-');
% title('force with time');

figure;
subplot(2,1,1);
plot(t,x,'r-');
hold on;
plot(t,LF_x,'b-');
hold on;
plot(t,x_e,'g-');
legend('simple scheme','leap frog scheme','ode45');
ylabel('deviation (m)');
title('simulating a spring damper system');
subplot(2,1,2);
relErr_x = 100*(x-x_e)./max(abs(x_e));
plot(t,relErr_x,'r-'); hold on;
LF_relErr_x = 100*(LF_x-x_e)./max(abs(x_e));
plot(t,LF_relErr_x,'b-');
legend('simple scheme','leap frog scheme');
xlabel('time (s)');
ylabel('relative error (%)');


% figure;
% plot(t,v,'r-');
% hold on;
% plot(t,v_e,'g-');
% relErr_v = 100*(v-v_e)/max(abs(v_e));
% plot(t,relErr_v,'g-');
% title('relative velocity with time');