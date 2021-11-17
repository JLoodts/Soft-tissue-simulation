close all;
if (length(tout)~=length(x1sim))
    tout = tout(3:length(x1sim)+2);
end;
plot(tout,x1sim,'r');
hold on;
plot(tout,x2sim,'g');
legend('x1','x2');
title('two masses and spring dampers');
xlabel('t (s)');
ylabel('position (m)');