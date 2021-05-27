%% Plots
figure(1) % Non linear Graph
subplot(611), plot(t,yNL(:,3)); grid;
title("Simulated Trajectories, "+ DR +"^{o} Deadrise ");
ylabel('\phi(^o)');
hold on
subplot(612), plot(t,yNL(:,4)); grid; 
ylabel('\psi(^o)');
hold on
subplot(613), plot(t,yNL(:,5)); grid; 
ylabel('u(ft/s)');
hold on
subplot(614), plot(t,yNL(:,6)); grid; 
ylabel('v(ft/s)');
hold on
subplot(615), plot(t,yNL(:,7)); grid; 
ylabel('p(^o/s)');
hold on
subplot(616), plot(t,yNL(:,8)); grid; 
ylabel('r(^o/s)');
hold on
xlabel('Time(s)')
% title("Simulated Trajectories, "+ DR +"^{o} Deadrise ");
% legend('baseline, No \phi','\tau=6^{o}, +\phi','\tau=3^{o}, +\phi')

figure(2) % x-y plot
title("Trajectories in x-y, "+ DR +"^{o} Deadrise")
plot(yNL(:,1),-yNL(:,2))
xlabel('ft')
ylabel('ft')
axis equal
hold on

figure(3)
plot(t,yNL(:,8)); grid;
title('Yaw rate')
xlabel('Time(s)');
ylabel('r');
hold on

figure(4)
subplot(411), plot(t,yNL(:,3)); grid; 
hold on
ylabel('\phi(^o)');
title('Nonlinear time domain simulation')
subplot(412), plot(t,yNL(:,7)); grid; 
ylabel('p(^o/s)');
hold on
subplot(413), plot(t,yNL(:,8)); grid; 
ylabel('r(^o/s)');
hold on
subplot(414), plot(t,yNL(:,10)); grid; 
ylabel('m_a velocity (ft/s)');
hold on
xlabel('Time(s)')

figure(5)
subplot(211)
plot(yNL(:,7),yNL(:,8)); hold on;
xlabel('p(^o/s)')
ylabel('r(^o/s)')
subplot(212)
plot(yNL(:,3),yNL(:,9)); hold on;
xlabel('\phi (^o)')
ylabel('d (ft)')

%legend
for i=1:5
    figure(i)
    legend('Open Loop, rudder & mass fixed','Single Actuator, AMA','Double Actuator, AMA + RRS','Single Actuator, RSS')
end
