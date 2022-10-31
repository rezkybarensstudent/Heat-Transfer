% Transient radiation heating of a plate in a furnace
function Transient_Radiation_Heating_of_a_Plate_in_a_Furnace
P1 = 1.67e-5; P2 = 8.8e-14; P3 = 6.3e-13;
Qguess = 1e5; Te = 1100; th = 600;
tend = 660; Two = 300; Tpo = 300;
Q = fzero(@QGen, Qguess, [], Te, th, Two, Tpo, tend, P1, P2, P3);
[t, T] = ode45(@RadTemp, [0, tend], [Two; Tpo], [], Q, Te, th, Two, Tpo, tend, P1, P2, P3);
plot(t, T(:,1), 'k-', t, T(:,2), 'k--')
z = axis;
hold on
plot([0, z(2)], [Te, Te], 'k.:', [th, th], [z(3), z(4)], 'k.:')
xlabel('Time (s)')
text(0.05*z(2), 0.85*z(4), ['Q = ' num2str(Q,6) ' W'])
ylabel('Temperature (K)')
legend('Wall temperature', 'Plate temperature', 'Location', 'NorthWest')
function dTdt = RadTemp(t, T, Q, Te, th, Two, Tpo, tend, P1, P2, P3)
dTdt = [P1*Q-P2*(T(1)^4-T(2)^4); -P3*(T(2)^4-T(1)^4)];
function PlateTempDev = QGen(Q, Te, th, Two, Tpo, tend, P1, P2, P3)
[t,T] = ode45(@RadTemp, [0, tend], [Two;Tpo], [], Q,Te, th,Two,Tpo, tend, P1, P2, P3);
PlateTempDev = Te-interp1(t, T(:,2), th, 'spline');