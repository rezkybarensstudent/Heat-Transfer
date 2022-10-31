function  Natural_Convection_along_a_Heated_Plate
Pr = [.07 .7 7]; etaMax = [11, 8, 8]; xm = [10, 5, 5]; ym = [2, 0.8, 0.5];
guess = [0 0 0 0 0];
for k = 1:3
figure(k)
solinit = bvpinit(linspace(0, etaMax(k), 5), guess);
sol = bvp4c(@NatConv, @NatConvBC, solinit, [], Pr(k));
eta = linspace(0, etaMax(k), 300);
y = deval(sol, eta);
subplot(2, 1, 1)
plot(eta, y(1,:), '-.k', eta, y(2,:), '-k', eta, y(3,:), '--k')
legend('Stream function, f = y_1', 'Velocity, df/d\eta = y_2',...
'Shear, d^2f/d\eta^2 = y_3')
axis([0 xm(k) -0.2 ym(k)])
xlabel('\eta')
ylabel('y_1, y_2, y_3')
subplot(2, 1, 2)
plot(eta, y(4,:), '-k', eta, y(5,:), '--k')
legend('Temperature, T^* = y_4', 'Heat flux, dT^*/d\eta = y_5')
axis([0 xm(k) -1.2 1])
xlabel('\eta')
ylabel('y_4, y_5')
end
function ff = NatConv(eta, y, Pr)
ff = [y(2); y(3); -3*y(1)*y(3)+2*y(2)^2-y(4); y(5); -3*Pr*y(1)*y(5)];
function res = NatConvBC(ya, yb, Pr)
res = [ya(1); ya(2); ya(4)-1; yb(2); yb(4)];