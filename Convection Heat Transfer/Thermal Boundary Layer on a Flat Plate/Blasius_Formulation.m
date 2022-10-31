% Heat Transfer from a Flat Plate: Blasius formulation
function Blasius_Formulation
Pr = [0.07, 0.7, 7.0]; etaMax = [15, 8, 8]; xm = [15, 5, 5];
for k=1:3
figure(k)
solinit = bvpinit(linspace(0, etaMax(k), 8), [0, 0, 0, 0, 0]);
sol = bvp4c(@BlasiusT, @BlasiusTbc, solinit, [], Pr(k));
eta = linspace(0, etaMax(k));
y = deval(sol, eta);
subplot(2, 1, 1)
plot(eta, y(1,:), '-.k', eta, y(2,:), '-k', eta, y(3,:), '--k')
xlabel('\eta')
ylabel('y_1, y_2, y_3')
legend('Stream function, f = y_1', 'Velocity, df/d\eta = y_2', 'Shear, d^2f/d\eta^2= y_3')
axis([0 xm(k) 0 2])
subplot(2,1,2)
plot(eta, y(4,:), '-k', eta, y(5,:), '--k')
axis([0 xm(k) 0 2])
legend('Temperature, T^* = y_4', 'Heat flux, dT^*/d\eta = y_5')
xlabel('\eta')
ylabel('y_4, y_5')
end
function F = BlasiusT(eta, y, Pr)
F = [y(2); y(3); -0.5*y(1)*y(3); y(5); -Pr*0.5*y(1)*y(5)];
function res = BlasiusTbc(ya, yb, Pr)
res = [ya(1); ya(2); ya(4); yb(2)-1; yb(4)-1];