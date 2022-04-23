% One-Dimensional Transient Heat Transfer with Source

% Write Function
function One_Dimensional_Transient_Heat_Transfer_with_Source
% Dimensional Source Strength
Bi = 0.1;
Tr = 0.55;
% Sigma
Sigma = 1;
% Generate Linearly Spaced Vector
xi = linspace(0, 1, 41);
tau = linspace(0, 1, 101);
% Solve Initial-Boundary Value Problems for Systems of Parabolic and Elliptic Partial Differential Equations (PDEs) in One Space Variable and Time
% @pde1D (The Partial Differential Equation)
% @pdeIC (The Initial Conditions)
% @pdeBC (The Boundary Conditions)
theta = pdepe(0, @pde1D, @pdeIC, @pdeBC, xi, tau, [], Bi, Tr, Sigma);
% The Temperature as a Function of Time at Five Locations in The Solid
z = 0:0.25:1;
% Chart 1
figure(1)
for k = 1:length(z)
    kk = find(xi == z(k));
    plot(tau, theta(:, kk), 'k-')
    hold on
    if k == 1
        text(0.5, 1.02*theta(end, kk), '\xi = 0.0 and 0.25')
    elseif k > 2
        text(0.5, theta(end, kk)+.02, ['\xi = ' num2str(xi(kk))])
    end
end
axis([0 1 0.5 1])
xlabel('\tau')
ylabel('\theta')
% Chart 2
figure(2)
[thmin imin] = min(theta(:,1));
plot(xi, theta(1,:), 'k-', 'LineWidth', 2)
hold on
plot(xi, theta(2,:) , 'k', xi, theta(imin,:), 'k:')
% Form Initial Guess for Boundary Value Problem Solver
solinit = bvpinit(linspace(0, 1, 20), [1 1]);
% Solve Boundary Value Problem — Fourth-Order Method
sol = bvp4c(@barode, @barbc, solinit, [], Bi, Sigma);
% Generate Linearly Spaced Vector
x = linspace(0, 1, 100);
% Evaluate differential equation solution structure
y = deval(sol, x);
plot(x, y(1,:), 'k--');
xlabel('\xi')
ylabel('\theta')
legend(['\tau = 0 (Initial condition)'], ['\tau = ' num2str(tau(2))],...
['\tau = ' num2str(tau(imin)) ' (Minimum at \xi = 0)'],...
'\tau > 2 (Steady state)', 'Location', 'SouthWest')
% Make Function
function dydx = barode(x, y, Bi, Sigma)
dydx = [y(2), -Sigma]';
function res = barbc(ya, yb, Bi, Sigma)
res = [ya(2)-Bi*ya(1), yb(1)-0.55]';
function [c, f, s] = pde1D(x, t, u, DuDx, Bi, Tr, Sigma)
c = 1;
f = DuDx;
s = Sigma;
function T0 = pdeIC(x, Bi, Tr, Sigma)
T0 = 1-0.45*x;
function [pl, ql, pr, qr] = pdeBC(xl, ul, xr, ur, t, Bi, Tr, Sigma)
pr = ur-Tr;
qr = 0;
pl = -Bi*ul;
ql = 1;