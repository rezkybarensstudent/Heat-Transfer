% One-Dimensional Transient Heat Transfer with Source

% Write Function
function One_Dimensional_Transient_Heat_Transfer_with_Source
% 
Bi = 0.1;
% 
Tr = 0.55;
% 
Sigma = 1;
% 
xi = linspace(0, 1, 41);
% 
tau = linspace(0, 1, 101);
% 
theta = pdepe(0, @pde1D, @pdeIC, @pdeBC, xi, tau, [], Bi, Tr, Sigma);
% 
z = 0:0.25:1;
% 
figure(1)
% 
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
%
axis([0 1 0.5 1])
%
xlabel('\tau')
%
ylabel('\theta')
%
figure(2)
%
[thmin imin] = min(theta(:,1));
%
plot(xi, theta(1,:), 'k-', 'LineWidth', 2)
%
hold on
%
plot(xi, theta(2,:) , 'k', xi, theta(imin,:), 'k:')
%
solinit = bvpinit(linspace(0, 1, 20), [1 1]);
%
sol = bvp4c(@barode, @barbc, solinit, [], Bi, Sigma);
%
x = linspace(0, 1, 100);
%
y = deval(sol, x);
%
plot(x, y(1,:), 'k--');
%
xlabel('\xi')
%
ylabel('\theta')
%
legend(['\tau = 0 (Initial condition)'], ['\tau = ' num2str(tau(2))],...
['\tau = ' num2str(tau(imin)) ' (Minimum at \xi = 0)'],...
'\tau > 2 (Steady state)', 'Location', 'SouthWest')
%
function dydx = barode(x, y, Bi, Sigma)
dydx = [y(2), -Sigma]';
%
function res = barbc(ya, yb, Bi, Sigma)
res = [ya(2)-Bi*ya(1), yb(1)-0.55]';
% 
function [c, f, s] = pde1D(x, t, u, DuDx, Bi, Tr, Sigma)
c = 1; f = DuDx; s = Sigma;
% 
function T0 = pdeIC(x, Bi, Tr, Sigma)
T0 = 1-0.45*x;
% 
function [pl, ql, pr, qr] = pdeBC(xl, ul, xr, ur, t, Bi, Tr, Sigma)
pr = ur-Tr; qr = 0;
pl = -Bi*ul; ql = 1;