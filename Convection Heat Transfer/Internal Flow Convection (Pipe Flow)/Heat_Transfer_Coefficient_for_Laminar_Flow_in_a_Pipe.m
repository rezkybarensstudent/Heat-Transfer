%  Heat Transfer Coefficient for Laminar Flow in a Pipe

% Write Function
function Heat_Transfer_Coefficient_for_Laminar_Flow_in_a_Pipe
% Wall Temperature
Tw = 40;
% The Boundary Condition for Constant Heat Flux
qw = 10;
% Laminar Flow
Re = 40;
Pr = 5;
% Pipe Radius
R = 0.01;
% Pipe Length
L = 0.5;
% The Comparison of Thermal Conductivity
k = 0.6;
Rt = 401;
zt = 50;
dxi = 1/(Rt-1);
xi = linspace(0, 1, Rt);
zeta = linspace(0, 1, zt);
% % Solve Initial-Boundary Value Problems for Systems of Parabolic and Elliptic Partial Differential Equations (PDEs) in One Space Variable and Time
solT = pdepe(1, @pdepde, @pdeic, @pdebcT, xi, zeta, [], Tw, qw, Re, Pr, R, L, k);
solF = pdepe(1, @pdepde, @pdeic, @pdebcF, xi, zeta, [], Tw, qw, Re, Pr, R, L, k);
NuT = zeros(zt,1);
NuF = NuT;
% Chart 1
figure(1)
for i = 1:zt
TmT = 4*trapz(xi, xi.*(1-xi.^2).*solT(i,:));
dThdxiT = (solT(i,Rt)-solT(i,Rt-1))/(dxi*(TmT-solT(i,Rt)));
NuT(i) = -2*dThdxiT;
TmF = 4*trapz(xi, xi.*(1-xi.^2).*solF(i,:));
dThdxiF = (solF(i,Rt)-solF(i,Rt-1))/(dxi*(TmF-solF(i,Rt)));
NuF(i) = -2*dThdxiF;
end
ThT = (solT(end,:)-ones(1,Rt)*Tw)/(TmT-Tw);
ThF = (solF(end,:)-ones(1,Rt)*solF(end,Rt))/(TmF-solF(end,Rt));
plot(xi, ThT, 'k-', xi, ThF, 'k--')
xlabel('\xi')
ylabel('\theta')
legend('Constant wall temperature', 'Constant wall heat flux')
% Chart 2
figure(2)
plot(zeta, NuT, 'k-', zeta, NuF, 'k--')
xlabel('\zeta')
ylabel('Nu')
ylim([0 6])
legend('Constant wall temperature', 'Constant wall heat flux')
% Write Function
function [c, f, s] = pdepde(xi, zeta, T, DTDxi, Tw, qw, Re, Pr, R, L, k)
c = Re*Pr*R/L*(1-xi^2);
f = DTDxi;
s = 0;
function T0 = pdeic(xi, Tw, qw, Re, Pr, R, L, k)
T0 = 20;
function [pl, ql, pr, qr] = pdebcT(xil, Tl, xir, Tr, zeta, Tw, qw, Re, Pr, R, L, k)
pl = 0;
ql = 1;
pr = Tr-Tw;
qr = 0;
function [pl, ql, pr, qr] = pdebcF(xil, Tl, xir, Tr, zeta, Tw, qw, Re, Pr, R, L, k)
pl = 0;
ql = 1;
pr = -qw;
qr = k/R;