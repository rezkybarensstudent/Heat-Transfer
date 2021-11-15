% Transient Heat Conduction in an Infinite Solid Cylinder With Convection

% Biot Number
Bi = 0.5;
% 
Nroot = 15;
CylinderRoots = inline('x.*besselj(1, x)-Bi*besselj(0, x)', 'x', 'Bi');
% Root of nonlinear function
r = FindZeros(CylinderRoots, Nroot, linspace(0, 50, 200), Bi);
% Generate linearly spaced vector
tau = linspace(0, 1.5, 20);
% 2-D and 3-D grids
[t, rt] = meshgrid(tau, r);
% Exponentialcollapse
Fn = exp(-t.*rt.^2);
cn = 2*besselj(1, r)./(r.*(besselj(0, r).^2+besselj(1, r).^2));
% 2-D and 3-D grids
ccn = meshgrid(cn, tau);
pro = ccn'.*Fn;
% Generate linearly spaced vector
rstar = linspace(0, 1, 20);
% 2-D and 3-D grids
[R, rx] = meshgrid(rstar, r);
% Bessel function of first kind
Jo = besselj(0, rx.*R);
the = Jo'*pro;
% 2-D and 3-D grids
[rr, tt] = meshgrid(rstar, tau);
% Mesh surface plot
mesh(rr, tt, the)
% \xi Label
xlabel('\xi')
% \tau label
ylabel('\tau')
% \theta label
zlabel('\theta')
% Camera line of sight
view(49.5, -34)