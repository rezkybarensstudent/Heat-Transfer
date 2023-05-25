% Transient Heat Conduction in an Infinite Solid Cylinder with Convection
Bi = .5;
Nroot = 15;
CylinderRoots = inline('x.*besselj(1, x)-Bi*besselj(0, x)', 'x', 'Bi');
r = FindZeros(CylinderRoots, Nroot, linspace(0, 50, 200), Bi);
tau = linspace(0, 1.5, 20);
[t, rt] = meshgrid(tau, r);
Fn = exp(-t.*rt.^2);
cn = 2*besselj(1, r)./(r.*(besselj(0, r).^2+besselj(1, r).^2));
ccn = meshgrid(cn, tau);
pro = ccn'.*Fn;
rstar = linspace(0, 1, 20);
[R, rx] = meshgrid(rstar, r);
Jo = besselj(0, rx.*R);
the = Jo'*pro;
[rr, tt] = meshgrid(rstar, tau);
mesh(rr, tt, the')
xlabel('\xi')
ylabel('\tau')
zlabel('\theta')
view(49.5, -34)