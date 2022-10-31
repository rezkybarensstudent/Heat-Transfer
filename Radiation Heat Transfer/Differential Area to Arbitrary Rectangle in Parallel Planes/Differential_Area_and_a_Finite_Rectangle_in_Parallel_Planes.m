function  Differential_Area_and_a_Finite_Rectangle_in_Parallel_Planes
Set1 = Fd1_2(-1, 0, -1, 0, 5)
Set2 = Fd1_2(-1, 1, -1, 1, 1)
N = 100; dz = linspace(0.1, 5, N);
Fd12 = zeros(N,1);
for i = 1:N
Fd12(i) = Fd1_2(0, 1, 0, 1, dz(i));
end
plot(dz, Fd12, 'k-')
xlabel('Separation distance of surfaces')
ylabel('View factor')
function F = Fd1_2(x_2a, x_2b, y_2a, y_2b, dz)
F = dblquad(@kernel2, x_2a, x_2b, y_2a, y_2b, [], [], dz)/pi;
function f = kernel2(x, y, dist)
L = length(x);
S = [x; repmat(y, 1, L); dist*ones(1, L)];
n = repmat([0, 0, 1]', 1, L);
f = dot(n, S).^2./dot(S, S).^2;