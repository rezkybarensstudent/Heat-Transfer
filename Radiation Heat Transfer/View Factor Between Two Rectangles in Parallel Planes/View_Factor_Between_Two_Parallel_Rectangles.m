function View_Factor_Between_Two_Parallel_Rectangles
Set1 = F1_2(-1, 1, -1, 1, -1, 1, -1, 1, 2)
Set2 = F1_2(-2, 0, -2, 0, 2, 0, 2, 0, 2)
function F12 = F1_2(x1a, x1b, y1a, y1b, x2a, x2b, y2a, y2b, dz)
A2 = abs(x1a-x1b)*abs(y1a-y1b);
F12 = dblquad(@OuterKernel, x1a, x1b, y1a, y1b, [], [], x2a, x2b, y2a, y2b, dz)/(A2*pi);
function f = OuterKernel(x1, y1, x2a, x2b, y2a, y2b, dz)
f = zeros(length(x1), 1);
for i = 1:length(x1)
f(i) = dblquad(@InnerKernel, x2a, x2b, y2a, y2b, [], [], dz, x1(i), y1);
end
function f = InnerKernel(x, y, dz, x2, y2)
L = length(x);
S = [x-x2*ones(1, L); (y-y2)*ones(1, L); dz*ones(1, L)];
n = repmat([0, 0, 1]', 1, L);
f = dot(n, S).^2./dot(S, S).^2;