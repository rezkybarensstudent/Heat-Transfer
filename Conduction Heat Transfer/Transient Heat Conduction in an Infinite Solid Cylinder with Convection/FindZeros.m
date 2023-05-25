function Rt = FindZeros(FunName, Nroot, x, w)
f = feval(FunName, x, w);
indx = find(f(1:end-1).*f(2:end)<0);
L = length(indx);
if L<Nroot
Nroot = L;
end
Rt = zeros(Nroot, 1);
for k = 1:Nroot
Rt(k) = fzero(FunName, [x(indx(k)), x(indx(k)+1)], [], w);
end