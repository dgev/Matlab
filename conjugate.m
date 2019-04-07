function F = conjugate(f,v,x0)
xk = x0;
g = gradient(f,v);
gk = subs(g,v,x0);
dk = -gk;
syms l
while gk ~= 0
alfa = [xk(1) + l.*dk(1), xk(2) + l.*dk(2)];
fk = matlabFunction(subs(f,v,alfa));
lk = double(fminsearch(fk,0));
xk1 = double(xk(1) + lk.*dk(1));
xk2 = double(xk(2) + lk.*dk(2));
xk = [xk1, xk2];
gk1 = double(subs(g,v,xk));
bk = double((mtimes(transpose(gk1), (gk1-gk)))./(mtimes(transpose(dk),(gk1-gk))));
dk1 = -gk1 + bk.*dk;
dk = double(dk1);
gk = double(gk1);
end
F = xk;
end