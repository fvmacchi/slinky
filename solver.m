function [x b] = solver(K, x, b, FIXED, FREE)
Kf = k(FREE, FREE);
Kfe = k(FREE, FIXED);

Bf = b([FREE]);
xe = x([FIXED]);

xf = gaussseidel(Kf, Bf - Kfe*xe);

x(FREE) = xf;

b = k*x;
end