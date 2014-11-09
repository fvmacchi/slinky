function [x b] = solver(K, x, b, FIXED, FREE)
Kf = K(FREE, FREE);
Kfe = K(FREE, FIXED);

Bf = b([FREE]);
xe = x([FIXED]);

xf = gaussseidel(Kf, Bf - Kfe*xe);

x(FREE) = xf;

b = K*x;
end