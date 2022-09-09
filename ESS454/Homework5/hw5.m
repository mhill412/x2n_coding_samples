k = [5e-5 0 0; 0 2.05e-3 0; 0 0 2.05e-3];
r = [5e-4; 8.7e-4; 0];
q = -k*r;
a = [cos(30) sin(30) 0; -sin(30) cos(30) 0; 0 0 1];
kprime = a*k*transpose(a)