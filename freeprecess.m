function [A, B] = freeprecess(t, T1, T2, df)
alpha = 2*pi*df*t/1000;
E1 = exp(-t/T1);
E2 = exp(-t/T2);
A = [E2 0 0; 0 E2 0; 0 0 E1] * zrot(alpha);
B = [0 0 1-E1]';