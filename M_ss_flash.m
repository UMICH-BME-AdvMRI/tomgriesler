function M = M_ss_flash(alpha, T1, T2, TE, TR, df)
R = yrot(alpha);
P = [0 0 0; 0 0 0; 0 0 1];
[A_TE, B_TE] = freeprecess(TR-TE, T1, T2, df);
[A_TR, B_TR] = freeprecess(TE, T1, T2, df);
M = inv(eye(3)-A_TE*R*A_TR*P) * (A_TE*R*B_TR+B_TE);
