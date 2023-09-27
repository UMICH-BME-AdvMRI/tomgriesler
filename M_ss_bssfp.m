function M = M_ss_bssfp(alpha, T1, T2, TE, TR, df)
R = yrot(alpha);
[A_TE, B_TE] = freeprecess(TR-TE, T1, T2, df);
[A_TR, B_TR] = freeprecess(TE, T1, T2, df);
M = inv(eye(3)-A_TE*R*A_TR) * (A_TE*R*B_TR+B_TE);