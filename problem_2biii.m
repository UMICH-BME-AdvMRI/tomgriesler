%%
clear;
T1 = 1000;
T2 = 100;
alpha = 10 * pi/180;
TE = 5;
TR = 10;
N_iso = 200;
N_iter = 500;

%%
deph = 2*pi;
M_ss = zeros(length(deph), 1);
count = 1;

for phase_inc=0:5:180   

    beta = linspace(deph/N_iso, deph, N_iso);
    M = zeros(3, N_iso);
    M(3, :) = 1/N_iso;
    phase = 0;
    [A_TE, B_TE] = freeprecess(TE, T1, T2, 0);
    [A_TR, B_TR] = freeprecess(TR, T1, T2, 0);
    
    for ii=1:N_iter
    
        % excitation
        for jj=1:N_iso
            M(:, jj) = rot(alpha, phase) * M(:, jj);
        end
  
        % update magnetization
        for jj=1:N_iso
            M(:, jj) = zrot(beta(jj)) * (A_TR * M(:, jj) + B_TR/N_iso);
        end

        phase = mod(phase + phase_inc*pi/180, 2*pi);

    end
    
    M_sum = sum(M, 2);
    M_ss(count) = abs(M_sum(1) + 1j*M_sum(2));
    count = count + 1;
 
end

figure
plot([0:5:180], M_ss)
xlabel('RF phase')
ylabel('Steady state amplitude')