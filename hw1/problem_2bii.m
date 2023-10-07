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
for deph=[2*pi 4*pi 8*pi 16*pi]   
    beta = linspace(deph/N_iso, deph, N_iso);
    M = zeros(3, N_iso);
    M(3, :) = 1/N_iso;
    signal = zeros(N_iter, 1);
    phase = 0;
    phase_inc = pi;
    [A_TE, B_TE] = freeprecess(TE, T1, T2, 0);
    [A_TR, B_TR] = freeprecess(TR, T1, T2, 0);
    
    for ii=1:N_iter
    
        % excitation
        for jj=1:N_iso
            M(:, jj) = rot(alpha, phase) * M(:, jj);
        end
    
        % signal
        M_TE = M;
        for jj=1:N_iso
            M_TE(:, jj) = zrot(beta(jj)) * (A_TE * M(:, jj) + B_TE/N_iso);
        end
        M_TE_summed = sum(M_TE, 2);
        signal(ii) = (M_TE_summed(1) + 1j*M_TE_summed(2)) * exp(-1j * phase);
    
        % update magnetization
        for jj=1:N_iso
            M(:, jj) = zrot(beta(jj)) * (A_TR * M(:, jj) + B_TR/N_iso);
        end
        
        phase = mod(phase + phase_inc, 2*pi);
    end
    
    figure
    plot(real(signal))
    hold on
    plot(imag(signal))
    hold on
    plot(abs(signal))
    legend('real', 'imag', 'abs')
    xlabel('Iteration')
    ylabel('Signal intensity')
    title(sprintf('%i cycles dephasing', deph/2/pi))
end
