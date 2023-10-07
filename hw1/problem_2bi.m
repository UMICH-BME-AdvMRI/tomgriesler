%%
clear;
T1 = 1000;
T2 = 100;
alpha = 10 * pi/180;
TE = 5;
TR = 10;
df = [-200:200];

%% One solution
signal = zeros(length(df), 1);

for jj = 1:length(df)
    M = M_ss_flash(alpha, T1, T2, TE, TR, df(jj));
    signal(jj) = M(1) + 1j*M(2);
end

figure
for ii=1:3
    plot(df, abs(signal));
    hold on
end
xlabel('Frequency [Hz]')
ylabel('Steady state signal')
ylim([0 0.07])

%% alternative solution
N = 500;
signal = zeros(N, 1);
M = [0 0 1]';
phase = 0;
phase_inc = pi;
TE = 5;
TR = 10;
[A_TE, B_TE] = freeprecess(TE, T1, T2, 0);
[A_TR, B_TR] = freeprecess(TR, T1, T2, 0);

for ii=1:N
    M = rot(alpha, phase) * M;
    M_TE = (A_TE * M + B_TE);
    signal(ii) = (M_TE(1) + 1j*M_TE(2)) * exp(-1j*phase);
    M = A_TR * M + B_TR;
    M = [0 0 M(3)]';
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




