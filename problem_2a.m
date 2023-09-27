%% Problem 2a
T1 = 1000;
T2 = 100;
TE = [2.5 5 10];
TR = [5 10 20];
alpha = 60*pi/180;

df = [-200:200];

signal = zeros(length(df), length(TE));

for ii = 1:length(TE)
    for jj = 1:length(df)
        M = M_ss_bssfp(alpha, T1, T2, TE(ii), TR(ii), df(jj));
        signal(jj, ii) = M(1) + 1j*M(2);
    end
end

figure
for ii=1:3
    plot(df, abs(signal(:, ii)));
    hold on
end
legend('TE=2.5ms, TR=5ms', 'TE=5ms, TR=10ms', 'TE=10ms, TR=20ms')
title('T1=1000ms, T2=100ms')
xlabel('Frequency [Hz]')
ylabel('Steady state signal')