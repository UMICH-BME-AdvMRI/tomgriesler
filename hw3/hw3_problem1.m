%%
clear

%%
load('Data_Assignment3_Problem1.mat')

%%
kspace_partialFourier = zeros(size(kspaceData_SingleCoil));
kspace_partialFourier(1:round(5/8*length(kspaceData_SingleCoil)), :) = kspaceData_SingleCoil(1:round(5/8*length(kspaceData_SingleCoil)), :);

figure

tiledlayout(1, 2)

nexttile
imagesc(log(abs(kspaceData_SingleCoil)))
axis equal
axis tight
xticks([])
yticks([])
title("Fully sampled")

nexttile
imagesc(log(abs(kspace_partialFourier)))
axis equal
axis tight
xticks([])
yticks([])
title("Partial Fourier")

sgtitle("k space", 'FontWeight','bold')
colormap("gray")


%%
image_fullySampled = ifftshift(ifft2(kspaceData_SingleCoil));
image_partialFourier = ifftshift(ifft2(kspace_partialFourier));

figure

tiledlayout(2, 3, 'TileSpacing','tight', 'Padding','tight')

nexttile
imagesc(abs(image_fullySampled), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('Fully sampled')
ylabel('Magnitude', 'FontWeight','bold')

nexttile
imagesc(abs(image_partialFourier), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('Partial Fourier')

nexttile
imagesc(10*abs(image_fullySampled-image_partialFourier), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('Difference x10')

nexttile
imagesc(mod(angle(image_fullySampled), pi))
axis equal
axis tight
xticks([])
yticks([])
ylabel('Phase', 'FontWeight','bold')
title('Fully sampled')

nexttile
imagesc(mod(angle(image_partialFourier), pi))
axis equal
axis tight
xticks([])
yticks([])
title('Partial Fourier')

nexttile
imagesc(mod(angle(image_fullySampled-image_partialFourier), pi))
axis equal
axis tight
xticks([])
yticks([])
title('Difference')

colormap("gray")

%%
kspace_phaseEstimate = zeros(size(kspaceData_SingleCoil));
k1 = round(3/8*length(kspaceData_SingleCoil));
k2 = round(5/8*length(kspaceData_SingleCoil));

kspace_phaseEstimate(k1:k2, :) = kspaceData_SingleCoil(k1:k2, :);

figure
tiledlayout(1, 2)

nexttile
imagesc(log(abs(kspace_phaseEstimate)))
axis equal
axis tight
xticks([])
yticks([])
title('w/o hanning filter')

filter = hanning(size(kspace_phaseEstimate, 1)) * hanning(size(kspace_phaseEstimate, 2))';
kspace_phaseEstimate = filter .* kspace_phaseEstimate;

nexttile
imagesc(log(abs(kspace_phaseEstimate)))
axis equal
axis tight
xticks([])
yticks([])
title('w/ hanning filter')

image_phaseEstimate = ifftshift(ifft2(kspace_phaseEstimate));

%%
figure
tiledlayout(2, 2, 'TileSpacing','tight', 'Padding','tight')

nexttile
imagesc(abs(image_partialFourier))
axis equal
axis tight
xticks([])
yticks([])
title('Magnitude')
ylabel('Partial Fourier', 'FontWeight','bold')

nexttile
imagesc(mod(angle(image_partialFourier), pi))
axis equal
axis tight
xticks([])
yticks([])
title('Phase')

nexttile
imagesc(abs(image_phaseEstimate))
axis equal
axis tight
xticks([])
yticks([])
ylabel('Low-Resolution Estimate', 'FontWeight','bold')

nexttile
imagesc(mod(angle(image_phaseEstimate), pi))
axis equal
axis tight
xticks([])
yticks([])

colormap('gray')

%% Iterative POCS
N = 10;

magnitude_temp = abs(image_partialFourier);
phase_temp = angle(image_phaseEstimate);
image_temp = magnitude_temp .* exp(1j * phase_temp);

for ii=1:N

    kspace_temp = fft2(image_temp);

    k = round(5/8*length(kspaceData_SingleCoil));

    kspace_temp(1:k) = kspaceData_SingleCoil(1:k);

    image_temp = ifft2(kspace_temp);

end

figure
tiledlayout(2, 3, 'TileSpacing','tight', 'Padding','tight')

nexttile
imagesc(abs(image_fullySampled), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('Fully sampled')
ylabel('Magnitude', 'FontWeight','bold')

nexttile
imagesc(abs(image_temp), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('POCS (N=10)')

nexttile
imagesc(10*abs(image_fullySampled-image_temp), [0 7e-3])
axis equal
axis tight
xticks([])
yticks([])
title('Difference x10')

nexttile
imagesc(mod(angle(image_fullySampled), pi))
axis equal
axis tight
xticks([])
yticks([])
ylabel('Phase', 'FontWeight','bold')
title('Fully sampled')

nexttile
imagesc(mod(angle(image_temp), pi))
axis equal
axis tight
xticks([])
yticks([])
title('POCS (N=10)')

nexttile
imagesc(mod(angle(image_fullySampled-image_temp), pi))
axis equal
axis tight
xticks([])
yticks([])
title('Difference')

colormap("gray")