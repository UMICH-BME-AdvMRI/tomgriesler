clear
load('Data_Assignment3_Problem2.mat')

%%
image_combined = zeros(length(kspaceData));

for ii=1:size(kspaceData, 3)
   image_temp = ifftshift(ifft2(kspaceData(:, :, ii)));
    image_combined = image_combined + coilmaps(:, :, ii) .* image_temp;
end

figure
imagesc(abs(image_combined))
axis equal
axis tight
xticks([])
yticks([])
colormap('gray')

%%
kspace_undersampled = kspaceData;

for ii=1:length(kspaceData)
    if mod(ii, 2)==0
        kspace_undersampled(ii, :, :) = 0;
    end
end

image_combined_undersampled = zeros(length(kspaceData));

for ii=1:size(kspaceData, 3)
    image_temp = ifftshift(ifft2(kspace_undersampled(:, :, ii)));
    image_combined_undersampled = image_combined_undersampled + coilmaps(:, :, ii) .* image_temp;
end

figure
imagesc(abs(image_combined_undersampled))
axis equal
axis tight
xticks([])
yticks([])
colormap('gray')

%%
kspace_undersampled_2 = kspaceData(1:2:end, :, :);