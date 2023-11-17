clear
load('Data_Assignment3_Problem2.mat')

%%
image_combined = zeros(length(kspaceData));

for ii=1:size(kspaceData, 3)
   image_temp = ifftshift(ifft2(kspaceData(:, :, ii)));
    image_combined = image_combined + image_temp ./ coilmaps(:, :, ii) ; % coilmaps(:, :, ii) .* 
end

figure
imagesc(abs(image_combined), [0 0.05])
axis equal
axis tight
xticks([])
yticks([])
colormap('gray')
colorbar

%%
kspace_undersampled = kspaceData;

for ii=1:200
    if mod(ii, 2)==0
        kspace_undersampled(ii, :, :) = 0;
    end
end

image_combined_undersampled = zeros(200);

for ii=1:size(8)
    image_temp = ifftshift(ifft2(kspace_undersampled(:, :, ii)));
    image_combined_undersampled = image_combined_undersampled + image_temp; % ./ coilmaps(:, :, ii);
end

figure
imagesc(abs(image_combined_undersampled))
axis equal
axis tight
xticks([])
yticks([])
colormap('gray')
colorbar

%%
kspace_undersampled_2 = kspaceData(1:2:end, :, :);

images_undersampled = zeros(100, 200);

for ii=1:8
    images_undersampled(:, :, ii) = ifftshift(ifft2(kspace_undersampled_2(:, :, ii)));
end

%%
figure
imagesc(abs(images_undersampled(:, :, 2)))
axis equal
axis tight
xticks([])
yticks([])
colormap('gray')

%%
image_sense = zeros(200);

for ii=1:100
    for jj=1:200

        C = zeros(8, 2);
        I = zeros(8, 1);

        for kk=1:8
            C(kk, 1) = coilmaps(ii, jj, kk);
            C(kk, 2) = coilmaps(ii+100, jj, kk);
            I(kk) = images_undersampled(ii, jj, kk);
        end
        
        rho = inv(C' * C) * C' * I;
        
        image_sense(ii, jj) = rho(1);
        image_sense(ii+100, jj) = rho(2);
    end
end

%%
figure
imagesc(abs(image_sense))
axis equal
axis tight
colormap('gray')
xticks([])
yticks([])

