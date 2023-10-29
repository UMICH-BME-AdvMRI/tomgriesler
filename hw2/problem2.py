#%%
import h5py
import matplotlib.pyplot as plt
import numpy as np

#%%
from tools import sesignal, fsesignal

#%%
brain_maps = h5py.File('brain_maps.mat')

t1map = np.array(brain_maps['T1map'])
t2map = np.array(brain_maps['T2map'])
m0map = np.array(brain_maps['M0map'])

#%%
def create_se_image(t1map, t2map, m0map, te, tr):

    result = np.zeros_like(t1map, dtype=complex)

    for ii in range(np.shape(t1map)[0]):
        for jj in range(np.shape(t1map)[1]):
            if t1map[ii, jj] == 0 or t2map[ii, jj] == 0:
                continue
            result[ii, jj] = m0map[ii, jj] * sesignal(t1map[ii, jj], t2map[ii, jj], te, tr)

    return result


def create_fse_images(t1map, t2map, m0map, esp, tr, etl):
    
    result = np.zeros((np.shape(t2map) + (etl,)), dtype=complex)

    for ii in range(np.shape(t1map)[0]):
        for jj in range(np.shape(t1map)[1]):
            if t1map[ii, jj] == 0 or t2map[ii, jj] == 0:
                continue
            result[ii, jj] = m0map[ii, jj] * fsesignal(t1map[ii, jj], t2map[ii, jj], esp, tr, etl)

    return result


#%%
def plot_image(result, title, compare, comparetitle):

    plt.figure()

    plt.subplot(1, 2, 1)
    plt.imshow(result, cmap='gray')
    plt.title(title)
    plt.axis('off')

    plt.subplot(1, 2, 2)
    plt.imshow(compare, cmap='gray')
    plt.colorbar(orientation='horizontal', pad=0.02)
    plt.title(comparetitle)
    plt.axis('off')

    plt.tight_layout()
    plt.show()

#%% PD weighted SE
te = 15
tr = 4000
result = np.abs(create_se_image(t1map, t2map, m0map, te, tr))
plot_image(result, f'PD weighted\nTE={te}ms, TR={tr}ms', m0map, 'Proton Density')

#%% T1 weighted SE
te = 15
tr = 500
result = np.abs(create_se_image(t1map, t2map, m0map, te, tr))
plot_image(result, f'T1 weighted\nTE={te}ms, TR={tr}ms', t1map, 'T1')

#%% T2 weighted SE
te = 100
tr = 6000
result = np.abs(create_se_image(t1map, t2map, m0map, te, tr))
plot_image(result, f'T2 weighted\nTE={te}ms, TR={tr}ms', t2map, 'T2')

#%% FSE
tr = 3000
esp = 5
etl = 128
te_eff = 80

result = create_fse_images(t1map, t2map, m0map, esp, tr, etl)

# %%
fts = np.zeros_like(result, dtype=complex)

for ii in range(etl):
    fts[..., ii] = np.fft.fftshift(np.fft.fft2(result[..., ii]))

# %%
final_kpace = np.zeros_like(t1map, dtype=complex)

index = int((etl/2-te_eff/esp)%etl)

for ii in range(etl):

    print(index)

    k1 = int(index*256/etl)
    k2 = int((index+1)*256/etl)

    final_kpace[k1:k2] = fts[k1:k2, :, ii]

    index = int((index+1)%etl)


plt.imshow(np.abs(np.fft.ifft2(final_kpace)))
# %%
