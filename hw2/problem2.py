#%%
import h5py
import matplotlib.pyplot as plt
import numpy as np

#%%
from tools import sesignal

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

    return np.abs(result)

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

#%% PD weighted
te = 15
tr = 4000
result = create_se_image(t1map, t2map, m0map, te, tr)
plot_image(result, f'PD weighted\nTE={te}ms, TR={tr}ms', m0map, 'Proton Density')

#%% T1 weighted
te = 15
tr = 500
result = create_se_image(t1map, t2map, m0map, te, tr)
plot_image(result, f'T1 weighted\nTE={te}ms, TR={tr}ms', t1map, 'T1')

#%%
te = 100
tr = 6000
result = create_se_image(t1map, t2map, m0map, te, tr)
plot_image(result, f'T2 weighted\nTE={te}ms, TR={tr}ms', t2map, 'T2')

#%%
images = {teff: {'image': create_se_image(t1map, t2map, m0map, teff, 3000)} for teff in np.arange(5, 165, 5)}
# %%
for key in images.keys():
    images.