#%%
import h5py
import matplotlib.pyplot as plt
import numpy as np

#%%
from epg import se_signal

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
            result[ii, jj] = se_signal(t1map[ii, jj], t2map[ii, jj], m0map[ii, jj], te, tr)

    return -np.imag(result)

#%%
te = 10
tr = 1000
plt.figure()

plt.subplot(1, 2, 1)
plt.imshow(create_se_image(t1map, t2map, m0map, te, tr), cmap='gray')
plt.title(f'PD weighted\nTE={te}ms, TR={tr}ms')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(m0map, cmap='gray')
plt.colorbar(orientation='horizontal', pad=0.02)
plt.title('Proton Density')
plt.axis('off')

plt.tight_layout()
plt.show()

#%%
te = 10
tr = 100
plt.figure()

plt.subplot(1, 2, 1)
plt.imshow(create_se_image(t1map, t2map, m0map, te, tr), cmap='gray')
plt.title(f'T1 weighted\nTE={te}ms, TR={tr}ms')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(t1map, cmap='gray')
plt.colorbar(orientation='horizontal', pad=0.02, label='T1 [ms]')
plt.title('T1')
plt.axis('off')

plt.tight_layout()
plt.show()

#%%
te = 100
tr = 1000
plt.figure()

plt.subplot(1, 2, 1)
plt.imshow(create_se_image(t1map, t2map, m0map, te, tr), cmap='gray')
plt.title(f'T2 weighted\nTE={te}ms, TR={tr}ms')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(t2map, cmap='gray')
plt.colorbar(orientation='horizontal', pad=0.02, label='T2 [ms]')
plt.title('T2')
plt.axis('off')

plt.tight_layout()
plt.show()