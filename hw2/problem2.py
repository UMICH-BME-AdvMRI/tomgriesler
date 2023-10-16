#%%
import h5py
import matplotlib.pyplot as plt
import numpy as np

#%%
from tools import sesignal, xrot, yrot, freeprecess, fsesignal

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

# %%
def fse_signal_5TR(T1, T2, TE):

    R90 = xrot(np.pi/2)
    R180 = yrot(np.pi)

    Ate2, Bte2 = freeprecess(TE/2, T1, T2, 0)

    M = np.array([[0], [0], [1]])

    # 90° excitation
    M = R90 @ M

    # 180° excitations and relaxation
    for _ in range(5): 
        M = Ate2 @ R180 @ (Ate2 @ M + Bte2) + Bte2

    return M[0] + 1j*M[1]

#%%
te_arr = np.linspace(0, 100, 1000)

for T1, T2 in [(1000, 50), (1000, 100), (2000, 50), (2000, 100)]:

    signal = [np.squeeze(fse_signal_5TR(T1, T2, te)) for te in te_arr]

    plt.plot(te_arr, -np.imag(signal), label=f'T1=%.0f, T2=%.0f' %(T1, T2))

plt.legend()
plt.xlabel('TE [s]')
plt.ylabel('signal')
plt.tight_layout()

# %%
def create_fse_image(t1map, t2map, m0map, esp, te_eff, etl, tr):

    k = int(te_eff/esp)

    result = np.zeros_like(t1map, dtype=complex)

    for ii in range(np.shape(t1map)[0]):
        for jj in range(np.shape(t1map)[1]):
            if t1map[ii, jj] == 0 or t2map[ii, jj] == 0:
                continue
            result[ii, jj] = m0map[ii, jj] * fsesignal(t1map[ii, jj], t2map[ii, jj], esp, tr, etl)[k-1]

    return np.abs(result)

# %%
tr = 3000
esp = 5

#%%
result = create_fse_image(t1map, t2map, m0map, esp, 80, 32, tr)
plt.imshow(np.abs(result), cmap='gray')
plt.axis('off')
plt.title('ETL=32, TE$_{eff}=80\,$ms')

# %%
result = create_fse_image(t1map, t2map, m0map, esp, 40, 32, tr)
plt.imshow(np.abs(result), cmap='gray')
plt.axis('off')
plt.title('ETL=32, TE$_{eff}=40\,$ms')

# %%
result = create_fse_image(t1map, t2map, m0map, esp, 120, 32, tr)
plt.imshow(np.abs(result), cmap='gray')
plt.axis('off')
plt.title('ETL=32, TE$_{eff}=120\,$ms')

# %%
results = np.zeros((t1map.shape[0], t1map.shape[1], 4))

for ii, etl in enumerate([16, 32, 64, 128]):
    results[..., ii] = create_fse_image(t1map, t2map, m0map, esp, 80, etl, tr)
    print(ii)
#%%
for ii, etl in enumerate([16, 32, 64, 128]):
    plt.subplot(2, 2, ii+1)
    plt.imshow(np.abs(results[..., ii]), cmap='gray')
    plt.axis('off')
    plt.title('ETL=%i, TE$_{eff}=80\,$ms' %etl)
plt.tight_layout()
# %%
