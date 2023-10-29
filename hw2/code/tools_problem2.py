import numpy as np

def xrot(phi):

    return np.array([[1, 0, 0],
                     [0, np.cos(phi), -np.sin(phi)],
                     [0, np.sin(phi), np.cos(phi)]])


def yrot(phi): 

    return np.array([[np.cos(phi), 0, np.sin(phi)],
                     [0, 1, 0],
                     [-np.sin(phi), 0, np.cos(phi)]])


def zrot(phi):
    
    return np.array([[np.cos(phi), -np.sin(phi), 0],
                     [np.sin(phi), np.cos(phi), 0], 
                     [0, 0, 1]])


def freeprecess(t, T1, T2, df): 

    phi = 2*np.pi*df*t/1000

    E1 = np.exp(-t/T1)
    E2 = np.exp(-t/T2)

    A = np.array([[E2, 0, 0], [0, E2, 0], [0, 0, E1]]) * zrot(phi)
    
    B = np.array([[0], [0], [1-E1]])
    
    return A, B


def sesignal(T1, T2, TE, TR):
    """
    Steady state Spin Echo signal.
    """

    R90 = xrot(np.pi/2)
    R180 = yrot(np.pi)

    Atr, Btr = freeprecess(TR-TE, T1, T2, 0)
    Ate2, Bte2 = freeprecess(TE/2, T1, T2, 0)

    Atr = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]]) @ Atr

    Mss = np.linalg.inv(np.eye(3) + Ate2 @ R180 @ Ate2 @ R90 @ Atr) @ (Bte2 + Ate2 @ R180 @ (Bte2 + Ate2 @ R90 @ Btr))

    return Mss[0] + 1j*Mss[1]


def fsesignal(T1, T2, ESP, TR, ETL):
    """
    Steady state Fast Spin Echo signal. 
    """

    R90 = xrot(np.pi/2) 
    R180 = yrot(np.pi)

    Atr, Btr = freeprecess(TR-ETL*ESP, T1, T2, 0) 
    Aesp2, Besp2 = freeprecess(ESP/2, T1, T2, 0) 

    Atr = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]]) @ Atr  

    A = np.eye(3)
    B = np.array([[0], [0], [0]])

    for k in range(ETL):
        A = Aesp2 @ R180 @ Aesp2 @ A 
        B = Besp2 + Aesp2 @ R180 @ (Besp2 + Aesp2 @ B)

    A = R90 @ Atr @ A
    B = R90 @ (Btr + Atr @ B)

    M = np.linalg.inv(np.eye(3) - A) @ B 

    Msig = np.zeros(ETL, dtype=complex)
    for k in range(ETL):
        M = Aesp2 @ R180 @ Aesp2 @ M + Besp2 + Aesp2 @ R180 @ Besp2
        Msig[k] = M[0] + 1j * M[1]

    return Msig


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