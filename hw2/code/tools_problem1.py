import numpy as np


def q(alpha, phi=np.pi/2):
    """
    Basically taken from Rad229.
    """
    mat = np.array([
        [(np.cos(alpha / 2.)) ** 2., np.exp(2. * 1j * phi) * (np.sin(alpha / 2.)) ** 2.,
         -1j * np.exp(1j * phi) * np.sin(alpha)],
        [np.exp(-2. * 1j * phi) * (np.sin(alpha / 2.)) ** 2., (np.cos(alpha / 2.)) ** 2.,
         1j * np.exp(-1j * phi) * np.sin(alpha)],
        [-1j / 2. * np.exp(-1j * phi) * np.sin(alpha), 1j / 2. * np.exp(1j * phi) * np.sin(alpha), np.cos(alpha)]
    ])

    return mat


def r(t1, t2, t):
    e1 = np.exp(-t / t1)
    e2 = np.exp(-t / t2)

    mat = np.array([
        [e2, 0, 0],
        [0, e2, 0],
        [0, 0, e1]
    ])

    return mat


def b(t1, t):

    return 1-np.exp(-t/t1)


def epg_grad(omega):

    omega = np.hstack((omega, np.zeros((3, 1))))

    omega[0, 1:] = omega[0, :-1]
    omega[1, :-1] = omega[1, 1:]
    omega[1, -1] = 0
    omega[0, 0] = np.conj(omega[1, 0])

    return omega


def se_signal(t1: float, t2: float, m0: float, te: float, tr: float):

    omega = np.vstack((0, 0, m0*b(t1, tr-te)))

    r_te_2 = r(t1, t2, te/2)
    b_te_2 = m0 * b(t1, te/2)

    omega = epg_grad(r_te_2 @ q(np.deg2rad(180)) @ epg_grad(r_te_2 @ q(np.deg2rad(90), 0) @ omega + b_te_2) + b_te_2)

    return omega[0, 0]



def fse_signal(alpha: float, t1: float, t2: float, tr: float, n_echoes: int):

    omega = np.vstack((0, 0, 1))

    # precompute some matrices to save computation time
    signal = np.empty(n_echoes, dtype=complex)
    r_tr2 = r(t1, t2, tr/2)
    b_tr2 = b(t1, tr/2)
    q_alpha = q(np.deg2rad(alpha))

    # 90° excitation and relaxation
    omega = r_tr2 @ q(np.deg2rad(90), 0) @ omega + b_tr2

    # second excitation and relaxation
    omega = r_tr2 @ epg_grad(q(np.deg2rad(alpha/2+90)) @ epg_grad(omega)) + b_tr2
    signal[0] = omega[0, 0]

    # all other excitations
    for i in range(1, n_echoes):
        omega = r_tr2 @ epg_grad(q_alpha @ r_tr2 @ epg_grad(omega)) + b_tr2
        signal[i] = omega[0, 0]

    return signal


