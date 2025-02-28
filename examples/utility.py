import collections
import numpy as np

import scipy as sp
import scipy.sparse

DiscretizationInfo = collections.namedtuple("DiscretizationInfo",
    ("t0 tf N h "
     "t_coarse t_ctr int_order "
     "lg_t lg_c lg_t_vec lg_c_vec "
     "lg_t_mat lg_c_mat"))


def create_discretization_info(t0, tf, N, int_order):
    h = (tf - t0) / N
    lg_t, lg_c = np.polynomial.legendre.leggauss(int_order)
    t_coarse = np.linspace(t0, tf, N + 1)
    t_ctr = .5 * (t_coarse[1:] + t_coarse[:-1])

    t_ctr_reps = np.repeat(t_ctr[np.newaxis, :], int_order, axis=0)
    lg_t_reps = np.repeat(lg_t[:, np.newaxis], N, axis=1) * .5 * h
    lg_t_mat = t_ctr_reps + lg_t_reps
    lg_t_vec = lg_t_mat.flatten('F')
    lg_c_vec = np.tile(lg_c, N)
    lg_c_mat = np.reshape(lg_c_vec, lg_t_vec.shape, 'F')

    return DiscretizationInfo(
        t0=t0, tf=tf, N=N, h=h,
        t_coarse=t_coarse, t_ctr=t_ctr, int_order=int_order,
        lg_t=lg_t, lg_c=lg_c, lg_t_vec=lg_t_vec, lg_c_vec=lg_c_vec,
        lg_t_mat=lg_t_mat, lg_c_mat=lg_c_mat)


class LgConvPwcMat:
    def __init__(self, lg_x, di):
        self.mat_short = self.lg_cm_create(lg_x, di.N)
        #self.mat_short = self.mat[:di.int_order * di.N]

    def lg_cm_create(self, x, Ny):
        x = np.asarray(x)
        assert x.ndim == 2
        Nx = x.shape[1]
        M = x.shape[0]
        Nzi = Nx + Ny - 1
        Nz = M * Nzi
        lo_xi_stack = []

        coo_i = []
        coo_j = []
        coo_v = []

        n = np.arange(Nzi)[:, np.newaxis]
        m = np.arange(Ny)        

        print("Building sparse convolution matrix")
        for i in range(M):
            xi_stack = np.hstack([x[i], np.zeros(Nzi)]).flatten()
            for j in range(Ny):# Nzi formerly
                xi_shifted = xi_stack[j - m]
                col_indices = np.where(xi_shifted != 0.)[0]
                values = xi_shifted[col_indices]
                coo_i.append(np.repeat(j * M + i, col_indices.shape))
                coo_j.append(col_indices)
                coo_v.append(values)
            # spmat = sp.sparse.csr_matrix(xi_stack[n - m])
            # lo_xi_stack.append(xi_stack[n - m])
        # M_dense = sp.hstack(lo_xi_stack).reshape((Nz, Nx))

        

        M = sp.sparse.csr_matrix(
            (np.concatenate(coo_v),
                (np.concatenate(coo_i),
                 np.concatenate(coo_j))),
            (M * Ny, Nx) ## (Nz, Nx) formerly
        )
        print("... done.")
        return M # M_dense

    def conv(self, y):
        return self.mat_short.dot(y)


A = 0.1
omega0 = np.pi

t_left = -1.125
t_ctr = -1.
t_right = -0.93875
factor = .37

def a_left(t):
    return factor * (
        (.5 * t**2 - t_left * t) / (t_ctr - t_left)
        - (-.5 * t_left**2) / (t_ctr - t_left)
    )

def a_right(t):
    return factor * (
        (t_right * t - .5 * t**2) / (t_right - t_ctr)
        - (t_right * t_ctr - .5 * t_ctr**2) / (t_right - t_ctr)
    )

def a(t):
    val = np.zeros(t.shape)
    val[t <= t_left] = 0.
    val[(t > t_left) & (t < t_ctr)] = a_left(t[(t > t_left) & (t < t_ctr)])
    a_ctr = a_left(t_ctr)
    val[t == t_ctr] = a_ctr
    val[(t > t_ctr) & (t < t_right)] = \
        a_ctr + a_right(t[(t > t_ctr) & (t < t_right)])
    a_ctr_plus_right = a_ctr + a_right(t_right)
    val[t >= t_right] = a_ctr_plus_right 
    return val

def f(t):
    return .125 * np.cos(3. * t * np.pi - 0.25) * np.exp(1.1 *t)

def lg_int_mda(di):
    t = np.linspace(di.t0 - di.h, di.tf, di.N + 2)
    t_ctr = .5 * (t[1:] + t[:-1])
    t_ctr_mat = np.repeat(t_ctr[np.newaxis, :], di.int_order, axis=0)
    lg_t_reps = np.repeat(di.lg_t[:, np.newaxis], di.N + 1, axis=1) * .5 * di.h
    t_mat = t_ctr_mat + lg_t_reps
    a0c_stack = a(t_mat)
    int_mda_stack = -(a0c_stack[:, 1:] - a0c_stack[:, :-1])

    import matplotlib.pyplot as plt
    plt.plot(t_mat[0], a0c_stack[0])
    plt.show()



    assert np.allclose(t_mat[:, 1:], di.lg_t_mat)
    return int_mda_stack


def lg_jacobian(x, lg_cm, di, f_vec):
    y_vec = lg_cm.conv(x)
    return .5 * .5 * di.h * lg_cm.mat_short.transpose().dot(np.multiply(di.lg_c_vec, 2. * (y_vec - f_vec)))


def lg_objective(x, lg_cm, di, f_vec):
    y_vec = lg_cm.conv(x)
    return .5 * .5 * di.h * np.inner(di.lg_c_vec, (y_vec - f_vec)**2)
