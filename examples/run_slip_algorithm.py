import copy
import numpy as np
import time
import matplotlib.pyplot as plt

from utility import *
import trs4slip


def eval_tv(x):
    return np.sum(np.abs(x[1:] - x[:-1]))

def slip(eval_f, eval_jac, x0, lo_bangs, alpha, h, Delta0, sigma, maxiter, use_astar=False):
    assert x0.ndim == 1
    assert lo_bangs.ndim == 1
    N, M = x0.shape[0], lo_bangs.shape[0]

    xn = copy.deepcopy(x0)
    vert_costs_buffer = np.empty(N*M*(Delta0 + 1) + 2)
    vert_layer_buffer = np.empty(N*M*(Delta0 + 1) + 2, dtype=np.int32)
    vert_value_buffer = np.empty(N*M*(Delta0 + 1) + 2, dtype=np.int32)
    vert_prev_buffer = np.empty(N*M*(Delta0 + 1) + 2, dtype=np.int32)
    vert_remcap_buffer = np.empty(N*M*(Delta0 + 1) + 2, dtype=np.int32)
    xnk = np.empty((N,), dtype=np.int32)
    pred_positive = True
    Delta = 0
    pred_Delta0 = 0

    print("Iter         obj   pred(Delta0)   Delta (acc)")
    for n in range(maxiter):
        fn = eval_f(xn)
        gn = eval_jac(xn)
        tvn = eval_tv(xn)

        print("%4u   %.3e      %.3e    %4u" % (n, fn + alpha * tvn, pred_Delta0, Delta))
        Delta, k = Delta0, 0
        accept = False
        while Delta >= 1 and not accept and pred_positive:
            if use_astar:
                trs4slip.run(
                    xnk, gn / alpha, xn, lo_bangs, Delta,
                    vert_costs_buffer,
                    vert_layer_buffer,
                    vert_value_buffer,
                    vert_prev_buffer,
                    vert_remcap_buffer,
                    False,
                    0,
                    0
                )
            else:
                trs4slip.run_top(
                    xnk, gn / alpha, xn, lo_bangs,
                    np.ones((N - 1,)), Delta,
                    False,
                    0.,
                    0.,
                    1.,
                    1.
                )
            fnk = eval_f(xnk)
            tvnk = eval_tv(xnk)

            ared = fn - fnk + alpha * tvn - alpha * tvnk
            pred = gn.dot(xn - xnk) + alpha * tvn - alpha * tvnk

            if Delta == Delta0:
                pred_Delta0 = pred

            pred_positive = pred > 0. 
            accept = ared >= sigma * pred and pred_positive

            if accept:
                xn[:] = xnk[:]
            elif pred_positive:
                Delta = Delta / 2

            k = k + 1

        if Delta < 1:
            print("Trust region contracted. Solution may be close to stationarity.")
            break
        if not pred_positive:
            print("Predicted reduction is nonpositive. Solution may be close to stationarity.")
            break
    if n == maxiter - 1:
        print("Iteration limit (%d) reached. Solution may be instationary." % (maxiter))

    return xn


def lg_objective_var(x, lg_cm, di, f_vec):
    return lg_objective(x, lg_cm, di, f_vec)

def lg_jacobian_var(x, lg_cm, di, f_vec):
    return lg_jacobian(x, lg_cm, di, f_vec)


N = 2**10
di = create_discretization_info(-1., 1., N, 5)

lo_bangs = np.array([-1, 0, 1], dtype=np.int32)

state_vec = np.zeros((N,))
control_vec = np.zeros((N,))

# == Setup convolution ==
# * at Legendre-Gauss points for exact evaluation
lg_int_mda_mat = lg_int_mda(di)
lg_cm = LgConvPwcMat(lg_int_mda_mat, di)


# == Setup optimization ==
# * regularization
alpha = 5e-5
# * initial value
x = np.zeros(N, dtype=np.int32)
# * desired state for tracking objective
f_vec = f(di.lg_t_vec)
# * function handles
eval_f = lambda x: lg_objective_var(x, lg_cm, di, f_vec)
eval_jac = lambda x: lg_jacobian_var(x, lg_cm, di, f_vec)
# * algorithm control
Delta0 = N // 16
sigma = 1e-3
maxiter = 500
h = 2./N


# == Optimization with convolution evaluated at Legendre-Gauss points ==
opt_start = time.time()
xs = slip(eval_f, eval_jac, x, lo_bangs, alpha, h, Delta0, sigma, maxiter)

state_vec = lg_cm.conv(xs)
control_vec = xs

print("final objective = %.2e" % eval_f(control_vec))

plt.subplot(1, 2, 1)
plt.step(np.linspace(0., 1., N), control_vec)
plt.subplot(1, 2, 2)
plt.plot(np.linspace(0., 1., state_vec.shape[0]), state_vec)

plt.show()

