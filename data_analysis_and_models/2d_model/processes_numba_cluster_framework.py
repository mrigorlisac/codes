# Snippet of the 2D model process used in the GPU cluster framework. Here, numba is used to speed up critical processes
# that were identified as being slow in execution (gamma functions...). The whole model uses many different connected
# processes that are dependent on each other. The processes represent e.g. soil response to rainfall, translation of rainfall
# to discharge etc...

from numba.extending import get_cython_function_address
from numba import vectorize, njit
import ctypes
import numpy as np
import xsimlab as xs

addr = get_cython_function_address("scipy.special.cython_special", "gammaincc")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
gammaincc_fn = functype(addr)


@vectorize("float64(float64, float64)")
def vec_gammaincc(x, y):
    return gammaincc_fn(x, y)


@njit
def gammaincc_in_njit(x, y):
    return vec_gammaincc(x, y)


addr = get_cython_function_address("scipy.special.cython_special", "__pyx_fuse_1gamma")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
gamma_fn = functype(addr)


@vectorize("float64(float64)")
def vec_gamma(x):
    return gamma_fn(x)


@njit
def gamma_in_njit(x):
    return vec_gamma(x)


@njit()
def solve_SPL_implicit(
    stack,
    variability,
    hf_out,
    zf_out,
    n,
    Ff,
    p,
    Gf,
    gamma,
    b,
    dh0f,
    SlopeM,
    dx,
    receivers,
):

    solution = np.zeros_like(zf_out)

    for i in stack:
        tol = 1e-3
        err = tol * 2  # initializes error
        variab = variability.flatten()
        while np.abs(err) > tol:  # Newton-Raphson iterative loop
            #            zk = NRiteration(zf_out[i], n, Ff[i], p, Gf[i], gamma, variab[i], b) # updates normalized height
            x = zf_out[i]
            nu = variab[i]
            qc = (Gf[i] / Ff[i]) ** (1 / gamma) * x ** (
                (p - n) / gamma
            )  # qc is normalized critical streamflow
            dxqc = (
                qc * (p - n) / gamma / x
            )  # dxqc is derivative of qc with respect to x

            if b == 1:
                # compute mue = mu_epsilon
                mue = (
                    gamma_in_njit(1 / nu + gamma)
                    / gamma_in_njit(1 / nu)
                    / nu ** (-gamma)
                    * gammaincc_in_njit(1 / nu + gamma, qc / nu)
                )
                # copmute lae = lambda_epsilon
                lae = gammaincc_in_njit(1 / nu, qc / nu)
                # compute dqmue = derivative of mu_epsilon with respect to qc
                dqmue = (
                    -np.exp(-qc / nu)
                    * (qc / nu) ** (1 / n + gamma - 1)
                    / gamma_in_njit(1 / nu)
                    / nu ** (-gamma + 1)
                )
                # compute dqlae = derivative of lambda_epsilon with respect to qc
                dqlae = (
                    -np.exp(-qc / nu)
                    * (qc / nu) ** (1 / n - 1)
                    / gamma_in_njit(1 / nu)
                    / nu
                )

            elif b == 2:
                # compute mue = mu_epsilon
                mue = (
                    gamma_in_njit(1 / nu + 1 - gamma)
                    / gamma_in_njit(1 / nu)
                    / nu ** (gamma - 1)
                    * gammaincc_in_njit(1 / nu + 1 - gamma, 1 / qc / nu)
                )
                # copmute lae = lambda_epsilon
                lae = gammaincc_in_njit(1 / nu + 1, 1 / qc / nu)
                # compute dqmue = derivative of mu_epsilon with respect to qc
                dqmue = (
                    -np.exp(-1 / qc / nu)
                    * (1 / qc / nu) ** (1 / nu - gamma + 2)
                    / gamma_in_njit(1 / nu)
                    / nu ** (gamma - 2)
                )
                # compute dqlae = derivative of lambda_epsilon with respect to qc
                dqlae = (
                    -np.exp(-1 / qc / nu)
                    * nu
                    * (1 / qc / nu) ** (1 / nu + 2)
                    / gamma_in_njit(1 / nu + 1)
                )

            dxmue = dqmue * dxqc  # dxmue is derivative of mu_epsilon with respect to x
            dxlae = (
                dqlae * dxqc
            )  # dxlae is derivative of lambda_epsilon with respect to x

            # compute f(x) the value of the function for which we search the root
            f = x - 1 + Ff[i] * mue * x ** n - Gf[i] * lae * x ** p
            # compute the derivative of f(x) with respect to x at x
            fp = (
                1
                + Ff[i] * (mue * n * x ** (n - 1) + dxmue * x ** n)
                - Gf[i] * (lae * p * x ** (p - 1) + dxlae * x ** p)
            )

            # result is an improved value of x, the root, following Newton-Raphson scheme
            zk = x - np.where(f > 0, f / fp, 0)

            err = np.abs(zf_out[i] - zk)  # computes error
            zf_out[i] = zk  # updates normalized height

        solution[i] = hf_out[receivers[i]] + min(
            zf_out[i] * dh0f[i], SlopeM * dx
        )  # updates height once convergence is reached

    return solution


from fastscape.processes.flow import DrainageArea, SingleFlowRouter
