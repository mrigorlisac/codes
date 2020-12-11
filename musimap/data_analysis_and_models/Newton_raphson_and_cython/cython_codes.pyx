import numpy as np
import scipy.special as sp
import scipy.special.cython_special as csc


def NRiteration_cy(float x, float n, float F, float p, float G, int gamma, float nu, float b):
    
    cdef float qc, dxqc, mue, lae, dqmue, dqlae, dxmue, dxlae, f, fp, res
    
    qc = (G/F)**(1/gamma)*x**((p-n)/gamma) # qc is normalized critical streamflow
    dxqc = qc*(p-n)/gamma/x # dxqc is derivative of qc with respect to x
    
    

    if b==1:
        # compute mue = mu_epsilon
        mue = (sp.gamma(1/nu+gamma)/sp.gamma(1/nu)/nu**(-gamma)*
                       sp.gammaincc(1/nu+gamma,qc/nu))
        # copmute lae = lambda_epsilon
        lae = (sp.gammaincc(1/nu,qc/nu))
        # compute dqmue = derivative of mu_epsilon with respect to qc
        dqmue = (-np.exp(-qc/nu)*(qc/nu)**(1/n+gamma-1)/sp.gamma(1/nu)/nu**(-gamma+1))
        # compute dqlae = derivative of lambda_epsilon with respect to qc
        dqlae = (-np.exp(-qc/nu)*(qc/nu)**(1/n-1)/sp.gamma(1/nu)/nu)

    elif b==2:
        # compute mue = mu_epsilon
        mue = (sp.gamma(1/nu+1-gamma)/sp.gamma(1/nu)/nu**(gamma-1)*
                       sp.gammainc(1/nu+1-gamma,1/qc/nu))
        # copmute lae = lambda_epsilon
        lae = (sp.gammainc(1/nu+1,1/qc/nu))
        # compute dqmue = derivative of mu_epsilon with respect to qc
        dqmue = (-np.exp(-1/qc/nu)*(1/qc/nu)**(1/nu-gamma+2)/sp.gamma(1/nu)/nu**(gamma-2))
        # compute dqlae = derivative of lambda_epsilon with respect to qc
        dqlae = (-np.exp(-1/qc/nu)*nu*(1/qc/nu)**(1/nu+2)/sp.gamma(1/nu+1))

    dxmue = dqmue*dxqc # dxmue is derivative of mu_epsilon with respect to x
    dxlae = dqlae*dxqc # dxlae is derivative of lambda_epsilon with respect to x

    # compute f(x) the value of the function for which we search the root
    f = x - 1 + F*mue*x**n - G*lae*x**p
    # compute the derivative of f(x) with respect to x at x
    fp = 1 + F*(mue*n*x**(n-1) + dxmue*x**n) - G*(lae*p*x**(p-1) + dxlae*x**p)

    # result is an improved value of x, the root, following Newton-Raphson scheme
    res = x - np.where(f>0, f/fp, 0)

    return res


def NRiteration_cy1(float x, float n, float F, float p, float G, float gamma, float nu, float b):
    
    cpdef float qc, dxqc, mue, lae, dqmue, dqlae, dxmue, dxlae, f, fp, res
    
    qc = (G/F)**(1/gamma)*x**((p-n)/gamma) # qc is normalized critical streamflow
    dxqc = qc*(p-n)/gamma/x # dxqc is derivative of qc with respect to x
    
    if b==1:
        # compute mue = mu_epsilon
        mue = (csc.gamma(1/nu+gamma)/csc.gamma(1/nu)/nu**(-gamma)*csc.gammaincc(1/nu+gamma,qc/nu))
        # compute lae = lambda_epsilon
        lae = (csc.gammaincc(1/nu,qc/nu))
        # compute dqmue = derivative of mu_epsilon with respect to qc
        dqmue = (-np.exp(-qc/nu)*(qc/nu)**(1/n+gamma-1)/csc.gamma(1/nu)/nu**(-gamma+1))
        # compute dqlae = derivative of lambda_epsilon with respect to qc
        dqlae = (-np.exp(-qc/nu)*(qc/nu)**(1/n-1)/csc.gamma(1/nu)/nu)

    elif b==2:
        # compute mue = mu_epsilon
        mue = (csc.gamma(1/nu+1-gamma)/csc.gamma(1/nu)/nu**(gamma-1)*
                       csc.gammainc(1/nu+1-gamma,1/qc/nu))
        # copmute lae = lambda_epsilon
        lae = (csc.gammainc(1/nu+1,1/qc/nu))
        # compute dqmue = derivative of mu_epsilon with respect to qc
        dqmue = (-np.exp(-1/qc/nu)*(1/qc/nu)**(1/nu-gamma+2)/csc.gamma(1/nu)/nu**(gamma-2))
        # compute dqlae = derivative of lambda_epsilon with respect to qc
        dqlae = (-np.exp(-1/qc/nu)*nu*(1/qc/nu)**(1/nu+2)/csc.gamma(1/nu+1))

    dxmue = dqmue*dxqc # dxmue is derivative of mu_epsilon with respect to x
    dxlae = dqlae*dxqc # dxlae is derivative of lambda_epsilon with respect to x

    # compute f(x) the value of the function for which we search the root
    f = x - 1 + F*mue*x**n - G*lae*x**p
    # compute the derivative of f(x) with respect to x at x
    fp = 1 + F*(mue*n*x**(n-1) + dxmue*x**n) - G*(lae*p*x**(p-1) + dxlae*x**p)

    # result is an improved value of x, the root, following Newton-Raphson scheme
    res = x - np.where(f>0, f/fp, 0)

    return res



def solve_SPL_implicit(int stack, float variability, float hf_out, float zf_out, float n, Ff, p, Gf, gamma, b, dh0f, SlopeM, dx, receivers):
    
    solution = np.zeros_like(zf_out)

    for i in stack:
        tol = 1e-3
        err = tol*2 # initializes error
        variab = variability.flatten()
        while np.abs(err)>tol: # Newton-Raphson iterative loop
#            zk = NRiteration(zf_out[i], n, Ff[i], p, Gf[i], gamma, variab[i], b) # updates normalized height
            x = zf_out[i]
            nu = variab[i]
            qc = (Gf[i]/Ff[i])**(1/gamma)*x**((p-n)/gamma) # qc is normalized critical streamflow
            dxqc = qc*(p-n)/gamma/x # dxqc is derivative of qc with respect to x

            if b==1:
                # compute mue = mu_epsilon
                mue = (csc.gamma(1/nu+gamma)/csc.gamma(1/nu)/nu**(-gamma)*
                               csc.gammaincc(1/nu+gamma,qc/nu))
                # copmute lae = lambda_epsilon
                lae = (csc.gammaincc(1/nu,qc/nu))
                # compute dqmue = derivative of mu_epsilon with respect to qc
                dqmue = (-np.exp(-qc/nu)*(qc/nu)**(1/n+gamma-1)/csc.gamma(1/nu)/nu**(-gamma+1))
                # compute dqlae = derivative of lambda_epsilon with respect to qc
                dqlae = (-np.exp(-qc/nu)*(qc/nu)**(1/n-1)/csc.gamma(1/nu)/nu)

            elif b==2:
                # compute mue = mu_epsilon
                mue = (csc.gamma(1/nu+1-gamma)/csc.gamma(1/nu)/nu**(gamma-1)*
                               csc.gammaincc(1/nu+1-gamma,1/qc/nu))
                # copmute lae = lambda_epsilon
                lae = (csc.gammaincc(1/nu+1,1/qc/nu))
                # compute dqmue = derivative of mu_epsilon with respect to qc
                dqmue = (-np.exp(-1/qc/nu)*(1/qc/nu)**(1/nu-gamma+2)/csc.gamma(1/nu)/nu**(gamma-2))
                # compute dqlae = derivative of lambda_epsilon with respect to qc
                dqlae = (-np.exp(-1/qc/nu)*nu*(1/qc/nu)**(1/nu+2)/csc.gamma(1/nu+1))

            dxmue = dqmue*dxqc # dxmue is derivative of mu_epsilon with respect to x
            dxlae = dqlae*dxqc # dxlae is derivative of lambda_epsilon with respect to x

            # compute f(x) the value of the function for which we search the root
            f = x - 1 + Ff[i]*mue*x**n - Gf[i]*lae*x**p
            # compute the derivative of f(x) with respect to x at x
            fp = 1 + Ff[i]*(mue*n*x**(n-1) + dxmue*x**n) - Gf[i]*(lae*p*x**(p-1) + dxlae*x**p)

            # result is an improved value of x, the root, following Newton-Raphson scheme
            zk = x - np.where(f>0, f/fp, 0)

            err = np.abs(zf_out[i] - zk) # computes error
            zf_out[i] = zk # # updates normalized height

        solution[i] = hf_out[receivers[i]] + min(zf_out[i]*dh0f[i], SlopeM*dx) # updates height once convergence is reached
        
    return solution






def some_function(some_variable):
    pass



