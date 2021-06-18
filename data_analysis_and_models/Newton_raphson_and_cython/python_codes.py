import numpy as np
import scipy.special as sp


def NRiteration_py(x, n, F, p, G, gamma, nu, b):

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