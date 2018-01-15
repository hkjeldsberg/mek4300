import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import quad

n = 100
w_1 = np.zeros(n)
w_2 = np.zeros(n)
v_1 = np.zeros(n)
v_2 = np.zeros(n)

eta_N = 2
alpha = 4
# Initial conditions
w_1[-1] = 0
v_1[-1] = np.exp(-alpha*eta_N)

w_2[-1] = 1
v_2[-1] = 0
dx = 2. / n

def U(eta):
    return 2*eta - 2*eta**3 + eta**4

def ddU(eta):
    return 2*eta*(eta-1)


def c_quad(func, a, b, **kwargs):
    if func is callable:
        def real_func(x):
            return scipy.real(func(x))
        def imag_func(x):
            return scipy.imag(func(x))

        real_integral = quad(real_func, a, b, **kwargs)
        imag_integral = quad(imag_func, a, b, **kwargs)
    else:
        print func
        real_integral = quad(func.real, a, b, **kwargs)
        imag_integral = quad(func.imag, a, b, **kwargs)

    return real_integral[0] + 1j*imag_integral[0]

cr = 2
ci = 0.01
for n in range(n-1, 0, -1):
    x = n*dx
    UC = (U(x) - cr)**2 + ci**2
    # Part 1
    vstar1 = -dx*w_1[n] + v_1[n]
    fn1 = ( alpha**2 + (ddU(x) * U(x) - ddU(x) * cr + 1j*ci*ddU(x)) / UC ) * v_1[n]
    fstar1 = ( alpha**2 + (ddU(x) * U(x) - ddU(x) * cr + 1j*ci*ddU(x)) / UC ) * vstar1
    wstar1 = -dx*fn1  + w_1[n]

    v_1[n-1] = v_1[n] - dx*0.5*(wstar1 + w_1[n])
    w_1[n-1] = w_1[n] - dx*0.5*(fstar1+ fn1)

    # Part 2
    vstar2 = -dx*w_2[n] + v_2[n]
    fn2 = ( alpha**2 + (ddU(x) * U(x) - ddU(x) * cr + 1j*ci*ddU(x)) / UC ) * v_2[n]
    fstar2 = ( alpha**2 + (ddU(x) * U(x) - ddU(x) * cr + 1j*ci*ddU(x)) / UC ) * vstar2
    wstar2 = -dx*fn2  + w_2[n]

    v_2[n-1] = v_2[n] - dx*0.5*(wstar2 + w_2[n])
    w_2[n-1] = w_2[n] - dx*0.5*(fstar2+ fn2)


    """
    A = c_quad(v_1[n-1]*np.conj(v_1[n-1]), 0,2)
    B = c_quad(ddU*U / UC, 0,2)
    C = c_quad(ddU / UC, 0,2)
    D = c_quad(w[n-1]*np.conj(w[n-1]), 0,2)
    c = (alpha**2*A + B - D) / C
    cr = c.real
    ci = c.imag
    """

sol = w_1 + w_1[0] / w_2[0] * w_2
t_ = np.linspace(0,2, 100)
plt.plot(t_, sol.real)
plt.show()


