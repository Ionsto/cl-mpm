import numpy as np
import matplotlib.pyplot as plt


def MC(eps1,eps2,E,nu):
    eps = [eps1,eps2,0]
    de3 = (E/((1+nu) * (1-(2*nu)))) * np.array([[1-nu,nu,nu],
                                                [nu,1-nu,nu],
                                                [nu,nu,1-nu]])
    sig = np.matmul(de3,eps)
    princ = np.sort(sig)[::-1]
    s1 = sig[0]
    s3 = sig[2]
    # print(princ)
    sn = princ[0] + princ[2]
    st = princ[0] - princ[2]
    angle = 50
    f = st - (sn * np.sin(angle * 3.14/180))
    return f


angle = 50
c = 1

eps1 = np.linspace(-10,10)
eps2 = np.linspace(-10,10)
samples = len(eps1)
# f = np.zeros((samples,samples))
# for i in range(samples):
#     for j in range(samples):
#         f[i,j] = MC(eps1[i],eps2[i],1,0.25)

f = np.zeros((samples))
for i in range(samples):
    f[i] = MC(eps1[i],0,1,0.25)

# plt.contour(eps1,eps2,f)
plt.scatter(eps1,f)
plt.show()
