import matplotlib.pyplot as plt
import numpy as np



def peerlings(k,beta,e0):
    return 0

def damage(p,q,damage,rt,rc,rs):
    pr = 0
    pr = np.copy(p)
    pr[p>0] = p[p>0] * (1 - (damage * rc))
    pr[p<=0] = p[p<=0] * (1 - (damage * rt))
    qr = q * (1 - (damage * rs))
    return pr,qr

xmin, xmax, ymin, ymax = -2, 5, 0, 5
ticks_frequency = 1

A = 1
B = 1
p = np.linspace(0,100,100) - A
q = np.linspace(0,100,100) * B


plt.plot(p, q,label="Undamaged yield surface")
for d in [0.25,0.5,0.9,0.99,0.99999999]:
    p_damaged,q_damaged = damage(p,q,d,1,0.5,0.9)
    plt.plot(p_damaged, q_damaged,label="d = {}".format(d))
plt.legend()
ax = plt.gca()
ax.set(xlim=(xmin-1, xmax+1), ylim=(ymin, ymax+1), aspect='equal')
#ax.set(aspect='equal')
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('$x$', size=14, labelpad=-24, x=1.02)
ax.set_ylabel('$y$', size=14, labelpad=-21, y=1.02, rotation=0)
ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

plt.show()
