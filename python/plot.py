import numpy as np
import matplotlib.pyplot as plt

def histogram(data,bins=100,range=None,weights=None):
  n, bin_edges = np.histogram(data,bins=bins,range=range,weights=weights)
  bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
  return n, bin_centers

NPOM = '2'
filepath = '../counted/S'+NPOM+'/'
fontsize = '25'
'''
nf = np.loadtxt(filepath+'dnf.out')
nbnf = np.loadtxt(filepath+'nbdnf.out')
nbdn_nf = np.loadtxt(filepath+'nbdnf_dnf.out')
#nf = np.ma.masked_equal(nf,0)
fig = plt.figure()
fig.suptitle('NPOMS = '+NPOM, fontsize=fontsize)
ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(122)

ax1.plot(nf)
ax1.set_ylabel('$\\frac{dn}{dn_F}$',fontsize=fontsize,rotation=0)
ax1.yaxis.labelpad = 30
ax1.set_xlabel('$n_F$',fontsize=fontsize)
ax1.set_yscale('log')

ax2.plot(nbnf)
ax2.set_ylabel('$n_B\\frac{dn}{dn_F}$',fontsize=fontsize,rotation=0)
ax2.yaxis.labelpad = 30
ax2.set_xlabel('$n_F$',fontsize=fontsize)
ax2.set_yscale('log')

ax3.plot(nbdn_nf)
ax3.set_ylabel('$\\frac{n_B\\frac{dn}{dn_F}}{\\frac{dn}{dn_F}}$',fontsize=fontsize,rotation=0)
ax3.yaxis.labelpad = 30
ax3.set_xlabel('$n_F$',fontsize=fontsize)
#ax3.set_yscale('log')

#plt.tight_layout()
plt.show()
'''
fig,ax = plt.subplots()
fig.suptitle('NPOMS = '+NPOM,fontsize=fontsize)
#ax.set_xlim(0,20)
#ax.set_ylim(10,15)
ax.set_xlabel('$n_F$',fontsize=fontsize)
ax.set_ylabel('$<n_B(n_F)>$',fontsize=fontsize,rotation=90)
ax.yaxis.labelpad = 0
for H in range(5):
  nbdn_nf = np.loadtxt(filepath+'H'+str(H)+'nbdnf_dnf.out')
  ax.plot(nbdn_nf,label='H='+str(H))
tot = np.loadtxt(filepath+'tot_sum.out')
ax.plot(tot,label='tot sum',marker='o',color='black')
plt.legend()
plt.show()

