import matplotlib.pyplot as plt
import numpy as np
import time

path = '../'
event_file_name = 'finalpr.data'
B_MULT = 'B_MULT'

event_file = open(path+event_file_name,'r').readlines()
list_nr_events = np.loadtxt(path+B_MULT,dtype=int,usecols=(2,))

#counters
in_elastic_event = 0
elastic_event = 0
next_index = 0
total_counted = 0
pt          = []
pt_weight = []
rapidity    = []
charged     = []
events      = np.asarray([])

for i,number in enumerate(list_nr_events[::]):
    #print('-'*100)
    if number == 2:
        in_elastic_event += 1
        #print('in-elastic event') 

    else:
        #temp_average = []
        elastic_event += 1
        #E, px, py, pz, C = np.loadtxt(event_file[next_index:next_index+number],
        #        unpack=True, usecols=(4,5,6,7,13))

        #p_transverse = np.sqrt(px**2+py**2)
        #pt_average = np.sum(p_transverse)/len(p_transverse)
        #pt_weight.extend([pt_average]*number)
        #pt.extend(p_transverse)
        #rapidity.extend(0.5*np.log((E + pz)/(E - pz)))
        #charged.extend(C)
        
        for particle in event_file[next_index:next_index+number]:
            particle = list(map(float,particle.strip().split()))

            E = particle[4]
            px = particle[5]
            py = particle[6]
            pz = particle[7]
            C = particle[13]
            if C != 0:
                charged.append(C)
                p_transverse = (np.sqrt(px**2 + py**2))
                pt_weight.append(1./(2*np.pi*p_transverse))
                pt.append(p_transverse)
                 
                if (E-pz)==0:
                    rapidity.append(10)
                elif (E+pz)==0:
                    rapidity.append(0)
                else:
                    rapidity.append(0.5*np.log((E + pz)/(E - pz)))
        

        #pt_average = (np.sum(temp_average)/len(temp_average))
        #pt_weight.extend([pt_average]*number) 
    
    if i%200==0:
        print(i)
        #print(np.sum(temp_average),len(temp_average))
        #print(next_index, next_index+number, number)
    
    next_index += number 
    #print('-'*100)

def get_hist (x, bins=100,ranges=None, weights=None):
    #fig, ax = plt.subplots()

    n, bin_edges = np.histogram(x,bins=bins,range=ranges,weights=weights)
    #ax_pt.hist(x,bins=bins,histtype=histtype,log=log)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

    #print(len(bin_centers),len(n))
    #ax.plot(bin_centers,n,'-')
    #if log:
    #    ax_ptweight.set_yscale('log')

    #ax_ptweight.errorbar(bin_centers,
    #        n,
    #        yerr=np.std(x),
    #        linestyle='None',
    #        ecolor='red')

    return bin_centers, n
    #plt.savefig(filename)

################################################################################
#HISTOGRAMS
################################################################################
pt          = np.asarray(pt) 
pt_weight   = np.asarray(pt_weight) 
rapidity    = np.asarray(rapidity) 
charged     = np.asarray(charged) 
print('in-elastic events = ', in_elastic_event)
print('elastic events = ', elastic_event)
#pt_weight = pt/pt_weight
pt_range = (0,10)

pt_max = np.max(pt)
pt_min = np.min(pt)

rapidity_limits = np.where(np.abs(rapidity) < 2.5)
#pt_in_range = np.zeros(len(rapidity_limits))
#pt_weight_in_range = np.zeros(len(rapidity_limits))
#for i,index in enumerate(rapidity_limits):
#    pt_in_range[i] = pt[index]
#    pt_weight_in_range[i] = pt_weight[index]


pt = np.array([pt[i] for i in rapidity_limits])
pt_weight = np.array([pt_weight[i] for i in rapidity_limits])

bins = np.logspace(-1,np.log10(np.max(pt)),num=101)
pt_x, pt_y = get_hist(pt,bins)

pt_weighted_x, pt_weighted_y = get_hist(pt,bins,weights=pt_weight)
        
#pt_weighted_sum = np.array([ np.sum(pt_weight    for i in pt_weighted_x[1::]])
#rapidity_x, rapidity_y = get_hist(x=rapidity,bins=500)

################################################################################
#PLOTS
################################################################################
fig_pt, ax_pt = plt.subplots()

ax_pt.plot(pt_x,pt_y,'-')#,label='Pt')
ax_pt.set_yscale('log')
ax_pt.set_xlim(pt_range)
ax_pt.set_title(r'$p_t$, charged particles with abs(rapidity) < 2.5 ')
#ax_pt.legend(loc=0)
ax_pt.set_ylabel(r'$\frac{dN}{dPt}$')#,rotation=0)
ax_pt.text(0.5,1,'pt max = %f\npt min = %f'%(pt_max,pt_min))

step = 5 
ax_pt.errorbar(pt_x[::step],
        pt_y[::step],
        yerr=1./np.sqrt(pt_y[::step]),
        linestyle='None',
        ecolor='red')

plt.savefig('pt.png')
fig_ptweight, ax_ptweight = plt.subplots()

ax_ptweight.plot(pt_weighted_x,pt_weighted_y,'-')#,label='Pt weighted')
#ax_ptweight.set_title('Pt weighted')
ax_ptweight.set_yscale('log')
ax_ptweight.set_xlim(pt_range)
#ax_ptweight.legend(loc=0)
ax_ptweight.set_title('Weighted $p_t$ with abs(rapidity) < 2.5')
ax_ptweight.set_ylabel(r'$\frac{1}{2\pi Pt}\frac{dN}{dPt}$')

ax_ptweight.errorbar(pt_weighted_x[::step],
        pt_weighted_y[::step],
        yerr=1./np.sqrt(pt_weighted_y[::step]),#*pt_weight[::step],
        linestyle='None',
        ecolor='red')

plt.savefig('pt_weight.png')
#ax[2].plot(rapidity_x,rapidity_y,'-',label='Rapidity')
##ax[2].set_title('Rapidity')
#ax[2].set_yscale('log')
#ax[2].set_ylabel(r'$\frac{dN}{dy}$',rotation=0)
#ax[2].legend(loc=0)

#plt.savefig('all.png')

#fig2, ax2 = plt.subplots(1)
#
#ax2.hist(charged)
#
#plt.savefig('charge.png')



