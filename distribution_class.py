import matplotlib.pyplot as plt
import numpy as np

class Distribution_plot:

    def __init__(
            self, path_to_event=None, 
            path_to_B_MULT=None,
            part_start=0,
            part_stop=0,
            part_step=1,
            print_step=500
            ):
        

        self.event_file  = open(path_to_event,'r').readlines()
        self.nr_of_particles = np.loadtxt(path_to_B_MULT,dtype=int,usecols=(2,))

        if part_stop == 0:
            self.part_stop  = len(self.nr_of_particles)
        else:
            self.part_stop = part_stop
        self.part_start     = part_start
        self.part_step      = part_step 
        self.print_step     = print_step

        self.in_elastic_event   = 0
        self.elastic_event      = 0
        self.pt                 = []
        self.pt_weight          = []
        self.rapidity           = []
        self.charge             = []
        self.pt_rap_lim         = []
        self.pt_weight_rap_lim  = []

    def fill_arrays(self):
        event_file = self.event_file
        nr_particles = self.nr_of_particles
        part_start = self.part_start
        part_stop = self.part_stop
        part_step = self.part_step
        print_step = self.print_step
        
        previous_index = 0

        for i,number in enumerate(nr_particles[part_start:part_stop:part_step]):
            if number == 2:
                self.in_elastic_event += 1
            else:
                self.elastic_event += 1

                for particle in event_file[previous_index : previous_index+number]:
                    particle = list(map(float,particle.strip().split()))

                    E = particle[4]
                    px = particle[5]
                    py = particle[6]
                    pz = particle[7]
                    C = particle[13]
                    if C != 0:
                        self.charge.append(C)
                        p_transverse = (np.sqrt(px**2 + py**2))
                        self.pt_weight.append(1./(2*np.pi*p_transverse))
                        self.pt.append(p_transverse)
                         
                        if (E-pz)==0:
                            self.rapidity.append(10)
                        elif (E+pz)==0:
                            self.rapidity.append(0)
                        else:
                            self.rapidity.append(0.5*np.log((E + pz)/(E - pz)))

            previous_index += number

            if i%print_step==0:
                print(i)

        self.pt = np.asarray(self.pt) 
        self.pt_weight = np.asarray(self.pt_weight)
        self.rapidity = np.asarray(self.rapidity) 
        self.charge = np.asarray(self.charge) 
    
    def find_rap_lim(self, absolute_limit = 1):
    
        rap_lim = np.where(np.abs(self.rapidity < absolute_limit))
        self.pt_rap_lim = np.array([self.pt[i] for i in rap_lim]) 
        self.pt_weight_rap_lim = np.array([self.pt_weight[i] for i in rap_lim]) 
    
    def create_hist(self, x, bins=100, ranges=None, weights=None):
        
        n, bin_edges = np.histogram(
                x,
                bins=bins,
                range=ranges,
                weights=weights,
                )
         
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

        return bin_centers, n, bin_edges
    
    def create_plot(self,x,y):
        fig, ax = plt.subplots()
        
        ax.plot(x,y,'-')
    
        return ax
         
if __name__=='__main__':
    path                    = '../build/data/'
    event_file_name         = 'finalpr.data'
    B_MULT_file             = 'B_MULT'

    EVENT = Distribution_plot(
            path_to_event   = path+event_file_name,
            path_to_B_MULT  = path+B_MULT_file,
            part_start      =0,
            part_stop       =0,
            part_step       =1,
            print_step      =500
            )

    EVENT.fill_arrays()
    EVENT.find_rap_lim(absolute_limit=2.5) 

    bins = np.logspace(-1,np.log10(np.max(EVENT.pt)),num=501)

    pt_x, pt_N, pt_bin_edges = EVENT.create_hist(
            x=EVENT.pt_rap_lim,
            bins=bins)

    pt_weight_x, pt_weight_y, pt_weight_bin_edges = EVENT.create_hist(
            x=EVENT.pt_rap_lim,
            bins=bins,
            weights=EVENT.pt_weight_rap_lim)

    fig_pt, ax_pt = plt.subplots()
    ax_pt.plot(pt_weight_x, pt_weight_y, '-')
    ax_pt.set_yscale('log')
    ax_pt.set_xlim((0,10))
    ax_pt.set_title(r'charged particles with $|y|< 2.5$ ')
    ax_pt.set_ylabel(r'$\frac{1}{2\pi p_T}\frac{dN}{p_T}$')
    ax_pt.set_xlabel(r'$p_T$')

    step = 35 
    pt_N[pt_N==0] = 1000000

    ax_pt.errorbar(pt_weight_x[::step],
            pt_weight_y[::step],
            yerr=(1./np.sqrt(pt_N[::step]))*pt_weight_y[::step],
            linestyle='None',
            ecolor='red')
    
    plt.savefig('class_pt.png') 

    print('DONE')
