import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import itertools


class Bacterial_population:
    def __init__(self, initial_pop_size, loci=1, psi_max_s=1, psi_min=-1, gamma=0.2, MIC = 1, cost_vector=[0], benefit_vector=[0], mutation_rate=0.001, living_style=0, biofilm_benefit=1, biofilm_cost=0.1, release_rate=0, adhesion_rate=0):

        ### define basic genetic parameters

        self.loci=loci
        self.psi_max_s=psi_max_s
        self.psi_min=psi_min
        self.gamma=gamma

        self.cost_vector=cost_vector.copy()
        self.benefit_vector=benefit_vector.copy()
        self.mutation_rate=mutation_rate
        self.MIC=MIC
        self.population_sizes=[]

        ### lifestyle associated parameters
        self.living_style=living_style
        self.biofilm_benefit=biofilm_benefit
        self.biofilm_cost=biofilm_cost
        self.release_rate=release_rate
        self.adhesion_rate=adhesion_rate

        ## define initial pop sizes:
        if len(initial_pop_size)==2**loci:
            self.initial_pop_size=initial_pop_size.copy()
        elif  len(initial_pop_size)==1:
            self.initial_pop_size=np.zeros(2**loci)
            self.initial_pop_size[0]=initial_pop_size[0]
        else:
            print ('problem with initial sizes')

        ### determine variables, genotypes and dependencies
        self.number_of_genotypes=2**loci
        self.genotypes=self.generate_genotypes()
        self.gen_shape = np.shape(self.genotypes)
        self.genotype_costs = 1 - np.prod(1 - (self.genotypes * (self.cost_vector)), axis=1)
        ### additive benefits:
        #self.genotype_benefits = np.dot(self.genotypes, self.benefit_vector)
        ## multiplicative  benefits
        self.genotype_benefits=np.prod((self.genotypes*(np.array(self.benefit_vector)-1))+1, axis=1)
        self.deps = self.define_dependencies()
        self.neg_deps = self.define_negative_dependencies()

        ### determine benefits
        if self.living_style==0: #this is for plankton
            #self.z_mic = self.MIC * (1 + (self.genotype_benefits) )  #this is for addititve
            self.z_mic = self.MIC * ( self.genotype_benefits )  ## this is for multiplicative
            self.total_cost=self.genotype_costs
            self.psi_max=self.psi_max_s*(1-self.total_cost)
        elif self.living_style==1: #this is for biofilm
            #self.z_mic = self.MIC * (1 + (self.genotype_benefits) + self.biofilm_benefit) #this is for addititve
            self.z_mic = self.MIC * ((self.genotype_benefits) + self.biofilm_benefit)## this is for multiplicative
            self.total_cost=1-(1-self.genotype_costs)*(1-self.biofilm_cost)
            self.psi_max=self.psi_max_s*(1-self.total_cost)
            self.psi_min=self.psi_min*(1-self.biofilm_cost)
            ### here i should also define reduction in psi min!
        else:
            print ('error in defining living style')
        self.population_sizes.append(self.initial_pop_size)

    def generate_genotypes(self):
        ### generating all possible diploid genotypes (0-1 number of alleles) from the number of loci.
        genotypes = np.empty([0, self.loci])
        for seq in itertools.product("01", repeat=self.loci):
            s = np.array(seq)
            s = list(map(int, s))
            genotypes = np.vstack([genotypes, s])
        return genotypes.astype(int)

    def define_dependencies(self):
        dependencies = []
        x_x = [2**i for i in range(self.gen_shape[1])]
        for i in range(self.gen_shape[0]):
            dependencies.append([])
            for j in x_x:
                sused = i ^ j
                if i > sused:
                    dependencies[i].append(sused)
        return (dependencies)

    def define_negative_dependencies(self):
        dependencies = []
        x_x = [2**i for i in range(self.gen_shape[1])]
        for i in range(self.gen_shape[0]):
            dependencies.append([])
            for j in x_x:
                sused = i ^ j
                if i < sused:
                    dependencies[i].append(sused)
        return (dependencies)


# In[ ]:


class Stochastic_utility():
    def __init__(self, population_pla, population_bio, my_treatment):
        print('initiated')
        self.tau=0.1
        self.stoich_mat=self.build_stoichiometric_matrix(population_pla.deps, population_pla.genotypes)
        current_population=np.concatenate((population_pla.initial_pop_size,population_bio.initial_pop_size))
        current_concentration=0
        self.propensities=self.calculate_propensities(current_population, current_concentration, my_treatment,population_pla, population_bio)

    def GillespieTauLeap(self, my_treatment, tau, V, y0, a_conc,population_pla, population_bio):
        '''N is the stoichiometry matrix
        rate function V. function that gives us the propensities as a function of the current state.
        V is the k reactivities or the equiibrium associations. Note V(y) gives the matrix we are looking for when starting the beginning of the algorithm.
        y0 is the initial condition/state vector.
        tlen is max length of time. We will build up to this.
        '''
        tlen=my_treatment.cycle_length
        N= self.stoich_mat
        t = 0.0  # starting time
        ts = [0.0]  # reporting array, to save information
        y = np.copy(y0)  # using the given initial condition
        # lists because these will be dynamically resizing bcs we will be randomly
        # choosing our time intervals. We could pre-allocate these in a np.array
        # to make more efficient.
        res = [list(y)]
        a_conc_array = [a_conc]
        while True:  # just continuously looping until there is a break from within
            #print ('time',t)
            a_conc_t = a_conc*np.exp(my_treatment.degradation_rate*t)  # is this the place where concentration is degraded?!
            #print ('this is concentration', a_conc_t)

            prop = V(y, a_conc_t)  # propensities
            #print ('these are prop', prop)
            a0 = sum(prop)  # to see total propensity
            if a0 == 0:
                print('propensities are 0')
                break
            dt = tau

            if a0<0:
                print ('propensities are negative', a0)


            if t + dt > tlen:  # if waiting time will exceed time limit
                break
            for idx in range(len(prop)):
                if prop[idx]<0:
                    print ('propensities', prop[idx], prop)
                if prop[idx]>=0:
                    how_many=np.random.poisson(prop[idx]*dt)  #this needs to be changed to poisson
                else:
                    print('this is wrong', prop)

                    print ('at index', idx)
                    print ('this', prop[idx])
                    print ('dt', dt)
                    print ('weird', prop[idx]*dt)
                    print ('random poisson', np.random.poisson(prop[idx]*dt))
                    break
                how_many=np.random.poisson(prop[idx]*dt)  #this needs to be changed to poisson
                change_to_apply = how_many*N[:, idx]
                change_to_apply.shape = len(change_to_apply)  # converting to 1D array
                y += change_to_apply  # this is a np.array
            if y.any()<0:
                print ('negative y')  #anchor negative y
                y[y<0]=0

            # Adding the time
            t += dt
            # saving the time and results so that we can use it later
            ts.append(t)
            res.append(list(y))
            a_conc_array.append(a_conc_t)

        ts=np.array(ts)
        a_conc_array=np.array(a_conc_array)
        return(ts, np.array(res), a_conc_array)

    def build_stoichiometric_matrix(self, deps, genotypes):
        no_genotypes = np.size(genotypes, 0)
        no_deps = int(np.sum(genotypes))
        mut_matrix = np.zeros([no_genotypes, no_deps])
        counter = 0
        for i in range(no_genotypes):
            #print (deps[i])
            for j in range(len(deps[i])):
                mut_matrix[deps[i][j], counter] = -1
                mut_matrix[i, counter] = 1
                counter += 1
        biofilm_release_matrix = np.zeros([2 * no_genotypes, no_genotypes])
        biofilm_adhesion_matrix = np.zeros([2 * no_genotypes, no_genotypes])

        for i in range(no_genotypes):
            biofilm_release_matrix[i, i] = 1
            biofilm_release_matrix[i + no_genotypes, i] = -1
            biofilm_adhesion_matrix[i, i] = -1
            biofilm_adhesion_matrix[i + no_genotypes, i] = 1

        stoich_matrix_simple = np.concatenate(
            (np.eye(no_genotypes), mut_matrix, -np.eye(no_genotypes)), axis=1)  # this is just
        # a very simple matrix of growth+mutation+death, for one type: either plankton or biofilm. To make the complete matrix,
        # we need to make it 4x larger, where 1_1 part and 2_2 quarants are these
        # simple ones, and the other quadrants are empty
        complete_stoich_top = np.concatenate(
            (stoich_matrix_simple, 0 * stoich_matrix_simple), axis=1)
        complete_stoich_bottom = np.concatenate(
            (0 * stoich_matrix_simple, stoich_matrix_simple), axis=1)
        complete_stoich = np.concatenate(
            (complete_stoich_top, complete_stoich_bottom), axis=0)
        # adding biofilm release
        complete_stoich = np.concatenate(
            (complete_stoich, biofilm_release_matrix, biofilm_adhesion_matrix), axis=1)
        return complete_stoich

#    def calculate_propensities(self,y, my_treatment, population_pla, population_bio, current_concentration):
    def calculate_propensities(self,y,a_conc, my_treatment, population_pla, population_bio):
        K=my_treatment.car_cap
        kappa=my_treatment.kappa
        #a_conc=current_concentration
        #print (a_conc)
        #y = np.copy(y0)
        #print ('this y shape', np.shape(y))
        y[y < 0] = 0
        #print ('first y', y)
        all_y = np.sum(y)  # this is whole pop size
        reactions = []

        # Dealing with plankton
        gamma=population_pla.gamma
        psi_min=population_pla.psi_min
        mu=population_pla.mutation_rate
        deps=population_pla.deps
        z_mic=population_pla.z_mic

        ### some shit happening with psi, psi max aand gamma definitions!
        psi=population_pla.psi_max
        psi_max = psi * (1 - all_y / K) - gamma
        #print ('this is pla', z_mic, psi)

        # start with growth rates
        for i in range(int(len(y) / 2)):
            reactions.append(y[i] * psi[i] * (1 - all_y / K))  # seems ok


        # mutations
        for i in range(int(len(y) / 2)):
            for j in (deps[i]):
                reactions.append(mu * y[j] * psi[j] * (1 - all_y / K))  # seems ok

        # deaths
        for i in range(int(len(y) / 2)):
            reactions.append(y[i] *
                             (gamma +
                              (psi_max[i] -
                               psi_min) *
                              (a_conc /
                               z_mic[i]) ** kappa /
                              ((a_conc /
                                z_mic[i])**kappa -
                                  psi_min /
                                  psi_max[i])))
        # Dealing with biofilm
        gamma=population_bio.gamma
        psi_min=population_bio.psi_min
        mu=population_bio.mutation_rate
        deps=population_bio.deps
        z_mic=population_bio.z_mic

        ### some shit happening with psi, psi max aand gamma definitions!
        psi=population_bio.psi_max
        psi_max = psi * (1 - all_y / K) - gamma

        add_for_bio=population_pla.number_of_genotypes
        # start with growth rates
        for i in range(int(len(y) / 2)):
            reactions.append(y[i+add_for_bio] * psi[i] * (1 - all_y / K))  # seems ok

        # mutations
        for i in range(int(len(y) / 2)):
            for j in (deps[i]):
                reactions.append(mu * y[j+add_for_bio] * psi[j] * (1 - all_y / K))  # seems ok
        # deaths
        for i in range(int(len(y) / 2)):
            reactions.append(y[i+add_for_bio] *
                             (gamma +
                              (psi_max[i] -
                               psi_min) *
                              (a_conc /
                               z_mic[i]) ** kappa /
                              ((a_conc /
                                z_mic[i])**kappa -
                                  psi_min /
                                  psi_max[i])))
        # and no dealing with biofilm releasing cells and adhesion
        for i in range(int(len(y) / 2), len(y)):
            reactions.append(y[i] * population_pla.release_rate)  # seems ok
        for i in range(0,int(len(y) / 2)):
            reactions.append(y[i] * population_bio.adhesion_rate)  # seems ok

        #print ('all reactions', reactions)
        #if np.isnan([reactions]).any():
        #    print ('reactions', reactions)
        #    print ('last y',  y)

        return np.array(reactions)  ##anchor




# In[ ]:


class Treatment():
    ### definition of the treatment
    def __init__(self, concentration_gradient, cycle_length=24*60, cycle_number=1, car_cap=10**9,kappa=1.5, degradation_rate=0, dilution_factor_pla=1, dilution_factor_bio=1 ):

        self.concentration_gradient=concentration_gradient
        self.cycle_length=cycle_length
        self.cycle_number=cycle_number

        self.simulated_population_sizes=[]
        self.dilution_factor_record=[]
        self.car_cap=car_cap
        self.kappa=kappa
        self.degradation_rate=degradation_rate
        self.dilution_factor_pla=dilution_factor_pla
        self.dilution_factor_bio=dilution_factor_bio

        #self.simulated_population_sizes.append(population.population_sizes)
        #self.dilution_factor_record.append(0)
        self.time=[]#np.arange(self.cycle_number+1)*self.cycle_length

    def run_deterministic_simulation(self, population_pla, population_bio):
        print ('running deterministic simulation')
        self.concentration_record=[]
        self.pla_population_record=np.empty([0, population_pla.number_of_genotypes])
        self.bio_population_record=np.empty([0, population_bio.number_of_genotypes])

        starting_pla_pop=population_pla.population_sizes[0]
        starting_bio_pop=population_bio.population_sizes[0]
        self.time=[]
        for i in range (self.cycle_number):
            concentration=self.concentration_gradient[i]
            starting_pla_pop[starting_pla_pop<1]=0
            starting_bio_pop[starting_bio_pop<1]=0
            #self.simulated_population_sizes= np.append(a_conc, y_0)
            z_0=np.append(concentration,starting_pla_pop )
            z_0=np.append(z_0,starting_bio_pop)
            t = np.linspace(0, self.cycle_length, 1000)
            self.time=np.append(self.time, t+i*self.cycle_length)
            current_conc =concentration
            #print ('we are in cycle', i, 'at concentration', concentration, 'with z0',z_0)
            [my_z,infodict] = odeint(self.calculate_rates, z_0,t,
                                args=(population_pla, population_bio, current_conc), full_output = 1)
            this_conc=my_z[:,0]

            self.concentration_record=np.concatenate((self.concentration_record,this_conc ))
            this_pop_pla=my_z[:,1:population_pla.number_of_genotypes+1]
            self.pla_population_record=np.concatenate((self.pla_population_record,this_pop_pla))
            this_pop_bio=my_z[:,population_bio.number_of_genotypes+1:]
            self.bio_population_record=np.concatenate((self.bio_population_record,this_pop_bio))

            #dilution:
            dilution_coef_pla=np.sum(this_pop_pla[-1,:])/np.sum(starting_pla_pop)
            dilution_coef_bio=np.sum(this_pop_bio[-1,:])/np.sum(starting_bio_pop)

            if dilution_coef_pla<1:
                dilution_coef_pla=1
            if dilution_coef_bio<1:
                dilution_coef_bio=1

            starting_pla_pop = this_pop_pla[-1,:]/dilution_coef_pla
            starting_bio_pop = this_pop_bio[-1,:]/dilution_coef_bio


        self.final_populations=np.append(self.pla_population_record[-1,:],self.bio_population_record[-1,:])

    def calculate_rates(self, z,t,population_pla, population_bio, current_conc ) : ### still missing release rate!!!
        ############################## Antibiotic concentration #############
        current_conc=z[0]
        if current_conc<0:
            #print (current_conc)
            current_conc=0

        dcdt=[current_conc*self.degradation_rate]
        ############################## Population 0 (plankton) ###############

        pla_pop=z[1:population_pla.number_of_genotypes+1]
        pla_pop[pla_pop<0]=0

        bio_pop=z[population_pla.number_of_genotypes+1:]
        bio_pop[bio_pop<0]=0
        total_pop= np.sum(pla_pop)+np.sum(bio_pop)
        dydt_plankton=self.population_rates(population_pla,pla_pop,current_conc,total_pop )
        new_release=population_bio.release_rate*bio_pop
        dydt_plankton=dydt_plankton+new_release

        dydt_biofilm=self.population_rates(population_bio,bio_pop,current_conc,total_pop )
        new_adhesion=population_pla.adhesion_rate*pla_pop
        dydt_biofilm=dydt_biofilm+new_adhesion

        dzdt= np.append(dcdt,dydt_plankton)
        dzdt=np.append(dzdt, dydt_biofilm)
        return dzdt

    def population_rates(self,population,this_pop,current_conc, total_pop):
                ############################## Population 1 (biofilm) ###############
        ############################## Psi a Psi max ###############
        #total_pop=np.sum(this_pop)

        # this is max growth rate

        pla_psi_max = population.psi_max_s * (1 - population.total_cost)*(1 - total_pop / self.car_cap) - population.gamma  # this is doubling rate
        growth = this_pop* pla_psi_max
        death = this_pop  * (population.gamma + (pla_psi_max - population.psi_min) * (current_conc / population.z_mic)
                                ** self.kappa / ((current_conc / population.z_mic)**self.kappa - population.psi_min /pla_psi_max))
        outgoing_mut = np.zeros(population.number_of_genotypes)
        incoming_mut = np.zeros(population.number_of_genotypes)

        this_neg_dep = np.array(population.neg_deps[0])
        outgoing_mut[0] = population.mutation_rate* growth[0] * len(this_neg_dep)

        for i in range(1, population.number_of_genotypes):
            this_dep = np.array(population.deps[i])
            this_neg_dep = np.array(population.neg_deps[i])
            outgoing_mut[i] = population.mutation_rate * growth[i] * len(this_neg_dep)
            incoming_mut[i] = population.mutation_rate* np.sum(growth[this_dep])

        dydt = growth - death + incoming_mut - outgoing_mut #+ release_rate * y_biofilm
        return dydt

    def plot_results(self, name=''):
        XX=name
        time=self.time/60
        plt.plot(time,self.concentration_record)
        plt.xlabel('time [h]')
        plt.ylabel('concentration')
        plt.title('this is concentration')
        plt.xlim([0,self.cycle_length*self.cycle_number/(60)])
        plt.show()

        plt.semilogy(time, self.pla_population_record)
        plt.xlabel('time [h]')
        plt.ylabel('pop size')
        plt.title('this is plankton')
        plt.xlim([0,self.cycle_length*self.cycle_number/(60)])
        plt.ylim([1,self.car_cap*1.5])
        #plt.savefig(XX+'Plankton.png',bbox_inches='tight')
        plt.show()

        plt.semilogy(time, self.bio_population_record)
        plt.title('this is biofilm')
        plt.xlabel('time [h]')
        plt.ylabel('pop size')
        plt.xlim([0,self.cycle_length*self.cycle_number/(60)])
        plt.ylim([1,self.car_cap*1.5])
        #plt.savefig(XX+'Biofilm.png',bbox_inches='tight')
        plt.show()

        plt.semilogy(time, self.pla_population_record)
        plt.semilogy(time, self.bio_population_record,'--')
        plt.title('this is both')
        plt.xlabel('time [h]')
        plt.ylabel('pop size')
        plt.xlim([0,self.cycle_length*self.cycle_number/(60)])
        plt.ylim([1,self.car_cap*1.5])
        plt.savefig(XX+'Biofilm.png',bbox_inches='tight')
        plt.show()

    def run_stochastic_simulation(self, population_pla, population_bio):
        print ('running stochastic simulation')
        stoch_sim=Stochastic_utility(population_pla, population_bio, self)
        print (stoch_sim.stoich_mat)
        tau=1

        self.concentration_record=[]
        self.pla_population_record=np.empty([0, population_pla.number_of_genotypes])
        self.bio_population_record=np.empty([0, population_bio.number_of_genotypes])

        starting_pla_pop=population_pla.population_sizes[0]
        starting_bio_pop=population_bio.population_sizes[0]
        self.time=[]

        for i in range (self.cycle_number):
            concentration=self.concentration_gradient[i]
            Vy = lambda y, a_conc: stoch_sim.calculate_propensities(y, a_conc, self, population_pla, population_bio)
            #Vy = lambda y, a_conc: stoch_sim.calculate_propensities(y, self, population_pla, population_bio, concentration)
            starting_pla_pop[starting_pla_pop<1]=0
            starting_bio_pop[starting_bio_pop<1]=0
            #creating starting point pla and biofilms
            z_0=np.append(starting_pla_pop, starting_bio_pop )
            #print ('we are in cycle', i, 'at concentration', concentration, 'with z0',z_0)
            (ts, my_z, a_conc_array2)= stoch_sim.GillespieTauLeap(self, tau, Vy, z_0, concentration, population_pla, population_bio)
            this_pop_pla = my_z[:, 0:population_pla.number_of_genotypes ]
            this_pop_bio = my_z[:, population_pla.number_of_genotypes:]
            this_conc=[concentration]

            ## check that stuff is non - zero
            if this_pop_pla.any()<0:
                print ('pla negative', this_pop_pla)
                this_pop_pla[this_pop_pla<0]=0
            if this_pop_bio.any()<0:
                print ('bio negative', this_pop_bio)
                this_pop_bio[this_pop_bio<0]=0



            self.time=np.append(self.time, ts+i*self.cycle_length)
            self.concentration_record=np.concatenate((self.concentration_record,a_conc_array2 ))
            self.pla_population_record=np.concatenate((self.pla_population_record,this_pop_pla))
            self.bio_population_record=np.concatenate((self.bio_population_record,this_pop_bio))


            #dilution:

            #print ('dilutions', self.dilution_factor_bio, self.dilution_factor_pla)
            if self.dilution_factor_pla==0:
                #print ('dilution pla 0: diluting to original population')
                if np.sum(this_pop_pla[-1,:])>0:
                    dil_fact_pla=np.sum(population_pla.initial_pop_size)/np.sum(this_pop_pla[-1,:])
                else:
                    dil_fact_pla=0

                #print ('dilution pla', dil_fact_pla, np.sum(this_pop_pla[-1,:]))

                if dil_fact_pla>1:
                    dil_fact_pla=1
                if dil_fact_pla<0:
                    print('negative dilution. Final pop size',np.sum(this_pop_pla[-1,:]) )
                    dil_fact_pla=0
                #print ('dilution pla 2', dil_fact_pla)
                starting_pla_pop =this_pop_pla[-1,:]*dil_fact_pla
                #print ('dilution is', dil_fact_pla)
            else:
                starting_pla_pop = this_pop_pla[-1,:]*self.dilution_factor_pla

            if self.dilution_factor_bio==0:
                #print ('dilution bio 0: diluting to original population')

                if np.sum(this_pop_bio[-1,:])>0:
                    dil_fact_bio=np.sum(population_bio.initial_pop_size)/np.sum(this_pop_bio[-1,:])
                else:
                    dil_fact_bio=0

                #print ('dilution bio', dil_fact_bio, np.sum(this_pop_bio[-1,:]))
                if dil_fact_bio>1:
                    dil_fact_bio=1
                if dil_fact_bio<0:
                    print('negative dilution. Final pop size',np.sum(this_pop_bio[-1,:]) )
                    dil_fact_bio=0
                #print ('dilution bio 2', dil_fact_bio)
                starting_bio_pop =this_pop_bio[-1,:]*dil_fact_bio
                #print ('dilution is', dil_fact_bio)
                #print ('new starting pop is', np.sum(starting_bio_pop))
            else:
                starting_bio_pop = this_pop_bio[-1,:]*self.dilution_factor_bio



            starting_pla_pop=np.round(starting_pla_pop)
            starting_bio_pop=np.round(starting_bio_pop)
            #print ('start pop pla', starting_pla_pop)
            #print ('start pop bio', starting_bio_pop)

        self.final_populations=np.append(self.pla_population_record[-1,:],self.bio_population_record[-1,:])
        #print (self.pla_population_record)

    def save_output(self, name='test0'):
        name_p=name+'_pla'
        name_bio=name+'_bio'
        np.savetxt(name_p+'_pop.txt', self.pla_population_record)
        np.savetxt(name_bio+'_pop.txt', self.bio_population_record)
        #np.savetxt(name+'_time.txt', self.time)
        #np.savetxt(name+'_conc.txt', self.concentration_record)

        print ('all saved as', name)


# In[ ]:


class Analysis_functions():
    def __init__(self, name='test', loci=4):
        self.active=1
        self.loci=loci
        self.name=name
        print ('analysis active')

    def generate_genotypes(self):
        ### generating all possible diploid genotypes (0-1 number of alleles) from the number of loci.
        self.genotypes = np.empty([0, self.loci])
        for seq in itertools.product("01", repeat=self.loci):
            s = np.array(seq)
            s = list(map(int, s))
            self.genotypes = np.vstack([self.genotypes, s])

    def genMeanMutNum(self,time, pop, k=4):
        self.generate_genotypes()
        sum_mut=np.sum(self.genotypes,1)
        weighted_pops=pop*sum_mut
        #mean_mut2=np.sum(weighted_pops,1)/np.sum(pop,1)
        a=np.sum(weighted_pops,1)
        b=np.sum(pop,1)
        #mean_mut=np.divide(a,b, out=np.zeros_like(a), where=b!=0)
        mean_mut=a/b#(a,b, out=np.zeros_like(a), where=b!=0)
        #plt.plot(time, mean_mut)
        return (mean_mut)

    def plot_total_populations(self, N=1):
        name0=self.name
        fin_pops=np.zeros(32)
        for my_id in range(N):
            name=name0+'pla'
            name_pla_pop=name+'pla_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            time=np.load(name_time)
            time=time/60
            pla_pops=np.load(name_pla_pop)
            plt.semilogy(time, np.sum(pla_pops,1), color='cyan', alpha=0.2)
            plt.ylim([1,2*10**9])
            #get_final_pops(self, -1, pla_pops,k=4)

        for my_id in range(N):
            name=name0+'bio'
            name_bio_pop=name+'bio_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            time=np.load(name_time)
            time=time/60
            bio_pops=np.load(name_bio_pop)
            #print (np.shape(bio_pops))
            #print (np.shape(np.sum(bio_pops,1)))
            plt.semilogy(time,np.sum(bio_pops,1) , ':m', linewidth=2, alpha=0.1)
            plt.ylim([1,2*10**9])
            #get_final_pops(self, -1, pla_pops,k=4)
        plt.xlim([180,500])

        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('population size', fontsize=15)
        plt.savefig(name0+'Popsize.png',dpi=200,bbox_inches='tight')

        plt.show()
        name_conc=name+'conc'+str(my_id)+'.npy'
        #time=np.load(name_time)
        conc=np.load(name_conc)
        plt.xlim([180,500])
        plt.plot(time, conc, color='red',linewidth=2 )
        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('concentration', fontsize=15)
        plt.savefig(name0+'Conc.png',dpi=200,bbox_inches='tight')

    def get_final_pops(self, time, pop,k=4):
        [t,phe]=np.shape(pop)
        sizes=np.zeros(phe)
        print ('this is last time', time[-1])

        if time[-1]<576:
            print ('dead')
        else:
            sizes=pop[-1,:]
        return(sizes)

    def get_composition(self, N=1):
        name0=self.name
        fin_pops=np.zeros(32)
        for my_id in range(N):
            name=name0+'bio'
            name_bio_pop=name+'bio_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            time=np.load(name_time)
            time=time/60
            bio_pops=np.load(name_bio_pop)
            final_pops=self.get_final_pops(time, bio_pops,k=4)
            print (final_pops)
            #print (np.shape(bio_pops))
            #print (np.shape(np.sum(bio_pops,1)))
            #plt.semilogy(time,np.sum(bio_pops,1) , color='green', linewidth=2, alpha=0.1)
            plt.ylim([1,2*10**9])
            #get_final_pops(self, -1, pla_pops,k=4)
        plt.xlim([180,500])

    def plot_average_mutation_number(self, N=1):
        fig1, ax1 = plt.subplots(dpi=2000)
        plt.axvline(x=12, color='red', linestyle='-', linewidth=84, alpha=0.05)
        plt.axvline(x=36, color='red', linestyle='-', linewidth=84, alpha=0.1)
        plt.axvline(x=60, color='red', linestyle='-', linewidth=84, alpha=0.15)
        plt.axvline(x=84, color='red', linestyle='-', linewidth=84, alpha=0.2)
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        name0=self.name
        fin_pops=np.zeros(32)
        for my_id in range(N):
            name=name0+'pla'
            name_pla_pop=name+'pla_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            pla_pops=np.load(name_pla_pop)
            time=np.load(name_time)
            time=time/60
            concentration=np.load(name_conc)
            pla_mean_mut_num=self.genMeanMutNum(time, pla_pops, k=4)
            plt.plot(time, pla_mean_mut_num,'c', linewidth=2, alpha=0.08)
        #plt.xlim([0,96])
        #plt.ylim([0,4.1])
        #plt.savefig(name0+'pla.png',dpi=200,bbox_inches='tight')
        #plt.show()
        for my_id in range(N):
            name=name0+'bio'
            name_bio_pop=name+'bio_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            bio_pops=np.load(name_bio_pop)
            time=np.load(name_time)
            time=time/60
            concentration=np.load(name_conc)
            bio_mean_mut_num=self.genMeanMutNum(time, bio_pops, k=4)
            plt.plot(time, bio_mean_mut_num, ':m', linewidth=2, alpha=0.08)


        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('number of mutations', fontsize=15)
        #ax1.yaxis.set_major_formatter(formatter)
        #plt.savefig(savename0+'MutNumT2.png', dpi=200,bbox_inches='tight')
        #plt.show()
        plt.xlim([0,96])
        plt.ylim([-0.1,4.1])
        plt.savefig(name0+'both.png',dpi=200,bbox_inches='tight')
        plt.show()

        fig1, ax1 = plt.subplots(dpi=2000)
        plt.axvline(x=12, color='red', linestyle='-', linewidth=84, alpha=0.05)
        plt.axvline(x=36, color='red', linestyle='-', linewidth=84, alpha=0.1)
        plt.axvline(x=60, color='red', linestyle='-', linewidth=84, alpha=0.15)
        plt.axvline(x=84, color='red', linestyle='-', linewidth=84, alpha=0.2)
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        for my_id in range(N):
            name=name0+'bio'
            print (my_id)
            name_bio_pop=name+'bio_pop'+str(my_id)+'.npy'
            print (name_bio_pop)
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            bio_pops=np.load(name_bio_pop)
            time=np.load(name_time)
            time=time/60
            concentration=np.load(name_conc)
            bio_mean_mut_num=self.genMeanMutNum(time, bio_pops, k=4)
            plt.plot(time, bio_mean_mut_num, ':m', linewidth=2, alpha=0.08)


        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('number of mutations', fontsize=15)
        #ax1.yaxis.set_major_formatter(formatter)
        #plt.savefig(savename0+'MutNumT2.png', dpi=200,bbox_inches='tight')
        #plt.show()
        plt.xlim([0,96])
        plt.ylim([-0.001,0.04])
        plt.savefig(name0+'bio.png',dpi=200,bbox_inches='tight')
        plt.show()

        print ('done')
        return(fin_pops/N)

    def plot_some_dynamics(self, N=1):
        fig, ax = plt.subplots(dpi=2000)
        #ax.set_prop_cycle('color', [plt.cm.winter(i) for i in np.linspace(0, 1, 2**4)])
        ax.set_prop_cycle('color', [plt.cm.jet(j) for j in np.linspace(0, 1, 2**4)])
        name0=self.name
        fin_pops=np.zeros(32)
        for my_id in range(N):
            name=name0+'pla'
            name_pla_pop=name+'pla_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            time=np.load(name_time)
            time=time/60
            pla_pops=np.load(name_pla_pop)
            plt.semilogy(time, pla_pops,  alpha=0.2)
            plt.ylim([1,2*10**9])
            #get_final_pops(self, -1, pla_pops,k=4)
        #plt.xlim([180,500])
        plt.xlim([0,500])
        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('population size', fontsize=15)
        plt.savefig(name0+'PopsizeAllPla.png',dpi=200,bbox_inches='tight')
        plt.show()
        fig, ax = plt.subplots(dpi=2000)
        ax.set_prop_cycle('color', [plt.cm.jet(j) for j in np.linspace(0, 1, 2**4)])
        for my_id in range(N):
            name=name0+'bio'
            name_bio_pop=name+'bio_pop'+str(my_id)+'.npy'
            name_time=name+'time'+str(my_id)+'.npy'
            name_conc=name+'conc'+str(my_id)+'.npy'
            time=np.load(name_time)
            time=time/60
            bio_pops=np.load(name_bio_pop)
            #print (np.shape(bio_pops))
            #print (np.shape(np.sum(bio_pops,1)))
            plt.semilogy(time,(bio_pops) , linewidth=1, alpha=0.1)
            plt.ylim([1,2*10**9])
            #get_final_pops(self, -1, pla_pops,k=4)
        plt.xlim([0,500])
        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('population size', fontsize=15)
        plt.savefig(name0+'PopsizeAllBio.png',dpi=200,bbox_inches='tight')
        plt.show()
        name_conc=name+'conc'+str(my_id)+'.npy'
        conc=np.load(name_conc)
        plt.xlim([0,500])
        plt.plot(time, conc, color='red',linewidth=1 )
        plt.xlabel('time [h]', fontsize=15)
        plt.ylabel('concentration', fontsize=15)
        plt.savefig(name0+'Conc.png',dpi=200,bbox_inches='tight')
