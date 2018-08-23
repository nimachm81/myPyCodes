## EvolAlgs.py
## Evolutionary algorithms tools and wrappers

__all__ = ["decodeBinToFloat", "codeFloatToBin", "decodeBinToFloat_multi",
            "codeFloatToBin_multi",
            "geneticAlg",
            "GAFloatArr", "GAFloatArrIndividual", "gaOptimizer"]


from ast import literal_eval
import numpy as np
import copy


def checkBinarySequence(seq_bin):
    for b in seq_bin:
        if not(b == 0 or b == 1):
            return False
    return True

def decodeBinToFloat(seq_bin, x_0, x_1):
    """
    It converts a sequence of zeros and ones to a float between x_0 and x_1
    """
    str1 = '0b'+''.join(map(str, seq_bin))
    b_val = float(literal_eval(str1))
    b_max = float(2**len(seq_bin) - 1)
    x = x_0 + b_val/b_max*(x_1 - x_0)
    return x
    
def codeFloatToBin(x, x_0, x_1, seq_len):
    """
    converts a float x between x_0 and x_1 to a binary sequence of length seq_len
    """
    if  not(x_0 <= x <= x_1):
        raise ValueError("x out of range")
    b_max = float(2**seq_len - 1)
    b_val = int((x - x_0)/(x_1 - x_0)*b_max)
    if b_val>b_max:
        b_val = int(b_max)
    seq_str = "{0:b}".format(b_val)
    chars = []
    chars.extend(seq_str)
    seq_bin = [literal_eval(c) for c in chars]
    seq_bin = [0]*(seq_len-len(seq_bin)) + seq_bin
    return seq_bin
    
    
def decodeBinToFloat_multi(seq_bin, intervals, seq_lens):
    if sum(seq_lens)!=len(seq_bin):
        raise ValueError("sum(seq_lens)!=len(seq_bin)")
    x = []
    start = 0
    for i in range(len(seq_lens)):
        x.append(decodeBinToFloat(seq_bin[start:start+seq_lens[i]],
                intervals[i][0], intervals[i][1]))
        start += seq_lens[i]
    return x
    
def codeFloatToBin_multi(x, intervals, seq_lens):
    b = []
    for i in range(len(x)):
        b = b + codeFloatToBin(x[i], intervals[i][0], intervals[i][1],
            seq_lens[i])
    return b
    

##--------------------------  Genetic algorithm ----------------
import random
from deap import base
from deap import creator
from deap import tools
import multiprocessing

def geneticAlg(f_eval, intervals, n_bits, pop_size=100, pr_cross=0.5, 
        pr_mut=0.2, n_gen=40, pr_bit_mut=0.05, tourno_size=3, 
        maximize=True, vbose=False, parallel=False, x_0=None):

    wghts = (1.0,)
    if not maximize:
        wghts = (-1.0,)
    creator.create("FitnessMax", base.Fitness, weights=wghts)
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    pool = None
    if parallel:
        print("initializing multiprocessing")
        pool = multiprocessing.Pool()
        toolbox.register("map", pool.map)

    # Attribute generator
    toolbox.register("attr_bool", random.randint, 0, 1)
    # Structure initializers
    toolbox.register("individual", tools.initRepeat, creator.Individual, 
        toolbox.attr_bool, sum(n_bits))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    def f_eval_bin(individual):
        x = decodeBinToFloat_multi(individual, intervals, n_bits)
        return [f_eval(x)]

    # Operator registering
    toolbox.register("evaluate", f_eval_bin)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=pr_bit_mut)
    toolbox.register("select", tools.selTournament, tournsize=tourno_size)

    random.seed()
    
    pop = toolbox.population(n=pop_size)
    CXPB, MUTPB, NGEN = pr_cross, pr_mut, n_gen
    
    if x_0!=None:
        pop[0].fitness.values = codeFloatToBin_multi(x_0, intervals, n_bits)
        
    print("Start of evolution")
    
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    
    if vbose:
        print("  Evaluated %i individuals" % len(pop))
    
    # Begin the evolution
    for g in range(NGEN):
        if vbose:
            print("-- Generation %i --" % g)
        
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
    
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
    
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        if vbose:
            print("  Evaluated %i individuals" % len(invalid_ind))
        
        # The population is entirely replaced by the offspring
        pop[:] = offspring
        
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        if vbose:
            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)
    print("-- End of (successful) evolution --")
    
    best_ind = tools.selBest(pop, 1)[0]
    x_best = decodeBinToFloat_multi(best_ind, intervals, n_bits)
    print("Best individual is %s, %s" % (x_best, best_ind.fitness.values))
    return [x_best, best_ind.fitness.values[0]]



##------------------------------------------------------------------------------



class EAindividual:
    def __init__(self, params):
        raise NotImplementedError()
    
    def SetContent(self, params):
        raise NotImplementedError()
        
    def Validate(self, params):
        raise NotImplementedError()
        
    def InitRandom(self, params):
        raise NotImplementedError()

    def ChangeIndivitual(self, params):
        raise NotImplementedError()
        
    def ChangeClone(self, params):
        raise NotImplementedError()

    def Evaluate(self, params):
        raise NotImplementedError()

    def GetValue(self, params):
        raise NotImplementedError()


class EApopulation:
    def __init__(self, popsize):
        self.pop = [None]*popsize

    def InitRandom(self, individualType, indiv_params=None, rand_params=None):
        """ individualType : the individual class name
            it will use this class to initialize the population: all members
            will be the same species
        """
        for i in range(len(self.pop)):
            self.pop[i] = individualType(indiv_params)
            self.pop[i].InitRandom(rand_params)
              
    def Evaluate(self, ev_params):       
        if "parallel" in ev_params and ev_params["parallel"]==True:
            raise NotImplementedError()
        else:
            if ev_params['pop']=="current-pop": 
                for i in range(len(self.pop)):
                    self.pop[i].Evaluate(ev_params)
            elif ev_params['pop']=="new-pop": 
                for i in range(len(self.pop_new)):
                    self.pop_new[i].Evaluate(ev_params)


    def SelectSubPopulation(self, sel_params):
        if sel_params['type']=='all':
            pop_size = len(self.pop)
            return range(pop_size)
        if sel_params['type']=='random':
            n_sel = sel_params['n_sel']   ## number of selections
            sel_inds = [None]*n_sel
            pop_size = len(self.pop)
            if sel_params['force-unique']==False:
                for i in range(n_sel):
                    sel_inds[i] = np.random.randint(0, pop_size)
            else:
                ##TODO: performance can be enhanced using a different selection algorithm
                assert sel_inds<0.5*pop_size   
                for i in range(n_sel):
                    sel_i = np.random.randint(0, pop_size)
                    while sel_i not in sel_inds:
                        sel_i = np.random.randint(0, pop_size)
                    sel_inds[i] = sel_i
            return sel_inds

    def ChangePopulation(self, pop_inds, change_params):
        pop_new = []
        for ind in pop_inds:
            indiv_new = self.pop[ind].ChangeClone(change_params)
            pop_new.append(indiv_new)
        self.pop_new = pop_new
    
    
    def TreatNewPopulation(self, select_params):
        raise NotImplementedError()
    

    def Run(self, alg_params):
        #print("Inside EApopulation")
        raise NotImplementedError()
        
##---------------------------------------------------------------

class GAFloatArrIndividual(EAindividual):
    def __init__(self, params):
        self.range = params["range"]
        Nx = len(self.range)
        self.x = np.zeros(Nx)
        self.fx = np.zeros(params['n_obj'])  ## fitness, n_obj:number of objectives
        self.valid = False
    
    def SetContent(self, params):
        x_0 = params['x_0']
        self.x = x_0
        self.fx = None
        
    def Validate(self, params=None):
        if params is None:
            valid = True
            for i in range(len(self.x)):
                x_0, x_1 = self.range[i]
                if self.x[i]<x_0 or self.x[i]>x_1:
                    valid = False
                    break
            self.valid = valid
            return valid
            
    def AdjustInvalid(self, params=None):
        if params is None:
            for i in range(len(self.x)):
                x_0, x_1 = self.range[i]
                if self.x[i]<x_0:
                    self.x[i] = x_0
                elif self.x[i]>x_1:
                    self.x[i] = x_1
            self.valid = True
                                                                                                                                  
        
    def InitRandom(self, params=None):
        if params is None:
            for i in range(len(self.x)):
                x_0, x_1 = self.range[i]
                self.x[i] = x_0 + np.random.rand()*(x_1 - x_0)
            

    def ChangeIndivitual(self, params):
        if params['op']=="mutate":
            mean = params['mean']
            sd = params['sd']  ## standard deviation
            for i in range(len(self.x)):
                x_0, x_1 = self.range[i]
                Dx_i = x_1-x_0
                self.x[i] += np.random.normal(loc=mean, scale=sd*Dx_i)
            if not self.Validate():
                self.AdjustInvalid()

    def ChangeClone(self, params):
        clone = copy.deepcopy(self)
        clone.ChangeIndivitual(params)
        return clone
        
    def Evaluate(self, params):
        raise NotImplementedError()
  

class GAFloatArr(EApopulation):
    
    def SetIndivitual(self, init_params):
        ind = init_params['ind']
        x_0 = init_params['x_0']
        assert ind < len(self.pop)
        self.pop[ind].SetContent({'x_0': x_0})

    def TreatNewPopulation(self, select_params):
        if select_params["method"]=="best-fit-all":
            ## select the best fit among all old+new population    
            pop_size = len(self.pop)
            pop_all = self.pop + self.pop_new
            ##single objective
            values_all = [pop_all[i].fx[0] for i in range(len(pop_all))]
            
            inds_sort = np.argsort(np.array(values_all))
            pop_next = [None]*pop_size
            if select_params["minmax"]=="max":
                for i in range(pop_size):
                    pop_next[i] = pop_all[inds_sort[-i-1]]
            else:
                for i in range(pop_size):
                    pop_next[i] = pop_all[inds_sort[i]]
            
            self.pop = pop_next


    def Run(self, alg_params):
        #print("Inside GAFloatArr")
        n_gen = alg_params["n_gen"]
        self.Evaluate(ev_params={'pop':"current-pop", "parallel":alg_params["parallel"]})       
        for i in range(n_gen):
            selection = self.SelectSubPopulation(sel_params={'type':'all'})
            change_params = {'op':"mutate", "mean":0.0, "sd":alg_params["sd"]}
            self.ChangePopulation(pop_inds=selection, change_params=change_params)
            self.Evaluate(ev_params={'pop':"new-pop", "parallel":alg_params["parallel"]})
            
            select_params = {"minmax":alg_params["minmax"], "method":alg_params["sel-method"]}
            self.TreatNewPopulation(select_params)
            
            print("{}/{}".format(i+1, n_gen), end=" ")
            
        return self.pop
            



def gaOptimizer(f_eval, intervals, pop_size=100, n_gen=40, mu_sd=0.2, 
        maximize=True, vbose=False, parallel=False, x_0=None):
    """ mu_sd: standard deviation for the mutation operator
        x_0: initial value (optional)
    """

    class funcOptimizerND(GAFloatArrIndividual):
        def Evaluate(self, params):
            fx = f_eval(self.x)
            self.fx = np.array([fx])
            return self.fx


    ga = GAFloatArr(popsize=pop_size)
    indiv_params = {"range":intervals, "n_obj":1}
    ga.InitRandom(individualType=funcOptimizerND, indiv_params=indiv_params)
    if x_0 is not None:
        ga.SetIndivitual(init_params={"ind":0, "x_0":x_0})
        
    minmax = "max"
    if not maximize:
        minmax = "min"
    alg_params = {"n_gen":n_gen, "minmax":minmax, "sel-method":"best-fit-all", "sd":mu_sd,
                  "parallel":parallel}
    pop = ga.Run(alg_params)

    opts = [(pop[i].x, pop[i].fx) for i in range(len(pop))]
    return opts






 
