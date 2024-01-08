# (spebnr.py) last update: 2023.05.18


import numpy.random as rnd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines
from statistics import mean



num_round = 6 # useful to round the distance


"""
species_levels: list of species (with min and max level)
process_name: name of a process
It finds the index of the species (process_name).
return: index (this is useful also to have an ordered list of parallel processes)
"""
def get_index( species_levels, process_name ):
    
    species = ""
    split_p_name = process_name.split('_')
    species = split_p_name[0]
    
    if len(split_p_name)>2:
        i = 1
        while i < len(split_p_name)-1:
            species += '_' + split_p_name[i]
            i += 1
            
    index = next((i for i, item in enumerate(species_levels) if item[0] == species), -1)
    
    return index

# get_index: end



### EVOLUTION SEQUENCES

"""
pdef: map of all processes
p: current process
It controls the map and all possible reactions that p can do.
return: temp_map (map of all possible processes),
    possible_reactions (possible reaction that p can do)
"""
def new_constr( pdef, p ):

    temp_map = []
    for process in pdef: # for each process
        if process['process_name'] in p:
            temp_map.append(process)

    # selects only possible reactions
    possible_reactions = set()
    is_first_set = True
    
    for process in temp_map:
        set_of_action = set()
        for action in process['reactions']:
            set_of_action.add(action['reaction_id'])
        if is_first_set:
            possible_reactions = set_of_action
            is_first_set = False
        else:
            possible_reactions=set_of_action.intersection(possible_reactions)
            
    return temp_map, possible_reactions

# new_constr: end


"""
pdef: map of all processes
p: set of parallel processes
r: list of reactions
It calls new_constr which returns all possible reactions to do.
return: all triple (possible_reaction, weight, next_process)
    possible_reactions (the name of possible reactions to to)
"""
def process( pdef, p, r ):
    
    temp_map, possible_reactions = new_constr(pdef, p)
    
    if len(possible_reactions) == 0: # if there is no possible reaction
        #print('Possible reactions: none') # possible reactions to do
        return [],possible_reactions

    result = []
    for process in temp_map:
        for reaction in process['reactions']:
            if reaction['reaction_id'] in possible_reactions:
                result.append(reaction)
                
    return result, possible_reactions

# process: end


"""
reaction_name: the name of a reaction 
r: list of reactions
return: the intrinsic value (k) of the reaction reaction_name
"""
def rate_of_reaction( reaction_name, r ):

    for i in range(0, len(r)):
        if list(r[i].keys())[0] == reaction_name:
            return r[i][reaction_name]
        
# rate_of_reaction: end


"""
list_prob: list of triples (probability, reaction_id, next_process)
It calculates the reaction to do.
return: the selected triple (probability, reaction_id, next_process)
"""
def sample_element_from_list( list_prob ):
    
    #print('list_prob: ',list_prob)
    
    # create a list with the sum of prob.
    # e.g. we have 3 objects:
    # 1) prob = 0.2
    # 2) prob = 0.3
    # 3) prob = 0.5
    # list = [ 0.2, 0.5, 1.0 ]
    total_prob = 0
    list_total_prob = []
    for i in range(0,len(list_prob)):
        if i == len(list_prob)-1:
            list_total_prob.append(1)
        else:
            total_prob += list_prob[i]['prob']
            list_total_prob.append(total_prob)

    # generate a random number
    u = rnd.random()

    # return the corresponding triple
    # e.g. u=0.15 --> return the triple at position 0
    for i in range(len(list_prob)):
        if u<=list_total_prob[i]:
            return list_prob[i]
        
# sample_element_from_list: end
        

"""
pdef: map of all processes
p: set of parallel processes
r: list of reaction
s: list of species (with min and max level)
It takes all triple - about each reaction - returned by process.
It calculates the probability of each reaction and chooses a reaction to do.
It returns the next reaction to do.
return: reaction_id (chosen reaction),
    next_processes (set of next parallel processes)
"""
def pstep( pdef, p, r, s ): 
    
    all_reactions, reactions_name = process(pdef, p, r)
    
    if len(reactions_name) == 0: # if there is no possible reaction
        return '',p
    
    # creates a list of dictionaries
    # each dictionary contains:
    #   the name of the reaction
    #   the total weight of the reaction (initialized to 1)
    #   the list of next processes
    list_actions = []
    for reaction_name in reactions_name:
        dictionary = {}
        dictionary['reaction_id'] = reaction_name
        dictionary['total_weight'] = rate_of_reaction(reaction_name, r) # equals to the rate 
        dictionary['next_processes'] = [ [] for i in s ] # creates an empty list to insert each specie
        list_actions.append(dictionary)
        
    # for each reaction in all_reactions
    #   calculates the right weight and inserts it in list_actions
    #   adds the next process which is in list_actions
    # does not control the role because the weight in pdef is equal to 1 (where role!='R')
    for reaction in all_reactions:
        reaction_id = reaction['reaction_id']
        weight = reaction['weight']
        next_process = reaction['next_process']

        # for each possible reaction, adds the weight and the next_processes
        index = next((i for i, item in enumerate(list_actions) if item["reaction_id"] == reaction_id), -1)
        list_actions[index]['total_weight'] = list_actions[index]['total_weight'] * weight
        index_next_process = get_index(s,next_process)
        list_actions[index]['next_processes'][index_next_process] = next_process

    # calculates the total weight
    total_weight_for_all_actions = 0
    for action in list_actions:
        total_weight_for_all_actions += action['total_weight']
       
    # return the triple (probability, reaction_id, next_processes)
    list_prob = []
    for reaction in list_actions:
        dictionary = {}

        # dictionary['prob'] = (weight of the reaction)/(sum of the weight of all reactions)
        dictionary['prob'] = reaction['total_weight']/total_weight_for_all_actions

        dictionary['reaction_id'] = reaction['reaction_id']
        dictionary['next_processes'] = reaction['next_processes']
        list_prob.append(dictionary)
        
    next_action = sample_element_from_list(list_prob)
    next_processes = next_action['next_processes']
    #print('Chosen action: ',next_action['reaction_id'])
    
    return reaction_id, next_processes

# pstep: end


"""
Used only for debugging.
"""
def add_occurrences( occurrences, key, element, n ):
    
    index = next((i for i, item in enumerate(occurrences) if item[key] == element), -1)
    if index == -1:
        dictionary = {}
        dictionary[key] = element
        dictionary['occurrences'] = 1
        occurrences.append(dictionary)
    else:
        occurrences[index]['occurrences'] += n
        
    return occurrences

# add_occurrences: end


"""
lastMax: an integer number
p: set of parallel processes
return: the maximum level between lastMax and the maximum level in p
"""
def maxLevelInP( lastMax, p ):
    
    newMax = lastMax
    for process in p:
        split = process.split('_')
        if int(split[len(split)-1]) > newMax:
            newMax = int(split[len(split)-1])
            
    return newMax

# maxLevelInP: end


"""
pdef: map of all processes
p: set of parallel processes
r: list of reaction
s: list of species (with min and max level)
h: number of steps
It simulates h steps for only one run.
For each step, it calls pstep which returns the next reaction to do.
return: result (the list of steps and the occurrences),
     occurrences (only for debugging),
     maxLevelInP (the maximum level of a species in p)
"""
def run( pdef, p, r, s, h ):
    
    result = []
    
    occurrences = []
    #occurrences = add_occurrences(occurrences, 'p', p, 1)
    
    _maxLevelInP = maxLevelInP(0, p)
    
    #print('\np 0: ',p,'\n')
    for i in range(1,h+1):

        #print('STEP ',i,' on p',i-1)
        action, p = pstep(pdef, p, r, s)
        if action == '':
            print('Run with no possible reactions for step = '+str(i))
            for j in range(i, h+1):
                result.append(p)

            #occurrences = add_occurrences(occurrences, 'p', p, (h+1-i))
                
            break
        result.append(p)
        #print('\np',i,': ',p,'\n')

        #occurrences = add_occurrences(occurrences, 'p', p, 1)

        # max level in p
        _maxLevelInP = maxLevelInP(_maxLevelInP, p)     

    """print('\nOCCURRENCES:')
    for occurence in occurrences:
        print('\np:',occurence['p'],'\nOccurrences:',occurence['occurrences'])"""
    
    return result, occurrences, _maxLevelInP

# run: end


"""
pdef: map of all processes
p: initial set of parallel processes
r: list of reaction
s: list of species (with min and max level)
h: number of steps
N: number of run
It call run for N times.
It generates a simulation from p, for h steps.
return: a sequence of h set (each set contains N dataspace),
    maxLevelInPInTheSimulation (maximum simulation, counting each run)
"""
def simulate( pdef, p, r, s, h, N ):

    maxLevelInPInTheSimulation = 0
    print('No. step: '+str(h)+'\nNo. run: '+str(N))
    
    # orders p (set of parallel processes)
    temp_p = [ [] for i in s ] # creates an empty list for each specie
    for species in p:
        index_next_process = get_index(s,species)
        temp_p[index_next_process] = species
    p = temp_p
    
    data = [ [] for i in range(h+1) ] # position 0 contains the initial p (p_0)

    occurrences = []
    
    # for each run, the first element must be p_0 (the initial set of parallel processes)
    for i in range(N):
        data[0].append(p)
    
    for i in range(N):
        print('\n\nRUN: ',i+1)
        sample, temp_occurrences, maxLevelInThisP = run(pdef, p, r, s, h)

        if maxLevelInThisP > maxLevelInPInTheSimulation:
            maxLevelInPInTheSimulation = maxLevelInThisP

        """if i==0:
            occurrences = temp_occurrences
        else:
            for o in temp_occurrences:
                occurrences = add_occurrences(occurrences, 'p', o['p'], o['occurrences'])"""

                
        for j in range(h):
            data[j+1].append(sample[j])
        
    """
    if n != 1:
        print('\n\n\nTOTAL OCCURRENCES:')
        for occurence in occurrences:
            print('\np:',occurence['p'],'\nOccurences:',occurence['occurrences'])
    """
              
    return data, maxLevelInPInTheSimulation

# simulate: end



### DISTANCES

"""
lst1: list of values (for the original system) used to cvalculate the input/output distance
lst2: list of values (for the perturbated system) used to cvalculate the input/output distance
N: number or run for the original system
N*l: number of run for the perturbated system
It calculates the Wasserstein distance as in Castiglioni-Loreti-Tini.
"""
def wasserstein( lst1, lst2, N, l ):
    
    lst1.sort()
    lst2.sort()
    _sum = 0.0
    
    for i in range(N):
        for j in range(l):
            _sum += abs(lst2[i*l+j]-lst1[i])
            
    return _sum/(N*l)

# wasserstein: end


"""
data1: the simulation of a system
data2: the simulation of a system (data1<>data2)
h: number of steps
N: number of run (original behaviour, data1)
N*l: number of run (perturbated behaviour, data2)
rho: penalty function (assign a value used to input or output species
     to calculate the input or output distance)
It calculates the distance between distributions over systems calling wasserstein.
return: max_ (the maximum distance in [0,h]),
    max_interval (the maximum distance in [h/2,h]),
    dist (the list of distances for each step)
"""
def calculate_distance( data1, data2, h, N, l, rho ):
    
    dist = [ 0 for i in range(h+1) ]
    #print('\n\nDISTANCES')
    
    for i in range(h+1): # simulate[0] = p; simulate[1..h+1] = new p
        #print('\n\ni:',i)
        #print('data1[',i,']: ',data1[i])
        lst1 = list(map(rho,data1[i]))
        #print('lst1: ',lst1)
        #print('data2[',i,']: ',data2[i])
        lst2 = list(map(rho,data2[i]))
        #print('lst2: ',lst2)

        # for element in (0..len(list1[0])
        #    wasserstein
        # calculates the weighted sum and inserts it in dist[i]
        
        n_specie = 0
        dist[i] = 0
        while n_specie<len(lst1[0]): #len(lst1[0]) is the total species for the rho
            
            l1 = []
            for el in lst1:
                l1.append(el[n_specie])
                
            l2 = []
            for el in lst2:
                l2.append(el[n_specie])

            dist[i] += wasserstein(l1,l2,N,l)
            
            n_specie += 1 
        
        dist[i] = dist[i]/n_specie
        #print('dist[',i,']:',dist[i])
    
    max_ = 0
    max_interval = 0
    for i in range(h+1):
        if max_ < max(dist[i:]):
            max_ = max(dist[i:])
        if i >= int(h/2):
            if max_interval < max(dist[i:]):
                max_interval = max(dist[i:])
          
    return max_, max_interval, dist

# calculate_distance: end


"""
Used only for debugging.
"""
def distance( pdef, p1, p2, r, s, h, N, l, rho ): 
    
    print('################\nFIRST SIMULATION\n################')
    data1, maxLevelInPInSimulation1 = simulate(pdef, p1, r, s, h, N)
    
    print('\n\n#################\nSECOND SIMULATION\n#################')
    data2, maxLevelInPInSimulation2 = simulate(pdef, p2, r, s, h, N*l)

    max_, max_interval, dist = calculate_distance(data1, data2, h, N, l, rho)
    
    return data1, data2, max_, max_interval, dist, [maxLevelInPInSimulation1, maxLevelInPInSimulation2]

# distance: end


"""
Used only for debugging.
"""
def systems_and_distance( pdef, p1, p2, r, s, h, N, l, rho ): 
    
    data1, data2, max_dist, max_interval, dist, listMaxLevelInP = distance(pdef, p1, p2, r, s, h, N, l, rho)

    print()
    
    data1_max_value = [ 0.0 for i in range(h+1) ]
    for i in range(h+1):
        lst = list(map(rho, data1[i]))
        data1_max_value[i] = max(lst)
    print('System 1 max values:',data1_max_value)

    data2_max_value = [ 0.0 for i in range(h+1) ]
    for i in range(h+1):
        lst = list(map(rho, data2[i]))
        data2_max_value[i] = max(lst)
    print('System 2 max values:',data2_max_value)
    print()
    
    return data1, data1_max_value, data2, data2_max_value, max_dist, dist

# systems_and_distance: end



### PLOT

"""
data, xlabel, ylabel: datas to show the plot
title: the list of titles for each plot
file: name of the file
color_line: color of the line
It shows a plot.
"""
def plot_histogram( data, xlabel, ylabel, title, file, linestyles, color_line ):

    #print('histogram'+str(data))
    #data=[0,2,5,2]
    #print(data)
    values = []
    i = 0
    for value in data:
        values.append([i,value])
        i += 1
    xvalues = [x for (x,y) in values]
    yvalues = [y for (x,y) in values]
    #print(xvalues)
    #print()
    #print(yvalues)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.plot(yvalues, linestyle=linestyles, color=color_line)
    plt.savefig(file+'.png')
    #plt.show()
    plt.close()
    
# plot_histogram: end


"""
datas, legends, xlabel, ylabel: datas to show the plot
title: the list of titles for each plot
file: name of the file
colors: (opt) list of colors.
It shows some plot.
"""
def plot_some_histograms( datas, legends, xlabel, ylabel, title, file, colors ):

    linestyles = ['dashdot', 'dashed', 'solid', 'dotted']
    #linestyles = ['-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']
    
    i = 0
    for data in datas:
        yvalues = np.array(data)
        if len(colors) != 0:
            plt.plot(yvalues,label=legends[i],linestyle=linestyles[i%len(linestyles)], color=colors[i])
        else:
            plt.plot(yvalues,label=legends[i],linestyle=linestyles[i%len(linestyles)])
        i += 1
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    #plt.ylim(0,100) # ylim((bottom, top))   # set the ylim to bottom, top
    plt.savefig(file+'.png')
    #plt.show()
    plt.close()
    
# plot_some_histograms: end



### ROBUSTNESS

"""
p1: initial set of parallel processes (original behaviour)
p2: initial set of parallel processes (perturbated behaviour)
rho_i: penalty function (assign a value used to input species to calculate the input distance)
It calculates the input distance at step 0.
"""
def calculate_distance_i_step0( p1, p2, rho ):

    lst1 = list(map(rho,[p1]))
    #print('lst1:',lst1)

    lst2 = list(map(rho,[p2]))
    #print('lst2:',lst2)

    distance = 0
    i = 0
    while i < len(lst1[0]):
        distance += abs(lst1[0][i]-lst2[0][i])
        i += 1
    distance /= len(lst1[0])
    
    return distance

# calculate_distance_i_step0: end


"""
data1: the simulation of a system
data2: the simulation of a system (data1<>data2)
h: number of steps
n: number of run (of data1)
M*l: number of run (of data2)
rho: penalty function (assign a value used to input species to calculate the input distance)
It calls calculate_distance.
calculate_distance calculates the input distance for each step (which write in dist)
 and the maximum input distance (a).
return: a (the maximum input distance in [0,h]),
    t (step when a occurs),
    dist (the input distance for each step)
"""
def calculate_distance_i( data1, data2, h, N, l, rho ):
    
    max_, max_interval, dist = calculate_distance(data1, data2, h, N, l, rho)
    a = max_
    t = dist.index(a)
    #print('\n\ndist_i:',dist)
    
    return a, t, dist

# calculate_distance_i: end


"""
data1: the simulation of a system
data2: the simulation of a system (data1<>data2)
h: number of steps
n: number of run (of data1)
M*l: number of run (of data2)
rho: penalty function (assign a value used to output species to calculate the output distance)
It calls calculate_distance.
calculate_distance calculates the output distance for each step (which write in dist)
 and the maximum output distance (a).
return: b (the maximum output distance in [0,h]),
    max_interval (the maximum output distance in [h/2, h]),
    t (when b occurs in dist),
    dist (the input distance for each step)
"""
def calculate_distance_o( data1, data2, h, N, l, rho ):
    
    max_, max_interval, dist = calculate_distance(data1, data2, h, N, l, rho)
    b = max_
    t = dist.index(b)

    #print('\n\ndist_o:',dist)
    
    return b, max_interval, t, dist

# calculate_distance_o: end


"""
data1: the simulation of a system
data2: the simulation of a system (data1<>data2)
h: number of steps
N: number of run (of data1)
N*l: number of run (of data2)
rho_i: penalty function (assign a value used to input species to calculate the input distance)
rho_o: penalty function (assign a value used to output species to calculate the output distance)
It calls calculate_distance_i and calculate_distance_o.
calculate_distance_o calculates the output distance for each step (which write in dist)
 and the maximum output distance (a).
calculate_distance_i does the same as calculate_distance_o but for the output.
return: a (the maximum input distance), 
    b (the maximum output distance), 
    max_interval_o (the maximum output distance in [h/2, h]),
    t1 (when a occurs in dist_i),
    t2 (when b occurs in dist_o),
    dist_i (the input distance for each step),
    dist_o (the output distance for each step)
"""
def calculate_robustness( data1, data2, h, N, l, rho_i, rho_o ): 
    
    a, t1, dist_i = calculate_distance_i(data1, data2, h, N, l, rho_i)
    b, max_interval_o, t2, dist_o = calculate_distance_o(data1, data2, h, N, l, rho_o)
    print('Maximum input distance at step #', str(t1), ':', a)
    print('Maximum output distance at step #', str(t2), ':',b)
    print('Maximum output distance in the interval:',max_interval_o)
    
    return a, b, max_interval_o, t1, t2, dist_i, dist_o

# calculate_robustness: end


"""
n: a number (int)
return: the string representing the cardinal number of n
"""
def ordinal( n ):
    
    suffix = {1:'st', 2:'nd', 3:'rd', 11:'th', 12:'th', 13:'th'}
    return str(n)+(suffix.get(n%100) or suffix.get(n%10, 'th'))

# ordinal: end


"""
pdef: map of all processes
other_p: a list of initial set of parallel processes (one element for each perturbated behaviour)
r: list of reaction
s: list of species (with min and max level)
h: number of steps
N: number of run (original behaviour)
N*l: number of run (perturbated behaviour)
rho_i: penalty function (assign a value used to input species to calculate the input distance)
rho_o: penalty function (assign a value used to output species to calculate the output distance)
eta_1: the maximum input distance
It calculates the robustness of the system which has p1 as the initial set of parallel processes.
It calls simulate (to simulate the original and the perturbated behaviours),
    calculate_robustness (to calculate the robustness of the original system
    compared to the perturbated ones),
    check_robustness (to control if the original system is robust).
return: rob_string (a string contains an information about the robustness),
    data0 (the simulation of a perturbated system),
    other_datas (the list of simulation of each perturbated system),
    dist_i_list (the list of input distance for each step),
    dist_o_list (the list of output distance for each step),
    listMaxLevelInP (list of max level in each simulation),
    maximum_input_distance (maximum input distance),
    np.mean(dist_o_list) (mean output distance),
    round(eta_2,num_round) (max output distance in [0,h] rounded to num_round decimal places),
    round(dist_o_interval,num_round) (max output distance in [h/2,h] rounded to num_round decimal places)
"""
def robustness_some_system( pdef, p1, other_p,
                            r, s, h, N, l,
                            rho_i, rho_o, eta_1 ): 

    # SIMULATE
    # original system
    print('\n\n\n##############################\nSIMULATING THE ORIGINAL SYSTEM\n##############################')
    data0, maxLevelInPInSimulation0 = simulate(pdef, p1, r, s, h, N) # our S
    # perturbated systems
    other_datas = []
    other_maxLevelInPInSimulations = []
    print('\n\n\n###############################\nSIMULATING THE PERTURBATED ONES\n###############################')
    no_perturbated_sys = 0
    while no_perturbated_sys < len(other_p):
        print('\n### '+
              ordinal(no_perturbated_sys+1)+
              ' PERTURBATED SYSTEM ###')
        temp_data, temp_maxLevelInPInSimulation = simulate(pdef, other_p[no_perturbated_sys], r, s, h, N*l) # S'_[no_perturbated_sys]
        other_datas.append(temp_data)
        other_maxLevelInPInSimulations.append(temp_maxLevelInPInSimulation)
        no_perturbated_sys += 1

    # COMPARISONS
    print('\n\n\n\n###########\nCOMPARISONS\n###########')
    
    # calculate the distances (dist_i and dist_o)
    data_robustness = []
    dist_i_list = []
    dist_o_list = []

    # for each perturbated system
    maximum_input_distance = 0 # input distance
    rob_string = ''
    no_perturbated_sys = 0
    eta_2 = 0.0
    dist_o_interval = 0.0
    system_with_max_eta_2 = 0
    while no_perturbated_sys < len(other_p):
        print('\n### COMPARE S-S\'_'+str(no_perturbated_sys+1)+' ###')
        a, b, max_interval_o, t1, t2, dist_i, dist_o = calculate_robustness(data0, other_datas[no_perturbated_sys], h, N, l, rho_i, rho_o)
        if a > maximum_input_distance:
            maximum_input_distance = a
        dist_i_list.append(dist_i)
        dist_o_list.append(dist_o)
        print('S - S\'_'+str(no_perturbated_sys+1)+
              ' is (d_i,d_o,[0,0],['+str(int(h/2))+','+str(h)+'],eta_1='+
              str(round(a,num_round))+',eta_2='+str(round(max_interval_o,num_round))+')-robust')
        the_robustness = 'S - S\'_'+str(no_perturbated_sys+1)+' is '
        the_robustness += '(d_i,d_o,[0,0],['+str(int(h/2))+','+str(h)+'],'+str(round(a,num_round))+','+str(round(max_interval_o,num_round))+')-robust'
        the_robustness += '\nThe maximum input distance is '+str(a)
        the_robustness += ' and it is observed at step #'+str(t1)
        if t1 != 0:
            the_robustness += '\nThe maximum input distance at step #0 is '+str(dist_i[0])
        the_robustness += '\nThe maximum output distance is '+str(b)
        the_robustness += ' and it is observed at step #'+str(t2)
        the_robustness += '\nThe maximum output distance in the interval '
        the_robustness += '['+str(int(h/2))+','+str(h)+'] is '+str(max_interval_o)
        rob_string += the_robustness + '\n'
        
        if eta_2 < b:
            eta_2 = b
        if dist_o_interval < max_interval_o:
            dist_o_interval = max_interval_o
            system_with_max_eta_2 = no_perturbated_sys + 1 # no_perturbated_sys=0..(no. sys-1)
        
        print('\n\n')
        no_perturbated_sys += 1

    _str = 'The system is (d_i,d_o,[0,0],['+str(int(h/2))+','+str(h)+'],eta_1='+str(round(eta_1[1],num_round))+',eta_2='+str(round(dist_o_interval,num_round))+')-robust'
    _str += '\nThe maximum output distance in the interval '
    _str += '[0,'+str(h)+'] is '+str(round(eta_2,num_round))
    _str += '\nThe maximum output distance in the interval '
    _str += '['+str(int(h/2))+','+str(h)+'] is '+str(round(dist_o_interval,num_round))
    _str += '\n\nThe system with the maximum output distance in ['+str(int(h/2))+','+str(h)+'] is S\'_'+str(system_with_max_eta_2)
    print('\n\n')
    print(_str)
    rob_string += '\n\n\n' + _str

    listMaxLevelInP = [maxLevelInPInSimulation0]+other_maxLevelInPInSimulations
    
    return rob_string, data0, other_datas, dist_i_list, dist_o_list, listMaxLevelInP, maximum_input_distance, np.mean(dist_o_list), round(eta_2,num_round), round(dist_o_interval,num_round)

# robustness_some_system: end


"""
data: the simulation of a system
from_level..to_level: steps to see the concentration
get_output_level: function which, given the parallel processes, returns the output level
It returns the concentration of output species in the interval [from_level,to_level]
"""
def get_concentration_system(data, from_level, to_level, get_output_level):

    count = from_level
    lists = []
    
    while count <= to_level:
        lst = list(map(get_output_level,data[count]))
        # e.g. lst: [[1], [1], [1], [1], [1]] (for one step, when h=2 and run=5 and 1! specie)
        #print('\n\ncount: '+str(count))
        #print(lst)
        lists += lst
        count += 1
        
    #print('\n\nConcentrations temp:',lists)
    # e.g. [[1], [1], [1], [2], [2], [1], [1], [1], [2], [2]]
    #print(lists == [[2], [1], [1], [1], [1], [3], [1], [2], [1], [2]])
    
    #lists = [[1,2], [1,2], [1,2], [1,2], [1,2]]
    # e.g. [[1], [1], [1], [1], [1], [2], [2], [2], [2], [1]]
    # for all steps (1..2), when h=2 (steps: 0..2) and run=5 and only 1 output specie
    no_output_species = len(lists[0]) # num of output species
    i = 0
    concentration = []
    while i < no_output_species:
        concentration_temp = [] # concentration of the i° output specie
        for levels in lists:
            """print('len:',len(levels))
            print('value:',levels)
            print('list:',list(levels))
            levels=list(levels)
            print('levels[0]:',levels[0])"""
            # levels: [level_1st_output_specie, ..., level_nth_output_specie]
            concentration_temp.append(levels[i])
        concentration.append(concentration_temp)
        i += 1
        
    # e.g. concentration: [[1, 1, 1, 2, 2, 1, 1, 1, 2, 2]]
    # because there is only one output specie
    return concentration

# get_concentration_system: end
