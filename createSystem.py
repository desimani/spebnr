# (createSystem.py) last update: 2023.05.22


from spebnr import *
import random
import os # directories
import scipy.special as spe # binomial coefficient



num_round = 6 # useful to round the distance
max_no_perturbated_sys = 100


"""
strings: list of string
file: name of the file
It writes something (strings) in a file file.txt
"""
def write( strings, file ):

    file = file+'.txt'
    f = open(file, 'w')  # write in text mode
    for string in strings:
        f.write(str(string))
    f.close()

    
"""
SPECIES: list of [species_name, min_level, max_level, roles for each reaction] for each specie
reactions: all possible reactions (list of reactions and their own rate)
It creates the map consisting of all possible processes.
return: map_processes (the map of all possible processes)
"""
def create_map( SPECIES, reactions ):
        
    reactions_name = []
    for i in range(0, len(reactions)):
        reactions_name.append(list(reactions[i].keys())[0])
 
    map_processes = []

    # for each specie
    for s in SPECIES:
        name = s[0]
        min_level = s[1]
        max_level = s[2]
        reactions_roles = s[3]
    
        # for each level (min..max), create a process and insert it in the map
        for i in range(min_level, 1+max_level): # range(a,b): a incluso, b escluso: a <= i < b
            process_name = name + '_' + str(i) # process_name = name_level

            set_reactions = []
        
            # read each reaction
            # write each reaction in a set
            # insert specie+reaction in the map
            for j in range(0, len(reactions_roles)):
                role = reactions_roles[j]
                num_molecules_r = role[0]
                num_molecules_p = role[1]

                condition1 = (num_molecules_r > 0) and (i-num_molecules_r < min_level)
                # condition1 is True if the role is 'reactant'
                # AND the specie has not enought level (i)

                condition2 = (num_molecules_p > 0) and (i-num_molecules_r+num_molecules_p > max_level)
                # condition2 is True if the role is 'product'
                # AND the specie will have too level (i)
            
                if not (condition1 or condition2):
                    weight = 1

                    # next process
                    _next = name + '_'
                    nextLevel = i
                    if num_molecules_r > 0: # reactant 
                        #weight = i ### ATTENTION: you can use it when all [r,p] have a number of r<=1 (e.g. [1,0])
                        weight = spe.binom(i,num_molecules_r) 
                        nextLevel -= num_molecules_r
                    if num_molecules_p > 0: # product
                        nextLevel += num_molecules_p
                    # a species can participates in a reaction as a reactant AND as a product
                    _next += str(nextLevel)

                    # reaction name
                    reaction_id = reactions_name[j]
               
                    set_reactions.append({
                        'reaction_id':reaction_id, #id: name
                        'role':role,
                        'weight':weight,
                        'next_process':_next,
                    })
        
            dictionary = {'process_name':process_name,
                'reactions':set_reactions, # reaction for the process process_name
            }
            map_processes.append(dictionary)
    
    return map_processes

# create_map: end
            

"""
species: list of [species_name, min_level, max_level, roles for each reaction] for each specie
It returns the list of [name, no_min_level, no_max_level] for each specie
return: result (list of [name, no_min_level, no_max_level] for each specie)
"""
def get_species_levels( species ):

    result = []
    for s in species:
        species_levels = []
        species_levels.append(s[0]) # name
        species_levels.append(s[1]) # no_min_level
        species_levels.append(s[2]) # no_max_level
        result.append(species_levels)
        
    return result

# get_species_levels: end


"""
species: list of [species_name, min_level, max_level, roles for each reaction] for each specie
It returns the minimum and the maximum level.
return: no_min_level (the minimal level of all species),
    no_max_level (the maximal level of all species)
"""
def find_min_max_levels( species ):

    no_min_level = species[0][1]
    no_max_level = species[0][2]
    for s in species:
        if s[1] < no_min_level:
            no_min_level = s[1] # no_min_level
        if s[2] > no_max_level:
            no_max_level = s[2] # no_max_level
            
    return no_min_level, no_max_level

# find_min_max_levels: end


"""
species: list of [species_name, min_level, max_level, roles for each reaction] for each specie
specie_name: the name of a species
It returns the minimum and the maximum level of a specie.
return: no_min_level (min level for specie_name), 
    no_max_level (max level for specie_name)
"""
def find_min_max_levels_of_a_specie( species, specie_name ):

    no_min_level = species[0][1]
    no_max_level = species[0][2]
    for s in species:
        if s[0] == specie_name:
            no_min_level = s[1] # no_min_level
            no_max_level = s[2] # no_min_level
            return no_min_level, no_max_level
        
# find_min_max_levels_of_a_specie: end


"""
It prints something only for show.
"""
def print_some_info( h, N, l, no_min_level, no_max_level,
                     p, other_p, eta_1,
                     example ):
    
    print('###################')
    print(example)
    print('###################')
    print('PARALLEL PROCESSES:')
    print('###################')
    print('Original System:', p)
    print('###################')
    count = 0
    print('Perturbated Systems:')
    while count < len(other_p):
        print('### p of S\'_'+str(count+1)+':', other_p[count])
        count += 1
    print('###################')
    print('Levels: from', no_min_level, 'to', no_max_level)
    print('###################')
    print('h:', h, 'N:', N, 'l:', l)
    print('###################')
    print('eta_1: ('+str(eta_1[0])+','+str(eta_1[1])+']')
    print('###################')
    
# print_some_info: end


"""
SPECIES: list of [species_name, min_level, max_level, roles for each reaction] for each specie
reactions: all possible reactions (list of reactions and their own rate)
p: original process
no_perturbated_sys: number of perturbated processes
INPUT_SPECIES: list of input species
OUTPUT_SPECIES: list of output species
h: number of steps
N: number of run for the original processes
N*l: number of run for each perturbated processes
eta_1: referred to the input threshold
example: list of string of the name of example
It checks the input.
return: False if there is at least one error; True otherwise.
"""
def check_input( SPECIES, reactions, p, no_perturbated_sys,
                 INPUT_SPECIES, OUTPUT_SPECIES,
                 h, N, l, eta_1, 
                 example ): 
    
    # h, N, l
    if type(h) != int:
        return False
    if h != int(h):
        return False
    if type(N) != int:
        return False
    if N != int(N):
        return False
    if type(l) != int:
        return False
    if l != int(l):
        return False
    #print('ok hnl')
    
    # at least two species
    if type(SPECIES) != list:
        return False
    if len(SPECIES) <= 1:
        return False
    #print('ok species')

    # eta_1
    if type(eta_1) != float: # then, (min, max] in a list
        return False
    if eta_1 > 1 or eta_1 < 0:
        return False
    #print('ok eta_1')

    # at least two reactions
    if type(reactions) != list:
        return False
    if len(reactions) <= 1:
        return False
    for r in reactions:
        value = r[list(r.keys())[0]]
        if not(type(value) != int or type(value) != float):
            return False
        if value < 0:
            return False
    #print('ok reactions')

    # input species
    name_of_species = set()
    for s in SPECIES:
        name_of_species.add(s[0]) # we'd better use the same name (not .upper() or .lower())
    for in_sp in INPUT_SPECIES:
        if in_sp not in name_of_species:
            return False
    #print('ok input species')

    # output species
    for out_sp in OUTPUT_SPECIES:
        if out_sp not in name_of_species:
            #print('out_sp: ',out_sp)
            return False
    #print('ok output species')
      
    # for each specie
    for s in SPECIES:
        #print(s,' ',)
        # levels (min<max)
        if not((type(s[1]) == int) and (type(s[2]) == int)):
            return False
        if s[1] > s[2]: # min can be equal to max
            return False
        # reactions
        if len(s[3]) != len(reactions):
            return False
        # [R,P]
        for role in s[3]:
            if not((type(role[0]) == int) and (type(role[1]) == int)):
                return False
            if role[0]<0 or role[1]<0:
                return False
    #print('ok each specie')

    # number of perturbated systems
    if type(no_perturbated_sys) != int:
        return False
    if no_perturbated_sys < 1 or no_perturbated_sys > max_no_perturbated_sys:
        return False
    #print('ok no. perturbated system')

    # example
    if type(example) != list:
        return False
    for ex in example:
        if type(ex) != str:
            return False
    #print('ok example name')

    # for each p
    try:
        for _p in p:
            
            #print(_p)
            if type(_p) != list:
                return False
                
            for sp_level in _p:
                    
                # specie name
                if type(sp_level) != str:
                    return False
                _split = sp_level.split('_')
                if len(_split) < 2:
                    return False
                    
                # species level: integer and in [min, max]
                index = get_index( SPECIES, sp_level)
                if index == -1:
                    return False
                level = _split[len(_split)-1]
                level = int(level)
                level_min = SPECIES[index][1]
                level_max = SPECIES[index][2]
                if level<level_min or level>level_max:
                    return False
        #print('ok each p')
   
    except:
        return False
    
    return True

# check_input: end


"""
s: species
original_p: original process
INPUT_SPECIES: list of input species
rho: penalty function
eta_1: referred to the input threshold
It generates a random process.
return: process (a random perturbated process)
"""
def generate_a_perturbated_process( s, original_p, INPUT_SPECIES, rho, eta_1 ):

    process = original_p.copy()

    while True:
        for specie_name in INPUT_SPECIES:
            min_lev, max_lev = find_min_max_levels_of_a_specie(s, specie_name)
            random_level = random.randrange(min_lev, max_lev+1) # min_lev <= random_level < max_lev+1
            specie_level = specie_name+'_'+str(random_level)
            index = get_index(s,specie_name) # in spear.py
            process[index] = specie_level
        distance = calculate_distance_i_step0(original_p, process, rho)
        if distance > eta_1[0] and distance <= eta_1[1]:
            return process
        
# generate_a_perturbated_process: end
    

"""
s: species
p: original process
no_perturbated_sys: number of perturbated process to generate
INPUT_SPECIES: list of input species
rho: penalty function
eta_1: referred to the input threshold
It generates the processes of (no_perturbated_sys) perturbated systems.
return: processes (a list of perturbated processes)
"""
def generate_some_perturbated_processes( s, p, no_perturbated_sys, INPUT_SPECIES, rho, eta_1 ):

    processes = []
    for n in range(1,no_perturbated_sys+1):
        processes.append(generate_a_perturbated_process(s, p, INPUT_SPECIES, rho, eta_1))

    return processes

# generate_some_perturbated_processes: end
    

"""
Writes lots of information in files and generates plots.
"""
def write_all_info_in_files_and_plot(h, N, l,
                            example,
                            INPUT_SPECIES, OUTPUT_SPECIES,
                            SPECIES, r, no_perturbated_sys,
                            get_input_level, get_output_level,
                            p, other_p,
                            data0, other_datas,
                            eta_1, eta_2,
                            dist_i_list, dist_o_list, dist_o_interval,
                            maximum_input_distance, mean_output_distance,
                            listMaxLevelInP, 
                            rob_string): 

    # create directories
    # main directory
    main_txt = example[0]
    main_directory = main_txt
    mdir = ''
    pdir = ''
    cdir = ''
    sdir = ''
    count = 0
    while True:
        if count != 0:
            main_directory = main_txt + '('+str(count)+')'
        mdir = os.path.join('',main_directory)
        if not os.path.exists(main_directory):
            os.mkdir(mdir)
            # plots
            pdir = os.path.join(main_directory,'plots')
            os.mkdir(pdir)
            # concentrations (only if comparing)
            if len(example) == 1:
                cdir = os.path.join(main_directory,'concentrations')
                os.mkdir(cdir)
            sdir = os.path.join(main_directory,'systems')
            os.mkdir(sdir)
            print('\n\n\nThe directories '+str(mdir)+', '+str(pdir)+' and '+str(sdir)+' have been created\n\n\n')
            break
        count += 1

    # write all information in robustness.txt
    title = 'robustness'
    list_to_write = ['p of S:',str(p)]
    count = 0
    while count < len(other_p):
        list_to_write += ['\np of S\'_'+str(count+1)+':'+str(other_p[count])]
        count += 1
    list_to_write += ['\n\nS:\n',str(data0)]
    count = 0
    while count < len(other_datas):
        list_to_write += ['\nS\'_'+str(count+1)+':'+str(other_datas[count])]
        count += 1
    count = 0
    while count < len(dist_i_list):
        list_to_write += ['\ndist_i S-S\'_'+str(count+1)+':'+str(dist_i_list[count])]
        list_to_write += ['\ndist_o S-S\'_'+str(count+1)+':'+str(dist_o_list[count])]
        count += 1
    list_to_write += ['\nMaximum input distance: ' + str(maximum_input_distance)]
    list_to_write += ['\nMean output distance: ' + str(mean_output_distance)]
    list_to_write += ['\n\n']
    list_to_write += ['\n\n',rob_string]
    write( list_to_write , os.path.join(mdir,title))

    # write the maximum input distance in maximum_input_distance.txt
    write(['The maximum input distance is ',str(maximum_input_distance),'\n',str(round(maximum_input_distance,num_round))],
          os.path.join(mdir,'maximum_input_distance'))

    # write the mean output distance in mean_output_distance.txt
    write(['The mean output distance is ',str(mean_output_distance),'\n',str(round(mean_output_distance,num_round))],
          os.path.join(mdir,'mean_output_distance'))

    # write results in results.txt
    write([rob_string], os.path.join(mdir,'results'))

    # write all p in each_p.txt
    list_to_write = ['p of S: '+str(p)]
    count = 0
    while count < len(other_datas):
        list_to_write += ['\np of S\'_'+str(count+1)+': '+str(other_p[count])]
        count += 1
    write(list_to_write, os.path.join(mdir,'each_p'))

    # write data0 in S.txt
    write(['S:\n',str(data0)], os.path.join(sdir,'S'))
    
    # write datax in S'_x.txt
    count = 0
    while count < len(other_datas):
        list_to_write += ['\nS\'_'+str(count+1)+':'+str(other_datas[count])]
        write(['S\'_'+str(count+1)+':\n',str(other_datas[count])], os.path.join(sdir,"S'_"+str(count+1)))
        count += 1
    
    # write the maximum level of species in p, in a file
    list_to_write = ['Maximum level in S: '+str(listMaxLevelInP[0])]
    count = 1
    while count < len(listMaxLevelInP):
        list_to_write += ['\nMaximum level in S\'_'+str(count)+': '+str(listMaxLevelInP[count])]
        count += 1
    write(list_to_write, os.path.join(mdir,'listMaxLevelInP'))

    # write dist_i in dist_i.txt
    write(str(dist_i_list), os.path.join(mdir,'dist_i'))

    # write dist_o in dist_o.txt
    write(str(dist_o_list), os.path.join(mdir,'dist_o'))
 
    # print only the most interesting info in readme.txt
    title = 'readme'
    list_to_write = ['Example: ', example[0]]
    list_to_write += ['\nh: ',h]
    list_to_write += ['\nN: ',N]
    list_to_write += ['\nl: ',l]
    list_to_write += ['\nNo. perturbated systems: ',no_perturbated_sys]
    list_to_write += ['\nNo. reactions: ',len(r)]
    list_to_write += ['\nNo. species: ',len(SPECIES)]
    list_to_write += ['\nNo. input species: ',len(INPUT_SPECIES)]
    list_to_write += ['\nNo. output species: ',len(OUTPUT_SPECIES)]
    list_to_write += ['\neta_1 (maximum input distance): ',str(round(maximum_input_distance,num_round))]
    list_to_write += ['\neta_2 (output distance in the interval [',int(h/2),',',h,']): ',round(dist_o_interval,num_round)]
    list_to_write += ['\nOutput distance in the interval [0,',h,']: ',str(round(eta_2,num_round))]
    list_to_write += ['\nThe system is (d_i,d_o,[0,0],['+str(int(h/2))+','+str(h)+'],eta_1='+
              str(round(eta_1[1],num_round))+',eta_2='+str(round(dist_o_interval,num_round))+')-robust']
    write( list_to_write , os.path.join(mdir,title))

    # plot: colors
    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'cyan']
    
    # plot distances
    legends = ['dist_i','dist_o']
    xlabel = 'step'
    ylabel = 'rho'
    count = 0
    while count < len(dist_i_list):
        title = 'robustness_S-S\'_'+str(count+1)
        plot_some_histograms([dist_i_list[count], dist_o_list[count]],
                             legends, xlabel, ylabel, title, os.path.join(pdir,title), colors=[])
        count += 1

    # plot concentration iff compare original-perturbated
    if len(example) > 1:
        # plot output (original)
        xlabel = 'step'
        ylabel = 'concentration'
        title = example[1]
        path = os.path.join(pdir,title)
        plot_concentrations(data0, INPUT_SPECIES, OUTPUT_SPECIES,
                        get_input_level, get_output_level,
                        xlabel, ylabel,
                        title,
                        path)
    
        # plot output (perturbated)
        xlabel = 'step'
        ylabel = 'concentration'
        title = example[2]
        path = os.path.join(pdir,title)
        plot_concentrations(other_datas[0], INPUT_SPECIES, OUTPUT_SPECIES,
                        get_input_level, get_output_level,
                        xlabel, ylabel,
                        title,
                        path)
    
# write_all_info_in_files_and_plot: end
    

"""
SPECIES: list of [species_name, min_level, max_level, roles for each reaction] for each specie
r: all possible reactions (list of reactions and their own rate)
p: original process
no_perturbated_sys: number of perturbated processes
INPUT_SPECIES: list of input species
OUTPUT_SPECIES: list of output species
h: number of steps
N: number of run for the original processes
N*l: number of run for each perturbated processes
rank_input: a function which returns a value used to calculate the output distance
rank_output: a function which returns a value used to calculate the input distance
get_output_level: a function which returns a list of the maximum output concentration
get_output_max_level: a function which returns a list of the output species name and the maximum output concentration
eta_1: referred to the input threshold
example: list of string of the name of example
It calls generate_some_perturbated_processes to generate () random processes having an input distance in (eta_1-1,eta_1].
    Then, it calls robustness_some_system to calculates the robustness of the original system.
    Finally, it calls write_all_info_in_files_and_plot to print information and generates plots.
"""
def execute_robustness( SPECIES, r, p, no_perturbated_sys,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_output_level, get_output_max_level,
                        eta_1,
                        example):
    
    if not (check_input(SPECIES, r, [p], no_perturbated_sys,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        eta_1, 
                        example)):
        print('InputError!')
        raise Exception("Exception!")

    eta_1 = [eta_1-1, eta_1] #interval eta_1
    
    pdef = create_map(SPECIES, r) #r: reactions
    
    s = get_species_levels(SPECIES)
    no_min_level, no_max_level = find_min_max_levels(SPECIES)
    other_p = generate_some_perturbated_processes(s, p, no_perturbated_sys, INPUT_SPECIES, rank_input, eta_1)
    
    print_some_info(h, N, l, no_min_level, no_max_level,
                    p, other_p, eta_1, example[0])
    print('ROBUSTNESS\n###################')

    # execute robustness
    rob_string, data0, other_datas, dist_i_list, dist_o_list, listMaxLevelInP, maximum_input_distance, mean_output_distance, eta_2, dist_o_interval = robustness_some_system(pdef, p, other_p, r, s, h, N, l, rank_input, rank_output, eta_1)

    write_all_info_in_files_and_plot(h, N, l,
                            example,
                            INPUT_SPECIES, OUTPUT_SPECIES,
                            SPECIES, r, no_perturbated_sys,
                            [], get_output_level, 
                            p, other_p,
                            data0, other_datas,
                            eta_1, eta_2,
                            dist_i_list, dist_o_list, dist_o_interval,
                            maximum_input_distance, mean_output_distance,
                            listMaxLevelInP, 
                            rob_string)
    
# execute_robustness: end


"""
s: species [name, min_level, max_level (approximated)]
data: all simulation
It finds the minimum and the maximum level of each species.
return: string_to_write (a string used to print relevant information)
"""
def find_min_max_level_each_species(s, data): 

    print('\nCALCULATING MIN e MAX LEVELS\n')
    print('data:')
    print(data)

    name_species_list = []
    for specie in s:
        name_species_list.append(specie[0]) # name
    print('name_species_list:')
    print(name_species_list)

    # creates a list with the maximum level for each specie
    max_levels_list = [ 0 for i in s ] 
    min_levels_list = [ 5000 for i in s ] 

    # find all process in p
    for list_processes in data:
        for p in list_processes:
            #print('p is: ',str(p))
            i = 0
            while i<len(p):
                specie_level = p[i] # e.g. 'X_1'
                #print('i: ',str(i),' -- specie_level: ',str(specie_level))
                split = specie_level.split('_')
                level = int(split[len(split)-1])
                if level > max_levels_list[i]:
                    max_levels_list[i] = level
                if level < min_levels_list[i]:
                    min_levels_list[i] = level
                i += 1
            #print('now max_levels_list is: ',max_levels_list)

    # write in a string
    string_to_write = 'Max:\n'
    for i in range(len(s)):
        string_to_write += str(name_species_list[i]) + ': ' + str(max_levels_list[i]) + '\n'
    string_to_write += '\n\nMin:\n'
    for i in range(len(s)):
        string_to_write += str(name_species_list[i]) + ': ' + str(min_levels_list[i]) + '\n'
    
    
    return string_to_write

# find_min_max_level_each_species: end

    
"""
ATTENTION: relevant species must be <= 6
It generates plots (input + output).
"""
def plot_concentrations(data, INPUT_SPECIES, OUTPUT_SPECIES,
                        get_input_level, get_output_level,
                        xlabel, ylabel,
                        title,
                        path): 

    # colors for plot
    colors = ['red', 'blue', 'green', 'orange', 'magenta', 'cyan'] ### ATTENTION: relevant species must be <= 6
    linestyles = ['dashdot', 'dashed', 'solid', 'dotted']

    # input + output
    tot_relevant_species = len(INPUT_SPECIES) + len(OUTPUT_SPECIES)
    concentrations = [ [] for i in range(tot_relevant_species) ] # e.g. [ [mean X each step], [mean Y each step], [mean YP each step] ]
    step = 0
    i = 0
    #print('concentrations:',concentrations)
    while i<len(data):
        #print('concentrations:',concentrations)
        conc_i = get_concentration_system(data, i, i, get_input_level) # e.g. [ [concentrations X], [concentrations Y] ]
        count = 0
        for conc in conc_i: # for each input species
            #print('conc_i:',conc_i)
            concentrations[count].append(mean(conc_i[count]))
            count += 1
        conc_o = get_concentration_system(data, i, i, get_output_level) # e.g. [ [concentrations YP] ]
        count_o = 0
        for conc in conc_o: # for each output species
            concentrations[count].append(mean(conc_o[count_o]))
            count += 1
            count_o += 1
        i += 1

    # plot
    legends = INPUT_SPECIES + OUTPUT_SPECIES
    plot_some_histograms(concentrations, legends,
                         xlabel, ylabel, title, path, colors[:len(legends)])

    # one plot for each specie
    i = 0
    while i < len(legends):
        plot_histogram( concentrations[i], 'step', 'concentration '+str(legends[i]),
                        title+' '+str(legends[i]), path+str(legends[i]), linestyles[i%len(linestyles)], colors[i] )
        i += 1
        
# plot_concentrations: end


"""
SPECIES: list of [species_name, min_level, max_level, roles for each reaction] for each specie
r: all possible reactions (list of reactions and their own rate)
p_original: original process
p_perturbated: perturbated process
INPUT_SPECIES: list of input species
OUTPUT_SPECIES: list of output species
h: number of steps
N: number of run for the original processes
N*l: number of run for the perturbated processes
rank_input: a function which returns a value used to calculate the output distance
rank_output: a function which returns a value used to calculate the input distance
get_output_level: a function which returns a list of the maximum output concentration
get_output_max_level: a function which returns a list of the output species name and the maximum output concentration
example: list of string of the name of example
Used to compare two systems: the original and the perturbated one. And to compare our model with Nasti-Gori-Milazzo one.
It calls calculate_distance_i_step0 to calculate the input distance between the original system and the perturbated one.
    Then, it calls robustness_some_system to calculates the robustness.
    Finally, it calls write_all_info_in_files_and_plot to print information and generates plots.
"""
def compare_two_systems( SPECIES, r, p_original, p_perturbated,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        example): 

    # only for debugging
    no_perturbated_sys = 1
    eta_1 = 0.3
    other_p = [p_perturbated]
    
    if not (check_input(SPECIES, r, [p_original, p_perturbated],
                        no_perturbated_sys,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        eta_1, 
                        example)):
        print('InputError!')
        raise Exception("Exception!")
    
    pdef = create_map(SPECIES, r) #r: reactions
    
    s = get_species_levels(SPECIES)
    no_min_level, no_max_level = find_min_max_levels(SPECIES)
    other_p = [p_perturbated]

    print('###################')
    print(example[0])
    print('###################')
    print('PARALLEL PROCESSES:')
    print('###################')
    print('Original System: '+str(p_original))
    print('###################')
    print('Perturbated System: '+str(p_perturbated))
    print('###################')
    print('h:', h, 'N:', N, 'l:', l)
    print('###################')
    print('ROBUSTNESS')
    print('###################')
    
    # eta_1
    distance_i = calculate_distance_i_step0(p_original, p_perturbated, rank_input)
    eta_1 = [0, distance_i]
    print('eta_1',round(distance_i,num_round))
    print('###################')
    
    # execute robustness
    rob_string, data0, other_datas, dist_i_list, dist_o_list, listMaxLevelInP, maximum_input_distance, mean_output_distance, eta_2, dist_o_interval = robustness_some_system(pdef, p_original, other_p, r, s, h, N, l, rank_input, rank_output, eta_1)

    write_all_info_in_files_and_plot(h, N, l,
                            example,
                            INPUT_SPECIES, OUTPUT_SPECIES,
                            SPECIES, r, no_perturbated_sys,
                            get_input_level, get_output_level,
                            p_original, other_p,
                            data0, other_datas,
                            eta_1, eta_2, 
                            dist_i_list, dist_o_list, dist_o_interval,
                            maximum_input_distance, mean_output_distance,
                            listMaxLevelInP,
                            rob_string)
    
# compare_two_systems: end


"""
Only for debugging
"""
def find_max_level( SPECIES, r, p, no_perturbated_sys,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output, 
                        example): 

    eta_1 = 0.7 # only for debugging
    l = 10 # only for debugging
    
    if not (check_input(SPECIES, r, [p], no_perturbated_sys,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        eta_1, 
                        example)):
        print('InputError!')
        raise Exception("Exception!")
    
    pdef = create_map(SPECIES, r) #r: reactions
    s = get_species_levels(SPECIES)
    no_min_level, no_max_level = find_min_max_levels(SPECIES)

    
    print('SIMULATION')

    # calculate simulate
    data, maxLevelInPInTheSimulation = simulate(pdef, p, r, s, h, N)
    
    # create directories
    # main directory
    main_txt = 'Max level h'+str(h)+' N'+str(N)
    main_directory = main_txt
    mdir = ''
    count = 0
    while True:
        if count != 0:
            main_directory = main_txt + ' ('+str(count)+')'
        mdir = os.path.join('',main_directory)
        if not os.path.exists(main_directory):
            os.mkdir(mdir)
            print('\n\n\nThe directory '+str(mdir)+' has been created\n\n\n')
            break
        count += 1

    # write the minimum and maximum level of each species in minMaxLevels.txt
    levels_to_write = find_min_max_level_each_species(s, data)
    print('Minimum and maximum level of each species:\n\n'+levels_to_write)
    list_to_write = 'Maximum level of each species:\n\n' + levels_to_write
    write(list_to_write, os.path.join(mdir,'minMaxLevels'))

    # write the minimum and maximum level of species in p, in a file
    print('maxLevelInPInTheSimulation: ',maxLevelInPInTheSimulation,'\n')
    list_to_write = ['Maximum level in S: ', str(maxLevelInPInTheSimulation)]
    write(list_to_write, os.path.join(mdir,'maxLevelInP'))

    # write the maximum level in input in a file
    list_to_write = 'Possible maximum level: ' + str(SPECIES[0][2])
    print(list_to_write)
    write(list_to_write, os.path.join(mdir,'possibleMaxLevel'))

    # write data0 in S.txt
    write(['S:\n',str(data)], os.path.join(mdir,'S'))
    
    # plot output
    xlabel = 'step'
    ylabel = 'concentration'
    title = example[0]
    path = os.path.join(mdir,title)
    plot_concentrations(data, INPUT_SPECIES, OUTPUT_SPECIES,
                        get_input_level, get_output_level,
                        xlabel, ylabel,
                        title,
                        path)
    
# find_max_level: end
