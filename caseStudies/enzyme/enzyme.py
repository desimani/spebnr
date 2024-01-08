# (enzyme.py) last update: 2023.07.14


from createSystem import *



example = 'Enzyme'


# reactions: a list of object {reaction name + reaction rates (k)}
reactions = [
    {'a1': 100},
    {'a2': 10},
    {'a3': 0.5},
    {'a4': 0.01},
    {'a5': 0.5}
    ]



#LEVELS: no_min_level...no_max_level
no_min_level = 0
no_max_level = int(15000*1.1)

SPECIES = [ # [species_name, min_level, max_level, roles for each reaction]
        # roles [no. molecules as R, no. molecules as P]
        # R --> [1,0]
        # P --> [0,1]
        # neither R nor P --> [0,0]
    
    ['R',no_min_level,int(1001*1.1),([1,0],[0,0],[0,1],[0,0],[0,0])],
    
    ['X',no_min_level,int(1032*1.1),([1,2],[1,0],[0,0],[1,0],[0,0])],
    
    ['Z',no_min_level,int(1000*1.1),([0,1],[0,0],[1,0],[0,0],[0,0])],
    
    ['P',no_min_level,int(20849*1.1),([0,0],[0,0],[0,0],[1,2],[1,0])],
    
    ['C',1,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,1])],
    
    ['W',10,no_max_level,([0,0],[0,1],[0,0],[0,0],[0,0])],
    
    ]



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the input distance
"""
# input: P
index_input = 3 # always
max_level_input = SPECIES[index_input][2]
def rank_input(p):
    
    # p is always ordered
    
    _split = p[index_input].split('_')
    level_input = int(_split[len(_split)-1])
    
    return [ level_input / max_level_input ] 



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the output distance
"""
# output: X
index_output = 1 # always
max_level_output = SPECIES[index_output][2]
def rank_output(p):
    
    # p is always ordered
    
    _split = p[index_output].split('_')
    level_output = int(_split[len(_split)-1])
    return [ level_output / max_level_output ] 



"""
p: the set of parallel processes in a specific time
return: the input concentration
"""
def get_input_level(p):
    
    input_split = p[index_input].split('_')
    level_input = int(input_split[len(input_split)-1])
    
    return [level_input] # 1 size because there is only one output species



"""
p: the set of parallel processes in a specific time
return: the output concentration
"""
def get_output_level(p):
    
    _split = p[index_output].split('_')
    level_output = int(_split[len(_split)-1])
    
    return [level_output] # 1 size because there is only one output species



"""
return: a list of the output species name and the maximum output concentration
"""
def get_output_max_level():
    
    return [max_level_output]





# input: P
# output: X
INPUT_SPECIES = ['P']
OUTPUT_SPECIES = ['X']

# initial set of parallel processes
# there is not the list of reactions
# spear.py reads the list of reaction via the map

p_original = ['R_1000','X_30','Z_0','P_1','C_1','W_10'] # P1

p_perturbated1 = ['R_1000','X_30','Z_0','P_1000','C_1','W_10'] # P1000

p_perturbated2 = ['R_1000','X_30','Z_0','P_20000','C_1','W_10'] # P20000





### FIND MAX LEVELS
"""
# original system
N=30
h=5000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\nOriginal system')
h=10000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=15000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=20000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=30000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])


print('\n\n\n1st perturbed system\n')



# 1st perturbed system
N=30
h=5000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=10000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=15000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=20000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=30000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])

                        
print('\n\n\n2nd perturbed system\n')

# 2nd perturbed system
N=30
h=5000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=10000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=15000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=20000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
print('\n\n\n')
h=30000
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
"""





### COMPARE ORIGINAL - PERTURBED
"""
h=30000
N=30
l=2

# perturbated 1
example = 'EnzymeCompare1st'
example0 = 'EnzymeP1'
example1 = 'EnzymeP1000'
compare_two_systems( SPECIES, reactions, p_original, p_perturbated1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example1])

# perturbated 2
example = 'EnzymeCompare2nd'
example2 = 'EnzymeP20000'
compare_two_systems( SPECIES, reactions, p_original, p_perturbated2,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example2])
"""





### ESTIMATING ROBUSTNESS

h = 30000
N = 20
l = 2
no_perturbated_systems = 10
p1 = p_original

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.3)
example = 'Enzyme03'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.3,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.4)
example = 'Enzyme04'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.4,
        [example])
 
# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.5)
example = 'Enzyme05'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.5,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.6)
example = 'Enzyme06'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.6,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.7)
example = 'Enzyme07'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.7,
        [example])
