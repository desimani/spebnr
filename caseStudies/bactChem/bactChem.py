# (bactChem.py) last update: 2023.07.02


from createSystem import *



example = 'BactChem'


# reactions: a list of object {reaction name + reaction rates (k)}
reactions = [
    {'a1': 1.15},
    {'a2': 0.25},
    {'a3': 0.1},
    {'a4': 10},
    {'a5': 0.002},
    {'a6': 1},
    {'a7': 1},
    {'a8': 80},
    {'a9': 0.01},
    {'a10': 0.2},
    {'a11': 1},
    {'a12': 0.25},
    {'a13': 1.15},
    {'a14': 0.18}
    ]



#LEVELS: no_min_level...no_max_level
no_min_level = 0
no_max_level = 50

SPECIES = [ # [species_name, min_level, max_level, roles for each reaction]
        # roles [no. molecules as R, no. molecules as P]
        # R --> [1,0]
        # P --> [0,1]
        # neither R nor P --> [0,0]
    
    ['L',no_min_level,int(100000*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],
                                    [1,1],[0,0],[0,0],[0,0],[0,0],
                                    [0,0],[0,0],[0,0],[0,0])],
    
    ['X',no_min_level,no_max_level,([1,0],[0,1],[0,1],[0,0],[0,0],
                                    [0,0],[0,0],[0,0],[0,0],[0,0],
                                    [0,0],[0,0],[0,0],[0,0])],
    
    ['X*',no_min_level,no_max_level,([0,1],[1,0],[1,0],[0,0],[0,0],
                                     [1,0],[0,0],[0,0],[0,0],[0,0],
                                     [0,1],[0,0],[0,0],[0,0])],
    
    ['CheY',no_min_level,no_max_level,([0,0],[0,0],[1,0],[0,1],[0,1],
                                       [0,0],[0,0],[0,0],[0,0],[0,0],
                                       [0,0],[0,0],[0,0],[1,0])],
    
    ['CheZ',no_min_level,no_max_level,([0,0],[0,0],[0,0],[1,1],[0,0],
                                       [0,0],[0,0],[0,0],[0,0],[0,0],
                                       [0,0],[0,0],[0,0],[0,0])],
    
    ['XL',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                     [0,0],[0,0],[0,1],[1,0],[0,0],
                                     [0,0],[0,0],[0,0],[0,0])],
    
    ['X*m',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                      [0,0],[0,0],[0,0],[0,1],[1,1],
                                      [1,0],[1,0],[0,1],[1,0])],
    
    ['CheR',no_min_level,int(1000*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],
                                       [0,0],[0,0],[0,0],[1,1],[0,0],
                                       [0,0],[0,0],[0,0],[0,0])],
    
    ['XY',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                     [0,1],[0,0],[1,0],[0,0],[0,0],
                                     [0,0],[0,0],[0,0],[0,0])],
    
    ['CheB',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                       [0,0],[0,1],[0,0],[0,0],[1,0],
                                       [0,0],[0,0],[0,0],[0,0])],
    
    ['CheBp',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                        [0,0],[1,0],[0,0],[0,0],[0,1],
                                        [1,1],[0,0],[0,0],[0,0])],
    
    ['Xm',no_min_level,no_max_level,([0,0],[0,0],[0,0],[0,0],[0,0],
                                     [0,0],[0,0],[0,0],[0,0],[0,0],
                                     [0,0],[0,1],[1,0],[0,1])],
    
    ['CheYp',no_min_level,int(10*1.1),([0,0],[0,0],[0,1],[1,0],[1,0],
                                        [0,0],[0,0],[0,0],[0,0],[0,0],
                                        [0,0],[0,0],[0,0],[0,1])]
    
    ]



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the input distance
"""
# input: L
max_level_L = SPECIES[0][2]
def rank_input_L(p):
    # input: L
    index_L = 0 # always
    # p is always ordered
    
    max_level_input = max_level_L
    l_split = p[index_L].split('_')
    level_input = int(l_split[len(l_split)-1])
    return [level_input / max_level_input] 



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the input distance
"""
# input: CheR
max_level_CheR = SPECIES[7][2]
def rank_input_CheR(p):
    
    index_CheR = 7 # always
    
    # p is always ordered
    
    max_level_input = max_level_CheR
    _split = p[index_CheR].split('_')
    level_input = int(_split[len(_split)-1])
    return [level_input / max_level_input] 



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the output distance
"""
# output: CheYp
index_output = 12 # always
max_level_CheYp = SPECIES[12][2]
def rank_output(p):
    
    # p is always ordered
    
    max_level_output = max_level_CheYp
    CheYp_split = p[index_output].split('_')
    level_output = int(CheYp_split[len(CheYp_split)-1])
    return [level_output / max_level_output] 



"""
p: the set of parallel processes in a specific time
return: the input concentration
"""
index_input = 0 ### L
def get_input_level(p):
    
    input_split = p[index_input].split('_')
    level_input = int(input_split[len(input_split)-1])
    
    return [level_input] # 1 size because there is only one output species



"""
p: the set of parallel processes in a specific time
return: the output concentration
"""
def get_output_level(p):
    
    CheYp_split = p[index_output].split('_')
    level_output = int(CheYp_split[len(CheYp_split)-1])
    
    return [level_output] # 1 size because there is only one output species




"""
return: a list of the output species name and the maximum output concentration
"""
def get_output_max_level():
    
    return [max_level_CheYp]



# input: L (or CheR)
# output: CheYp
INPUT_SPECIES = ['L'] # or INPUT_SPECIES = ['CheR']
OUTPUT_SPECIES = ['CheYp']

# initial set of parallel processes
# there is not the list of reactions
# spear.py reads the list of reaction via the map

p_original = ['L_1000','X_10','X*_10','CheY_10','CheZ_1',
              'XL_0','X*m_1','CheR_1000','XY_0','CheB_2',
              'CheBp_0','Xm_0','CheYp_1']

p_perturbated1 = ['L_100000','X_10','X*_10','CheY_10','CheZ_1',
              'XL_0','X*m_1','CheR_1000','XY_0','CheB_2',
              'CheBp_0','Xm_0','CheYp_1'] # L100,000, R1000

p_perturbated2 = ['L_1000','X_10','X*_10','CheY_10','CheZ_1',
              'XL_0','X*m_1','CheR_100','XY_0','CheB_2',
              'CheBp_0','Xm_0','CheYp_1'] # L1000, R100





### FIND MAX LEVELS
"""
h=30000

# original system
rank_input = rank_input_L
N=1
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
N=30
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])


# 1st perturbed system
N=1
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
N=30
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated1, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])


# 2nd perturbed system
rank_input = rank_input_CheR
INPUT_SPECIES = ['CheR']
N=1
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
N=30
# find min and max level of each species
find_max_level( SPECIES, reactions, p_perturbated2, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
"""





### COMPARE ORIGINAL - PERTURBED
"""
h=15000
N=50
l=2


# perturbated 1
example = 'BactChemCompare1st'
example0 = 'BactChemL1000CheR1000'
example1 = 'BactChemL100000CheR1000'
rank_input = rank_input_L
INPUT_SPECIES = ['L']
compare_two_systems( SPECIES, reactions, p_original, p_perturbated1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example1])

# perturbated 2
example = 'BactChemCompare2nd'
example2 = 'BactChemL1000CheR100'
rank_input = rank_input_CheR
INPUT_SPECIES = ['CheR']
compare_two_systems( SPECIES, reactions, p_original, p_perturbated2,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example2])

"""





### ESTIMATING ROBUSTNESS

h=15000
N=50
l=2
no_perturbated_systems = 20
p1 = p_original

# input L
rank_input = rank_input_L
INPUT_SPECIES = ['L']

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.3)
example = 'BactChemL03'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.3,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.4)
example = 'BactChemL04'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.4,
        [example])
 
# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.5)
example = 'BactChemL05'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.5,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.6)
example = 'BactChemL06'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.6,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.7)
example = 'BactChemL07'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.7,
        [example])




# input CheR
rank_input = rank_input_CheR
INPUT_SPECIES = ['CheR']
# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.3)
example = 'BactChemCheR03'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.3,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.4)
example = 'BactChemCheR04'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.4,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.5)
example = 'BactChemCheR05'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.5,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.6)
example = 'BactChemCheR06'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.6,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.7)
example = 'BactChemCheR07'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.7,
        [example])
