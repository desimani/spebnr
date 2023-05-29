# (EnvZOmpR.py) last update: 2023.05.18


from createSystem import *



example = 'EnvZOmpRsystem'


# reactions: a list of object {reaction name + reaction rates (k)}
reactions = [
    {'a1': 0.5},
    {'a2': 0.5},
    {'a3': 0.5},
    {'a4': 0.5},
    {'a5': 0.1},
    {'a6': 0.02},
    {'a7': 0.5},
    {'a8': 0.5},
    {'a9': 0.02},
    {'a10': 0.5},
    {'a11': 0.1}
    ]


# LEVELS: no_min_level...no_max_level
#e.g. no_max_level = int(2000*1.1)
no_min_level = 0
SPECIES = [ # [species_name, min_level, max_level, roles for each reaction]
        # roles [no. molecules as R, no. molecules as P]
        # R --> [1,0]
        # P --> [0,1]
        # N --> [0,0]
    ['X',no_min_level,int(252*1.1),([0,1],[1,0],[0,1],[1,0],[0,0],[0,0],[0,0],[0,1],[0,0],[0,0],[0,0])],
    ['Y',no_min_level,int(1003*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],[1,0],[0,1],[0,0],[0,0],[0,0],[0,1])],
    ['XT',no_min_level,int(114*1.1),([0,0],[0,0],[1,0],[0,1],[1,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0])],
    ['XP',no_min_level,int(21*1.1),([0,0],[0,0],[0,0],[0,0],[0,1],[1,0],[0,1],[0,0],[0,0],[0,0],[0,0])],
    ['XPY',no_min_level,int(33*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],[0,1],[1,0],[1,0],[0,0],[0,0],[0,0])],
    ['YP',no_min_level,int(53*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,1],[1,0],[0,1],[0,0])],
    ['XDYP',no_min_level,int(98*1.1),([0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,1],[1,0],[1,0])],
    ['XD',no_min_level,int(133*1.1),([1,0],[0,1],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[1,0],[0,1],[0,1])]
    ]



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the input distance
"""
# input: X, Y
index_X = 0 # always
index_Y = 1 # always
max_level_X = SPECIES[index_X][2]
min_level_X = SPECIES[index_X][1]
max_level_Y = SPECIES[index_Y][2]
min_level_Y = SPECIES[index_Y][1]
def rank_input( p ):
    
    # p is always ordered
    
    #print('p[index_X]: ',p[index_X],' p[index_Y]: ',p[index_Y])
    
    x_split = p[index_X].split('_')
    y_split = p[index_Y].split('_')
    
    level_input_X = int(x_split[len(x_split)-1])
    level_input_Y = int(y_split[len(y_split)-1])

    partX = (level_input_X / max_level_X) # if min_level_X == 0
    partY = (level_input_Y / max_level_Y) # if min_level_Y == 0

    #partX = ( level_input_X / (max_level_X-min_level_X) ) # if min_level_X != 0
    #partY = ( level_input_Y / (max_level_Y-min_level_Y) ) # if min_level_Y != 0

    return [partX, partY] 



"""
p: the set of parallel processes in a specific time
return: a value used to calculate the output distance
"""
# output: YP
index_output = 5 # always
max_level_YP = SPECIES[index_output][2]
min_level_YP = SPECIES[index_output][1]
def rank_output( p ):
    
    # p is always ordered
    
    #print('p[index_output]: ',p[index_output])
    
    yp_split = p[index_output].split('_')
    level_output = int(yp_split[len(yp_split)-1])
    
    return [level_output / (max_level_YP-min_level_YP) ] # if min_level_YP == 0



"""
p: the set of parallel processes in a specific time
return: a list of the maximum output concentration
"""
def get_input_level(p):

    # X
    x_split = p[index_X].split('_')
    level_x = int(x_split[len(x_split)-1])

    # Y
    y_split = p[index_Y].split('_')
    level_y = int(y_split[len(y_split)-1])
    
    return [level_x, level_y] # there are two input species



"""
p: the set of parallel processes in a specific time
return: a list of the maximum output concentration
"""
def get_output_level(p):
    
    yp_split = p[index_output].split('_')
    level_output = int(yp_split[len(yp_split)-1])
    
    return [level_output] # 1 size because there is only one output species



"""
return: a list of the output species name and the maximum output concentration
"""
def get_output_max_level():
    
    return [max_level_YP]




# parameters
h = 15000
N = 50
l = 2
no_perturbated_systems = 20 # 1 <= no_perturbated_systems <= 100


# CHECK ROBUSTNESS

# input: X, Y
# output: YP
INPUT_SPECIES = ['X', 'Y']
OUTPUT_SPECIES = ['YP']

# initial set of parallel processes
# there is not the list of reactions
# createSystem.py reads the list of reaction via the map

# original
p1 = ['X_25','Y_150','XT_0','XP_0','XPY_0','YP_10','XDYP_0','XD_50']

p_original = p1

"""
# find min and max level of each species
find_max_level( SPECIES, reactions, p_original, 1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N,
                        get_input_level, get_output_level,
                        rank_input, rank_output,
                        [example])
"""
# perturbated 1
example = 'EnvZOmpRcompare1st'
example0 = 'EnvZOmpRX25Y150'
example1 = 'EnvZOmpRsystemX10Y50'
p_perturbated1 = ['X_10','Y_50','XT_0','XP_0','XPY_0','YP_10','XDYP_0','XD_50']
compare_two_systems( SPECIES, reactions, p_original, p_perturbated1,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example1])

# perturbated 2
example = 'EnvZOmpRcompare2nd'
example2 = 'EnvZOmpRX250Y1000'
p_perturbated2 = ['X_250','Y_1000','XT_0','XP_0','XPY_0','YP_10','XDYP_0','XD_50']
compare_two_systems( SPECIES, reactions, p_original, p_perturbated2,
                        INPUT_SPECIES, OUTPUT_SPECIES, h, N, l,
                        rank_input, rank_output,
                        get_input_level, get_output_level, get_output_max_level,
                        [example, example0, example1])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.3)
example = 'EnvZOmpR03'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.3,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.4)
example = 'EnvZOmpR04'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.4,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.5)
example = 'EnvZOmpR05'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.5,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.6)
example = 'EnvZOmpR06'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.6,
        [example])

# compare the original system with |perturbated systems|=no_perturbated_systems (with eta_1=0.7)
example = 'EnvZOmpR07'
execute_robustness(SPECIES, reactions,
        p1, no_perturbated_systems,
        INPUT_SPECIES, OUTPUT_SPECIES,
        h, N, l,
        rank_input, rank_output,
        get_output_level, get_output_max_level,
        0.7,
        [example])


