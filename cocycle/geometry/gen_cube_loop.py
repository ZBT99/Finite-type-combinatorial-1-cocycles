import numpy as np
from .realisable_Gauss_signed import real_Gauss_signed
from .basic_tools_jit import R3, diagram_generate_param, infty_avancer
from .gen_tetra_loop import R3_index_sort
import math
from tqdm import tqdm

cube_table_matrix = np.zeros([8, 3, 48], int)

temp_matrix = np.array(
    [
        [1, -1, -1],
        [2, -1, -1],
        [3, -1, -1],
        [2, 1, -1],
        [4, -1, 1],
        [3, 1, -1],
        [1, 1, -1],
        [4, 1, 1],
    ],
    int,
)
cube_table_matrix[:, :, 0] = temp_matrix.copy()

temp_matrix = np.array(
    [
        [1, -1, 1],
        [2, -1, -1],
        [3, -1, 1],
        [2, 1, -1],
        [4, -1, -1],
        [4, 1, -1],
        [1, 1, 1],
        [3, 1, 1],
    ],
    int,
)
cube_table_matrix[:, :, 1] = temp_matrix.copy()

temp_matrix = np.array(
    [
        [1, -1, -1],
        [2, -1, 1],
        [3, -1, -1],
        [2, 1, 1],
        [4, -1, 1],
        [4, 1, 1],
        [1, 1, -1],
        [3, 1, -1],
    ],
    int,
)
cube_table_matrix[:, :, 2] = temp_matrix.copy()

temp_matrix = np.array(
    [[1, -1, 1], [2, -1, 1], [3, -1, 1], [2, 1, 1], [4, -1, -1], [3, 1, 1], [1, 1, 1], [4, 1, -1]],
    int,
)
cube_table_matrix[:, :, 3] = temp_matrix.copy()


# examiner(temp_matrix)
# draw(infty_avancer(temp_matrix,-3))
# for p in [2,4]:
#     temp_R3_ind = R3_index_sort(temp_matrix,np.array([0,1,p],int))
#     R3(temp_matrix,temp_R3_ind[0],temp_R3_ind[1],temp_R3_ind[2])


####4 diagrams in line 2
for s in range(4, 8):
    cube_table_matrix[:, :, s] = cube_table_matrix[:, :, s - 4].copy()
    for i in range(8):
        if cube_table_matrix[i, 0, s] in [3, 4]:
            cube_table_matrix[i, 1, s] = -cube_table_matrix[i, 1, s]
            cube_table_matrix[i, 2, s] = -cube_table_matrix[i, 2, s]
    # examiner(cube_table_matrix[:,:,s])

# s = 7
# draw(infty_avancer(cube_table_matrix[:,:,s],-3),s)

####

####4 diagrams in line 3
for s in range(8, 12):
    cube_table_matrix[:, :, s] = cube_table_matrix[:, :, s - 4].copy()
    for i in range(8):
        if cube_table_matrix[i, 0, s] == 1:
            cube_table_matrix[i, 1, s] = -cube_table_matrix[i, 1, s]
            cube_table_matrix[i, 2, s] = -cube_table_matrix[i, 2, s]
    # examiner(cube_table_matrix[:,:,s])
# s = 11
# draw(infty_avancer(cube_table_matrix[:,:,s],-3),s)
# print(cube_table_matrix[:,:,s])

#####

for s in range(12, 24):
    cube_table_matrix[:, :, s] = -cube_table_matrix[:, :, s - 12].copy()
    for i in range(8):
        cube_table_matrix[i, 0, s] = -cube_table_matrix[i, 0, s]
    # examiner(cube_table_matrix[:,:,s])

# s = 23
# draw(infty_avancer(cube_table_matrix[:,:,s],-3),s)
# print(cube_table_matrix[:,:,s])
####

##############
##############

for s in range(24, 48):
    for i in range(0, 2):
        cube_table_matrix[i, :, s] = cube_table_matrix[i, :, s - 24].copy()
    for i in range(2, 5):
        cube_table_matrix[i, :, s] = cube_table_matrix[i + 3, :, s - 24].copy()
    for i in range(5, 8):
        cube_table_matrix[i, :, s] = cube_table_matrix[i - 3, :, s - 24].copy()
    # examiner(cube_table_matrix[:,:,s])

cube_table_R3_index = np.zeros(
    [2, 3, cube_table_matrix.shape[2]], int
)  # [phases,3_index,case_number]
for s in range(cube_table_matrix.shape[2]):
    count = 0
    for p in [2, 4]:
        #         print(R3_index_sort(cube_table_matrix[:,:,s],np.array([0,1,p],int)))
        cube_table_R3_index[count, :, s] = R3_index_sort(
            cube_table_matrix[:, :, s], np.array([0, 1, p], int)
        )
        count += 1
# s = 5
for s in range(48):
    #     print('s=',s)
    #     print(cube_table_R3_index[:,:,s])
    temp0 = cube_table_R3_index[0, :, s]
    R3(cube_table_matrix[:, :, s], temp0)
    temp1 = cube_table_R3_index[1, :, s]
    R3(cube_table_matrix[:, :, s], temp1)

global_cube_cases = np.zeros([8, 3, 3 * cube_table_matrix.shape[2]], int)
global_cube_R3_index = np.zeros([2, 3, 3 * cube_table_R3_index.shape[2]], int)
s = 0
for w in range(cube_table_matrix.shape[2]):
    for r in range(3):
        if r == 0:
            global_cube_cases[:, :, s] = infty_avancer(cube_table_matrix[:, :, w], 0)
            global_cube_R3_index[:, :, s] = np.sort((cube_table_R3_index[:, :, w] - 0) % 8)
        if r == 1:
            global_cube_cases[:, :, s] = infty_avancer(cube_table_matrix[:, :, w], 2)
            global_cube_R3_index[:, :, s] = np.sort((cube_table_R3_index[:, :, w] - 2) % 8)
        if r == 2:
            global_cube_cases[:, :, s] = infty_avancer(cube_table_matrix[:, :, w], 5)
            global_cube_R3_index[:, :, s] = np.sort((cube_table_R3_index[:, :, w] - 5) % 8)
        s += 1


def generate_pre_cube_loops(cube_matrix, cube_R3_index, numb_arrows):
    """
    param: cube_matrix[8,3] (a Gauss matrix)
    param: cube_R3_index[2,3] [phases,3col]
    param: numb_arrows
    output: [cases,R3_index]
    """
    [chord, directed, directed_signed, numb_chord, count_directed, count_directed_signed] = (
        diagram_generate_param(numb_arrows)
    )
    # change names
    for s in range(count_directed_signed):
        for r in range(2 * numb_arrows):
            directed_signed[r, 0, s] += 4
    # name change finished
    # the diagram first generate all the n-arrow diagram, then put inside the 3 pieces for cube loops.
    # equivalent to choose 3 different elements from 2n+3 elements.
    cases = np.zeros(
        [
            2 * (numb_arrows + 4),
            3,
            2,
            3 * count_directed_signed * math.comb(2 * numb_arrows + 3, 3),
        ],
        int,
    )
    R3_index = np.zeros([2, 3, cases.shape[3]], int)
    s = 0
    for i in range(2 * numb_arrows + 1):
        for j in range(i, 2 * numb_arrows + 1):
            for k in range(j, 2 * numb_arrows + 1):
                for w in range(count_directed_signed):
                    for r in range(3):
                        if r == 0:
                            temp_matrix = infty_avancer(cube_matrix, 0)
                            temp_R3_index = np.sort((cube_R3_index - 0) % 8)
                            #                             R3(temp_matrix,temp_R3_index[0,0],temp_R3_index[0,1],temp_R3_index[0,2])
                            #                             R3(temp_matrix,temp_R3_index[1,0],temp_R3_index[1,1],temp_R3_index[1,2])
                            for t in range(2):
                                cases[:, :, t, s] = np.insert(
                                    directed_signed[:, :, w],
                                    [i, i, j, j, j, k, k, k],
                                    temp_matrix,
                                    axis=0,
                                )
                                index_old_to_new = np.zeros(8, int)
                                for ind in range(8):
                                    if ind < 2:
                                        index_old_to_new[ind] = i + ind
                                    if ind >= 2 and ind < 5:
                                        index_old_to_new[ind] = j + ind
                                    if ind >= 5 and ind < 8:
                                        index_old_to_new[ind] = k + ind
                                R3_index[t, :, s] = index_old_to_new[temp_R3_index[t, :]]

                        if r == 1:
                            temp_matrix = infty_avancer(cube_matrix, 2)
                            temp_R3_index = np.sort((cube_R3_index - 2) % 8)
                            #                             R3(temp_matrix,temp_R3_index[0,0],temp_R3_index[0,1],temp_R3_index[0,2])
                            #                             R3(temp_matrix,temp_R3_index[1,0],temp_R3_index[1,1],temp_R3_index[1,2])
                            for t in range(2):
                                cases[:, :, t, s] = np.insert(
                                    directed_signed[:, :, w],
                                    [i, i, i, j, j, j, k, k],
                                    temp_matrix,
                                    axis=0,
                                )
                                index_old_to_new = np.zeros(8, int)
                                for ind in range(8):
                                    if ind < 3:
                                        index_old_to_new[ind] = i + ind
                                    if ind >= 3 and ind < 3 + 3:
                                        index_old_to_new[ind] = j + ind
                                    if ind >= 3 + 3 and ind < 8:
                                        index_old_to_new[ind] = k + ind
                                R3_index[t, :, s] = index_old_to_new[temp_R3_index[t, :]]

                        if r == 2:
                            temp_matrix = infty_avancer(cube_matrix, 5)
                            temp_R3_index = np.sort((cube_R3_index - 5) % 8)
                            #                             R3(temp_matrix,temp_R3_index[0,0],temp_R3_index[0,1],temp_R3_index[0,2])
                            #                             R3(temp_matrix,temp_R3_index[1,0],temp_R3_index[1,1],temp_R3_index[1,2])
                            for t in range(2):
                                cases[:, :, t, s] = np.insert(
                                    directed_signed[:, :, w],
                                    [i, i, i, j, j, k, k, k],
                                    temp_matrix,
                                    axis=0,
                                )
                                index_old_to_new = np.zeros(8, int)
                                for ind in range(8):
                                    if ind < 3:
                                        index_old_to_new[ind] = i + ind
                                    if ind >= 3 and ind < 3 + 2:
                                        index_old_to_new[ind] = j + ind
                                    if ind >= 3 + 2 and ind < 8:
                                        index_old_to_new[ind] = k + ind
                                R3_index[t, :, s] = index_old_to_new[temp_R3_index[t, :]]
                        s += 1

    #                         for t in range(2):
    #                             cases[:,:,t,s] = np.insert(directed_signed[:,:,w],[i,i,j,j,j,k,k,k],
    #                                                            temp_matrix,
    #                                                            axis = 0)
    #                             index_old_to_new = np.zeros(8,int)
    #                             for ind in range(8):
    #                                 if ind<2:
    #                                     index_old_to_new[ind] = i+ind
    #                                 if ind>=2 and ind<5:
    #                                     index_old_to_new[ind] = j+ind
    #                                 if ind>=5 and ind<8:
    #                                     index_old_to_new[ind] = k+ind
    #                             R3_index[t,:,s] = index_old_to_new[temp_R3_index[t,:]]

    #                     s += 1
    return [cases, R3_index]


def generate_real_cube_loops(cube_matrix, cube_R3_index, numb_arrows):
    """
    param: cube_matrix[8,3]
    param: cube_R3_index[2,3] [phases,3col]

    This function generate all realisible cube loops adding n-number of arrows.
    """
    #     if numb_arrows == 0:
    #         if real_Gauss_signed(global_loop[:,:,0]) == True:
    #             case = np.zeros([12,3,8,1],int)
    #             R3_index = np.zeros([8,3,1],int)
    #             case[:,:,:,1] = copy.deepcopy(global_loop)
    #             R3_index[:,:,1] = copy.deepcopy(global_index)
    #             return [case,R3_index]
    #         else:
    #             return [np.zeros([12,3,8,0],int),np.zeros([8,3,0],int)]

    [pre_cases, pre_R3_index] = generate_pre_cube_loops(cube_matrix, cube_R3_index, numb_arrows)
    cases = np.zeros(pre_cases.shape, int)
    R3_index = np.zeros(pre_R3_index.shape, int)
    counter = 0
    for s in range(pre_cases.shape[3]):
        if real_Gauss_signed(pre_cases[:, :, 0, s]):
            cases[:, :, :, counter] = pre_cases[:, :, :, s]
            R3_index[:, :, counter] = pre_R3_index[:, :, s]
            counter += 1
    return [cases[:, :, :, 0:counter], R3_index[:, :, 0:counter]]


def real_cube_loops_collection(numb_arrows):
    case_blocks, index_blocks = [], []
    for s in tqdm(range(cube_table_matrix.shape[2])):
        [cases_temp, index_temp] = generate_real_cube_loops(
            cube_table_matrix[:, :, s], cube_table_R3_index[:, :, s], numb_arrows
        )
        case_blocks.append(cases_temp)
        index_blocks.append(index_temp)
    return np.concatenate(case_blocks, axis=3), np.concatenate(index_blocks, axis=2)


def virtual_cube_loops_collection(numb_arrows):
    case_blocks, index_blocks = [], []
    for s in tqdm(range(cube_table_matrix.shape[2])):
        (cases_temp, index_temp) = generate_pre_cube_loops(
            cube_table_matrix[:, :, s], cube_table_R3_index[:, :, s], numb_arrows
        )
        case_blocks.append(cases_temp)
        index_blocks.append(index_temp)
    return np.concatenate(case_blocks, axis=3), np.concatenate(index_blocks, axis=2)


