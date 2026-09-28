import math
import numpy as np
from tqdm import tqdm
from itertools import permutations
import copy
from .basic_tools_jit import R3, R3_table_matrix, diagram_generate_param, oppend
from .realisable_Gauss_signed import real_Gauss_signed


def index_outside(G, ind):
    if G[ind + 1, 0] == 0 and G[ind - 1, 0] != 0 and G[ind, 0] != 0:
        return ind + 1
    elif G[ind - 1, 0] == 0 and G[ind + 1, 0] != 0 and G[ind, 0] != 0:
        return ind - 1
    else:
        print("error")
        return


def initial_position(A):
    """
    param: A : R3_table_matrix[:,:,*]
    The function returns the "6" initial position matrix of a meridian along the tetrahedion equation
    """
    result = np.zeros([12, 3, 6], int)
    count = 0
    for i in range(3):  # i represent the side we choose
        #         print('i=',i)
        for j in range(2):
            #             print('j=',j)
            # 2*i+j is the crossing we choose
            if oppend(A, 2 * i + j) % 2 == 0:
                #                 print('ok')
                B = np.insert(A, [0, 2, 2, 4, 4, 6], np.zeros([6, 3]), axis=0)
                #                 print('B=')
                #                 print(B)
                ind_old_to_new = np.array([1, 2, 5, 6, 9, 10], int)
                signs = np.ones(3)
                index = np.zeros(3, int)
                index[0] = index_outside(B, ind_old_to_new[2 * i + j])
                index[1] = index_outside(B, ind_old_to_new[oppend(A, 2 * i + j)])
                index[2] = index_outside(B, ind_old_to_new[oppend(A, 2 * i + 1 - j)])
                #                 print('index=',index)
                # signs[0]=signs[1] always since we require that oppend(A,2*i+j) %2 == 0
                # i.e., the purple TAIL rule
                signs[1] = (1 - 2 * j) * signs[0]
                signs[2] = (
                    (1 - 2 * j)
                    * A[2 * i + j, 1]
                    * A[2 * i + j, 2]
                    * A[2 * i + 1 - j, 1]
                    * A[2 * i + 1 - j, 2]
                    * signs[0]
                )

                for s in [1, -1]:  # s represent the signe of the A point
                    signs_final = signs * s
                    piece = np.zeros([3, 3], int)
                    for k in range(3):
                        B[index[k], 0] = k + 4
                        B[index[k], 1] = 1
                        B[index[k], 2] = signs_final[k]
                        if (
                            A[2 * i + j, 1] * A[2 * i + j, 2] * s == 1
                        ):  # when s = 1, going left makes reversed order
                            piece[2 - k, 0] = k + 4
                            piece[2 - k, 1] = -1
                            piece[2 - k, 2] = signs_final[k]
                        else:
                            piece[k, 0] = k + 4
                            piece[k, 1] = -1
                            piece[k, 2] = signs_final[k]
                    #                     print('B=')
                    #                     print(B)
                    temp = np.array([], int)
                    for t in range(B.shape[0]):
                        if B[t, 0] == 0:
                            temp = np.append(temp, t)
                    C = np.delete(B, temp, axis=0)
                    result[:, :, count] = np.append(C, piece, axis=0)
                    count += 1
    #                 if 2*i+j %2 == 0:#the choosing point is the initial point on the side chosen
    #                     if oppend(A,2*i+1-j)%2 == 0:
    return result


def R3_index_sort(G, arr):
    """
    param: G: Gauss matrix
    param: arr: an array of index in G indicating the possible crossings to perform R3
    Output: the right index to perform R3 for the function R3()
    """
    arr_opp = np.zeros(len(arr), int)
    for i in range(len(arr)):
        arr_opp[i] = oppend(G, arr[i])
    temp = np.append(arr, arr_opp)
    temp = np.unique(temp)
    temp_sorted = np.sort(temp)
    #     print('temp_sorted=',temp_sorted)
    if len(temp_sorted) != 6:
        print("the arr cannot do R3")
    else:
        return temp_sorted[[0, 2, 4]]


def slide_direct(ind):
    if ind % 3 == 0:
        return 1
    if ind % 3 == 2:
        return -1


def third_index_find(G, ind_3, directions):
    ind_next = ind_3 + directions
    if (
        G[ind_next[0], 0] == G[ind_next[1], 0]
        and int(ind_next[0] / 3) == int(ind_3[0] / 3)
        and int(ind_next[1] / 3) == int(ind_3[1] / 3)
    ):
        return [ind_3[0], ind_3[1], ind_next[0]]
    elif (
        G[ind_next[0], 0] == G[ind_next[2], 0]
        and int(ind_next[0] / 3) == int(ind_3[0] / 3)
        and int(ind_next[2] / 3) == int(ind_3[2] / 3)
    ):
        return [ind_3[0], ind_3[2], ind_next[0]]
    elif (
        G[ind_next[1], 0] == G[ind_next[2], 0]
        and int(ind_next[1] / 3) == int(ind_3[1] / 3)
        and int(ind_next[2] / 3) == int(ind_3[2] / 3)
    ):
        return [ind_3[1], ind_3[2], ind_next[2]]
    else:
        print("error with finding third ind for R3")
        return


def R3_loop(initial):
    cases = np.zeros([12, 3, 8])
    R3_index_cases = np.zeros([8, 3])
    ind_3 = np.array([oppend(initial, 9), oppend(initial, 10), oppend(initial, 11)], int)
    directions = np.array([slide_direct(ind_3[0]), slide_direct(ind_3[1]), slide_direct(ind_3[2])])
    diag = initial.copy()
    for i in range(3):
        #         print('ind_3=',ind_3)
        #         print('directions=',directions)
        temp_ind_R3 = third_index_find(diag, ind_3, directions)
        ind_R3_sorted = R3_index_sort(diag, temp_ind_R3)
        #         print('temp_ind_R3=',temp_ind_R3)
        #         draw(diag,ind_R3_sorted,ind_R3_sorted,ind_3)
        #         print(diag)

        cases[:, :, i] = copy.deepcopy(diag)
        R3_index_cases[i, :] = copy.deepcopy(ind_R3_sorted)

        # [diag,_,index_old_to_new,_] = R3(diag,ind_R3_sorted[0],ind_R3_sorted[1],ind_R3_sorted[2])
        [diag, _, index_old_to_new, _] = R3(diag, ind_R3_sorted)
        ind_3 = index_old_to_new[ind_3]
    # middle move
    oppend_arr = np.array([oppend(diag, ind_3[0]), oppend(diag, ind_3[1]), oppend(diag, ind_3[2])])
    temp_arr = np.append(ind_3, oppend_arr, axis=0)
    ind_R3_middle_temp = np.array([], int)
    for s in range(12):
        if (s in temp_arr) == 0:
            ind_R3_middle_temp = np.append(ind_R3_middle_temp, [s])
    ind_R3_middle = R3_index_sort(diag, ind_R3_middle_temp)

    cases[:, :, 3] = copy.deepcopy(diag)
    R3_index_cases[3, :] = copy.deepcopy(ind_R3_middle)

    # [diag,_,index_old_to_new,_] = R3(diag,ind_R3_middle[0],ind_R3_middle[1],ind_R3_middle[2])
    [diag, _, index_old_to_new, _] = R3(diag, ind_R3_middle)
    ind_3 = index_old_to_new[ind_3]
    #     draw(diag,'',ind_3)
    # another half
    directions = -directions
    for i in range(3):
        temp_ind_R3 = third_index_find(diag, ind_3, directions)
        ind_R3_sorted = R3_index_sort(diag, temp_ind_R3)
        #         print('temp_ind_R3=',temp_ind_R3)
        #         draw(diag,ind_R3_sorted,ind_R3_sorted,ind_3)
        #         print(diag)
        cases[:, :, 4 + i] = copy.deepcopy(diag)
        R3_index_cases[4 + i, :] = copy.deepcopy(ind_R3_sorted)

        # [diag,_,index_old_to_new,_] = R3(diag,ind_R3_sorted[0],ind_R3_sorted[1],ind_R3_sorted[2])
        [diag, _, index_old_to_new, _] = R3(diag, ind_R3_sorted)
        ind_3 = index_old_to_new[ind_3]
    # last move
    oppend_arr = np.array([oppend(diag, ind_3[0]), oppend(diag, ind_3[1]), oppend(diag, ind_3[2])])
    temp_arr = np.append(ind_3, oppend_arr, axis=0)
    ind_R3_middle_temp = np.array([], int)
    for s in range(12):
        if (s in temp_arr) == 0:
            ind_R3_middle_temp = np.append(ind_R3_middle_temp, [s])
    ind_R3_middle = R3_index_sort(diag, ind_R3_middle_temp)

    cases[:, :, 7] = copy.deepcopy(diag)
    R3_index_cases[7, :] = copy.deepcopy(ind_R3_middle)

    # [diag,_,index_old_to_new,_] = R3(diag,ind_R3_middle[0],ind_R3_middle[1],ind_R3_middle[2])
    [diag, _, index_old_to_new, _] = R3(diag, ind_R3_middle)
    ind_3 = index_old_to_new[ind_3]
    #     draw(diag,'',ind_3)
    return [cases, R3_index_cases]


tetrahedron_intitial_positions = np.zeros([12, 3, 48], int)
for i in range(8):
    res_temp = initial_position(R3_table_matrix[:, :, 2 * i])
    for j in range(6):
        tetrahedron_intitial_positions[:, :, 6 * i + j] = copy.deepcopy(res_temp[:, :, j])
# for i in range(48):
#     draw(tetrahedron_intitial_positions[:,:,i],i)


local_tetrahedron_cases = np.zeros([12, 3, 8, 48], int)
local_tetrahedron_R3_index_cases = np.zeros([8, 3, 48], int)
for i in range(48):
    [local_tetrahedron_cases[:, :, :, i], local_tetrahedron_R3_index_cases[:, :, i]] = R3_loop(
        tetrahedron_intitial_positions[:, :, i]
    )

global_tetrahedron_cases = np.zeros([12, 3, 8, 48 * 24], int)
global_tetrahedron_R3_index_cases = np.zeros([8, 3, 48 * 24], int)

count_p = 0
for p in permutations([0, 1, 2, 3]):
    arr = np.zeros(12, int)
    for l in range(12):
        arr[l] = p[l // 3] * 3 + (l % 3)
    for s in range(48):
        for t in range(8):
            for i in range(12):
                global_tetrahedron_cases[arr[i], :, t, 48 * count_p + s] = copy.deepcopy(
                    local_tetrahedron_cases[i, :, t, s]
                )
                global_tetrahedron_R3_index_cases[t, :, 48 * count_p + s] = np.sort(
                    arr[local_tetrahedron_R3_index_cases[t, :, s]]
                )
    count_p += 1


def generate_pre_tetrahedron_loops(global_loop, global_index, numb_arrows):
    """
    param: global_loop[12,3,8]
    param: global_index[8,3]
    ouput:
    """
    [chord, directed, directed_signed, numb_chord, count_directed, count_directed_signed] = (
        diagram_generate_param(numb_arrows)
    )
    # change names
    for s in range(count_directed_signed):
        for r in range(2 * numb_arrows):
            directed_signed[r, 0, s] += 6
    # name change finished
    # the diagram first generate all the n-arrow diagram, then put inside the 4 pieces for tetrahedron.
    # equivalent to choose 4 different elements from 2n+4 elements.
    cases = np.zeros(
        [2 * (numb_arrows + 6), 3, 8, count_directed_signed * math.comb(2 * numb_arrows + 4, 4)],
        int,
    )
    R3_index = np.zeros([8, 3, cases.shape[3]], int)
    s = 0
    for i in range(2 * numb_arrows + 1):
        for j in range(i, 2 * numb_arrows + 1):
            for k in range(j, 2 * numb_arrows + 1):
                for l in range(k, 2 * numb_arrows + 1):
                    for w in range(count_directed_signed):
                        for t in range(8):
                            cases[:, :, t, s] = np.insert(
                                directed_signed[:, :, w],
                                [i, i, i, j, j, j, k, k, k, l, l, l],
                                global_loop[:, :, t],
                                axis=0,
                            )
                            index_old_to_new = np.zeros(12, int)
                            for ind in range(12):
                                if ind < 3:
                                    index_old_to_new[ind] = i + ind
                                if ind >= 3 and ind < 6:
                                    index_old_to_new[ind] = j + 3 + ind - 3
                                if ind >= 6 and ind < 9:
                                    index_old_to_new[ind] = k + 6 + ind - 6
                                if ind >= 9:
                                    index_old_to_new[ind] = l + 9 + ind - 9
                            #                             print('index_old_to_new=',index_old_to_new)
                            #                             print('global_index[t,:]=',global_index[t,:])
                            R3_index[t, :, s] = index_old_to_new[global_index[t, :]]
                        s += 1
    return [cases, R3_index]


def generate_real_tetrahedron_loops(global_loop, global_index, numb_arrows):
    """
    param: global_loop[12,3,8]
    param: global_index[8,3]

    This function generate all realisible tetrahedron loops adding n-number of arrows.
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

    [pre_cases, pre_R3_index] = generate_pre_tetrahedron_loops(
        global_loop, global_index, numb_arrows
    )
    cases = np.zeros(pre_cases.shape, int)
    R3_index = np.zeros(pre_R3_index.shape, int)
    counter = 0
    for s in range(pre_cases.shape[3]):
        if real_Gauss_signed(pre_cases[:, :, 0, s]):
            cases[:, :, :, counter] = pre_cases[:, :, :, s]
            R3_index[:, :, counter] = pre_R3_index[:, :, s]
            counter += 1
    return [cases[:, :, :, 0:counter], R3_index[:, :, 0:counter]]


def real_tetrahedron_loops_collection(numb_arrows):
    case_blocks, index_blocks = [], []
    for s in tqdm(range(global_tetrahedron_cases.shape[3])):
        [cases_temp, index_temp] = generate_real_tetrahedron_loops(
            global_tetrahedron_cases[:, :, :, s],
            global_tetrahedron_R3_index_cases[:, :, s],
            numb_arrows,
        )
        case_blocks.append(cases_temp)
        index_blocks.append(index_temp)
    return np.concatenate(case_blocks, axis=3), np.concatenate(index_blocks, axis=2)


def virtual_tetrahedron_loops_collection(numb_arrows):
    case_blocks, index_blocks = [], []
    for s in tqdm(range(global_tetrahedron_cases.shape[3])):
        [cases_temp, index_temp] = generate_pre_tetrahedron_loops(
            global_tetrahedron_cases[:, :, :, s],
            global_tetrahedron_R3_index_cases[:, :, s],
            numb_arrows,
        )
        case_blocks.append(cases_temp)
        index_blocks.append(index_temp)
    return np.concatenate(case_blocks, axis=3), np.concatenate(index_blocks, axis=2)
    return (cases_temp, index_temp)


