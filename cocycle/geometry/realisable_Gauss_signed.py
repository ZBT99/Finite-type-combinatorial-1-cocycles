#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 14:45:56 2024
Realisable signed Gauss diagram
@author: butianzhang
"""

import numpy as np
from numba import njit

# =============================================================================
# def index_names_trans(G):
#     """
#     param: G: Gauss matrix
#     output: [name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist]
#     ind_namelist_to_ind_matrix[ind_namelist]=[first ind in matrix, second ind in matrix]
#     """
#     name_list = getname(G)
#     ind_namelist_to_ind_matrix = np.zeros([name_list.shape[0],2],int)
#     ind_matrix_to_ind_namelist = np.zeros(G.shape[0],int)
#     for i in range(name_list.shape[0]):
#         occur_time = 0
#         for j in range(G.shape[0]):
#             if name_list[i] == G[j,0]:
#                 ind_namelist_to_ind_matrix[i,occur_time] = j
#                 occur_time += 1
#                 ind_matrix_to_ind_namelist[j] = i
#     return [name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist]
#
#
# =============================================================================


@njit(cache=True)
def index_names_trans(G):
    """
    param: G: Gauss matrix
    output: (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist)
    ind_namelist_to_ind_matrix[ind_namelist]=[first ind in matrix, second ind in matrix]
    """
    numb_cr = G.shape[0] // 2
    name_list = np.zeros(numb_cr, dtype=np.int64)
    ind_namelist_to_ind_matrix = np.zeros((numb_cr, 2), dtype=np.int64)
    ind_matrix_to_ind_namelist = np.zeros(G.shape[0], dtype=np.int64)
    count_name = 0
    for i in range(G.shape[0]):
        appeared = False
        for j in range(0, count_name):
            if name_list[j] == G[i, 0]:
                appeared = True
                ind_matrix_to_ind_namelist[i] = j
                ind_namelist_to_ind_matrix[j, 1] = i
                break
        if not appeared:
            name_list[count_name] = G[i, 0]
            ind_namelist_to_ind_matrix[count_name, 0] = i
            ind_matrix_to_ind_namelist[i] = count_name
            count_name += 1
    return (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist)


# =============================================================================
# # @njit(cache=True)
# def Gauss_to_band(G):
#     """
#     G: Gauss diagram
#     Output:
#     """
#     (G_name_list, G_ind_namelist_to_ind_matrix, G_ind_matrix_to_ind_namelist) = index_names_trans(G)
#     oppend_G = np.zeros(G.shape[0],dtype=np.int64)
#     len_name_list = len(G_name_list)
# # =============================================================================
# #     for i in range(len_name_list):
# #         oppend_G[G_ind_namelist_to_ind_matrix[i,0]] = G_ind_namelist_to_ind_matrix[i,1]
# #         oppend_G[G_ind_namelist_to_ind_matrix[i,1]] = G_ind_namelist_to_ind_matrix[i,0]
# # =============================================================================
#     oppend_G[G_ind_namelist_to_ind_matrix[:,0]] = G_ind_namelist_to_ind_matrix[:,1]
#     oppend_G[G_ind_namelist_to_ind_matrix[:,1]] = G_ind_namelist_to_ind_matrix[:,0]
#     crossings_check_list = np.zeros(len_name_list,dtype=np.int64)
#     arc_direction = np.ones(G.shape[0]+1,dtype=np.int64)
#     path = np.arange(G.shape[0]+1)# written for the arcs
#     name_crossings = np.zeros(G.shape[0],dtype=np.int64)#using the namelist index to call crossings
#     left_arcs_crossings = np.zeros((len_name_list,2),dtype=np.int64)#using the namelist index to call crossings
#     right_arcs_crossings = np.zeros((len_name_list,2),dtype=np.int64)#using the namelist index to call crossings
#
#
#
#     step = 0
#     step_counter = 0
#     while step != G.shape[0]:
#         if arc_direction[step] == 1:#anti_clockwise on the arc
#             name_crossings[step_counter] = G_ind_matrix_to_ind_namelist[step]
#             if crossings_check_list[G_ind_matrix_to_ind_namelist[step]] == 0: #first time come to a crossing
#                 crossings_check_list[G_ind_matrix_to_ind_namelist[step]] = 1
#                 if arc_direction[oppend_G[step]]==1: #otherside is anticlockwise
#                     left_arcs_crossings[G_ind_matrix_to_ind_namelist[step],:] = np.array([step,oppend_G[step]],dtype=np.int64)
#
#                     index_begin_change_in_path = np.where(path==oppend_G[step])[0][0]
#                     index_end_change_in_path = np.where(path==step+1)[0][0]
#                     temp_min = np.min([index_begin_change_in_path,index_end_change_in_path])
#                     temp_max = np.max([index_begin_change_in_path,index_end_change_in_path])
#                     temp_path = path.copy()
#                     for i in range(temp_min,temp_max+1):
#                         arc_direction[path[i]] = -arc_direction[path[i]]
#                         temp_path[i] = path[temp_min+temp_max-i]
#                     path = temp_path
#                     step = oppend_G[step].copy()
#
#                 else: #otherside is clockwise
#                     left_arcs_crossings[G_ind_matrix_to_ind_namelist[step],:] = np.array([step,1+oppend_G[step]],dtype=np.int64)
#
#                     index_begin_change_in_path = np.where(path==oppend_G[step]+1)[0][0]
#                     index_end_change_in_path = np.where(path==step+1)[0][0]
#                     temp_min = np.min([index_begin_change_in_path,index_end_change_in_path])
#                     temp_max = np.max([index_begin_change_in_path,index_end_change_in_path])
#                     temp_path = path.copy()
#                     for i in range(temp_min,temp_max+1):
#                         arc_direction[path[i]] = -arc_direction[path[i]]
#                         temp_path[i] = path[temp_min+temp_max-i]
#                     path = temp_path
#                     step = 1+oppend_G[step]
#
#             else: #second time to the crossing
#                 right_arcs_crossings[G_ind_matrix_to_ind_namelist[step],:] = np.array([step,path[1+step_counter]],dtype=np.int64)
#
#                 step = path[1+step_counter]
#
#
#         else: #clockwise on the arc
#             name_crossings[step_counter] = G_ind_matrix_to_ind_namelist[step-1]
#             if crossings_check_list[G_ind_matrix_to_ind_namelist[step-1]] == 0: #first time come to a crossing
#                 crossings_check_list[G_ind_matrix_to_ind_namelist[step-1]] = 1
#                 if arc_direction[oppend_G[step-1]]==-1:#otherside is clockwise
#                     left_arcs_crossings[G_ind_matrix_to_ind_namelist[step-1],:] = np.array([step,oppend_G[step-1]+1],dtype=np.int64)
#
#                     index_begin_change_in_path = np.where(path==oppend_G[step-1]+1)[0][0]
#                     index_end_change_in_path = np.where(path==step-1)[0][0]
#                     temp_min = np.min([index_begin_change_in_path,index_end_change_in_path])
#                     temp_max = np.max([index_begin_change_in_path,index_end_change_in_path])
#                     temp_path = path.copy()
#                     for i in range(temp_min,temp_max+1):
#                         arc_direction[path[i]] = -arc_direction[path[i]]
#                         temp_path[i] = path[temp_min+temp_max-i]
#                     path = temp_path
#                     step = oppend_G[step-1]+1
#
#                 else:#otherside is anticlockwise
#                     left_arcs_crossings[G_ind_matrix_to_ind_namelist[step-1],:] = np.array([step,oppend_G[step-1]],dtype=np.int64)
#
#                     index_begin_change_in_path = np.where(path==oppend_G[step-1])[0][0]
#                     index_end_change_in_path = np.where(path==step-1)[0][0]
#                     temp_min = np.min([index_begin_change_in_path,index_end_change_in_path])
#                     temp_max = np.max([index_begin_change_in_path,index_end_change_in_path])
#                     temp_path = path.copy()
#                     for i in range(temp_min,temp_max+1):
#                         arc_direction[path[i]] = -arc_direction[path[i]]
#                         temp_path[i] = path[temp_min+temp_max-i]
#                     path = temp_path
#                     step = oppend_G[step-1].copy()
#
#             else: #second time to the crossing
#                 right_arcs_crossings[G_ind_matrix_to_ind_namelist[step-1],:] = np.array([step,path[1+step_counter]],dtype=np.int64)
#
#                 step = path[1+step_counter]
#
#         step_counter += 1
#     both_sides_crossings = np.zeros(len_name_list,dtype=np.int64)#using index of namelist_G
#     for i in range(len_name_list):
#         if (left_arcs_crossings[i,1]-left_arcs_crossings[i,0])*(right_arcs_crossings[i,1]-right_arcs_crossings[i,0])<0:
#             both_sides_crossings[i] = 1
#
#
#     ############
#     ############the following part is for signed gauss diagram
#     epsilon_crossings = np.ones(len_name_list,dtype=np.int64)
#     above_down_crossings = np.zeros(len_name_list,dtype=np.int64)
#     for i in range(len_name_list):
#         temp_arr = np.append(left_arcs_crossings[i,:],right_arcs_crossings[i,:])
#         temp_arr_sorted = sorted(temp_arr)
#
#         if temp_arr[0] == temp_arr_sorted[2] and temp_arr[1] == temp_arr_sorted[0]:
#             epsilon_crossings[i] = -1
#         if temp_arr[0] == temp_arr_sorted[0] and temp_arr[1] == temp_arr_sorted[3]:
#             epsilon_crossings[i] = -1
#         if temp_arr[0] == temp_arr_sorted[3] and temp_arr[1] == temp_arr_sorted[1]:
#             epsilon_crossings[i] = -1
#         if temp_arr[0] == temp_arr_sorted[1] and temp_arr[1] == temp_arr_sorted[2]:
#             epsilon_crossings[i] = -1
#     for i in range(len_name_list):
#         if G[G_ind_namelist_to_ind_matrix[i,0],1]*G[G_ind_namelist_to_ind_matrix[i,0],2]==epsilon_crossings[i]:
#             above_down_crossings[i] = 1
#         else:
#             above_down_crossings[i] = -1
#     band_matrix = np.zeros([G.shape[0],4],dtype=np.int64)
#     band_matrix[:,0] = name_crossings
#     for i in range(G.shape[0]):
#         band_matrix[i,1] = both_sides_crossings[name_crossings[i]]
#         band_matrix[i,2] = above_down_crossings[name_crossings[i]]
#         band_matrix[i,3] = G[G_ind_namelist_to_ind_matrix[name_crossings[i],0],2]
#     return band_matrix
# =============================================================================


@njit(cache=True)
def Gauss_to_band(G):
    """
    G: Gauss diagram
    Output:
    """
    (G_name_list, G_ind_namelist_to_ind_matrix, G_ind_matrix_to_ind_namelist) = index_names_trans(G)
    oppend_G = np.zeros(G.shape[0], dtype=np.int64)
    len_name_list = len(G_name_list)
    oppend_G[G_ind_namelist_to_ind_matrix[:, 0]] = G_ind_namelist_to_ind_matrix[:, 1]
    oppend_G[G_ind_namelist_to_ind_matrix[:, 1]] = G_ind_namelist_to_ind_matrix[:, 0]
    crossings_check_list = np.zeros(len_name_list, dtype=np.int64)
    arc_direction = np.ones(G.shape[0] + 1, dtype=np.int64)
    path = np.arange(G.shape[0] + 1)  # written for the arcs
    name_crossings = np.zeros(
        G.shape[0], dtype=np.int64
    )  # using the namelist index to call crossings
    left_arcs_crossings = np.zeros(
        (len_name_list, 2), dtype=np.int64
    )  # using the namelist index to call crossings
    right_arcs_crossings = np.zeros(
        (len_name_list, 2), dtype=np.int64
    )  # using the namelist index to call crossings

    step = 0
    step_counter = 0
    while step != G.shape[0]:
        if arc_direction[step] == 1:  # anti_clockwise on the arc
            name_crossings[step_counter] = G_ind_matrix_to_ind_namelist[step]
            if (
                crossings_check_list[G_ind_matrix_to_ind_namelist[step]] == 0
            ):  # first time come to a crossing
                crossings_check_list[G_ind_matrix_to_ind_namelist[step]] = 1
                if arc_direction[oppend_G[step]] == 1:  # otherside is anticlockwise
                    left_arcs_crossings[G_ind_matrix_to_ind_namelist[step], :] = np.array(
                        [step, oppend_G[step]], dtype=np.int64
                    )

                    index_begin_change_in_path = np.where(path == oppend_G[step])[0][0]
                    index_end_change_in_path = np.where(path == step + 1)[0][0]
                    temp_min = np.min(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_max = np.max(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_path = path.copy()
                    for i in range(temp_min, temp_max + 1):
                        arc_direction[path[i]] = -arc_direction[path[i]]
                        temp_path[i] = path[temp_min + temp_max - i]
                    path = temp_path
                    step = oppend_G[step]

                else:  # otherside is clockwise
                    left_arcs_crossings[G_ind_matrix_to_ind_namelist[step], :] = np.array(
                        [step, 1 + oppend_G[step]], dtype=np.int64
                    )

                    index_begin_change_in_path = np.where(path == oppend_G[step] + 1)[0][0]
                    index_end_change_in_path = np.where(path == step + 1)[0][0]
                    temp_min = np.min(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_max = np.max(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_path = path.copy()
                    for i in range(temp_min, temp_max + 1):
                        arc_direction[path[i]] = -arc_direction[path[i]]
                        temp_path[i] = path[temp_min + temp_max - i]
                    path = temp_path
                    step = 1 + oppend_G[step]

            else:  # second time to the crossing
                right_arcs_crossings[G_ind_matrix_to_ind_namelist[step], :] = np.array(
                    [step, path[1 + step_counter]], dtype=np.int64
                )

                step = path[1 + step_counter]

        else:  # clockwise on the arc
            name_crossings[step_counter] = G_ind_matrix_to_ind_namelist[step - 1]
            if (
                crossings_check_list[G_ind_matrix_to_ind_namelist[step - 1]] == 0
            ):  # first time come to a crossing
                crossings_check_list[G_ind_matrix_to_ind_namelist[step - 1]] = 1
                if arc_direction[oppend_G[step - 1]] == -1:  # otherside is clockwise
                    left_arcs_crossings[G_ind_matrix_to_ind_namelist[step - 1], :] = np.array(
                        [step, oppend_G[step - 1] + 1], dtype=np.int64
                    )

                    index_begin_change_in_path = np.where(path == oppend_G[step - 1] + 1)[0][0]
                    index_end_change_in_path = np.where(path == step - 1)[0][0]
                    temp_min = np.min(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_max = np.max(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_path = path.copy()
                    for i in range(temp_min, temp_max + 1):
                        arc_direction[path[i]] = -arc_direction[path[i]]
                        temp_path[i] = path[temp_min + temp_max - i]
                    path = temp_path
                    step = oppend_G[step - 1] + 1

                else:  # otherside is anticlockwise
                    left_arcs_crossings[G_ind_matrix_to_ind_namelist[step - 1], :] = np.array(
                        [step, oppend_G[step - 1]], dtype=np.int64
                    )

                    index_begin_change_in_path = np.where(path == oppend_G[step - 1])[0][0]
                    index_end_change_in_path = np.where(path == step - 1)[0][0]
                    temp_min = np.min(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_max = np.max(
                        np.array([index_begin_change_in_path, index_end_change_in_path])
                    )
                    temp_path = path.copy()
                    for i in range(temp_min, temp_max + 1):
                        arc_direction[path[i]] = -arc_direction[path[i]]
                        temp_path[i] = path[temp_min + temp_max - i]
                    path = temp_path
                    step = oppend_G[step - 1]

            else:  # second time to the crossing
                right_arcs_crossings[G_ind_matrix_to_ind_namelist[step - 1], :] = np.array(
                    [step, path[1 + step_counter]], dtype=np.int64
                )

                step = path[1 + step_counter]

        step_counter += 1
    both_sides_crossings = np.zeros(len_name_list, dtype=np.int64)  # using index of namelist_G
    for i in range(len_name_list):
        if (left_arcs_crossings[i, 1] - left_arcs_crossings[i, 0]) * (
            right_arcs_crossings[i, 1] - right_arcs_crossings[i, 0]
        ) < 0:
            both_sides_crossings[i] = 1

    ############
    ############the following part is for signed gauss diagram
    epsilon_crossings = np.ones(len_name_list, dtype=np.int64)
    above_down_crossings = np.zeros(len_name_list, dtype=np.int64)
    for i in range(len_name_list):
        len_left_arcs = len(left_arcs_crossings[i, :])
        len_right_arcs = len(right_arcs_crossings[i, :])
        len_total = len_left_arcs + len_right_arcs

        temp_arr = np.zeros(len_total, dtype=np.int64)
        temp_arr[0:len_left_arcs] = left_arcs_crossings[i, :]
        temp_arr[len_left_arcs:len_total] = right_arcs_crossings[i, :]
        # temp_arr = np.append(left_arcs_crossings[i,:],right_arcs_crossings[i,:])
        # print('temp_arr=',temp_arr)
        temp_arr_sorted = sorted(temp_arr)

        if temp_arr[0] == temp_arr_sorted[2] and temp_arr[1] == temp_arr_sorted[0]:
            epsilon_crossings[i] = -1
        if temp_arr[0] == temp_arr_sorted[0] and temp_arr[1] == temp_arr_sorted[3]:
            epsilon_crossings[i] = -1
        if temp_arr[0] == temp_arr_sorted[3] and temp_arr[1] == temp_arr_sorted[1]:
            epsilon_crossings[i] = -1
        if temp_arr[0] == temp_arr_sorted[1] and temp_arr[1] == temp_arr_sorted[2]:
            epsilon_crossings[i] = -1

    for i in range(len_name_list):
        if (
            G[G_ind_namelist_to_ind_matrix[i, 0], 1] * G[G_ind_namelist_to_ind_matrix[i, 0], 2]
            == epsilon_crossings[i]
        ):
            above_down_crossings[i] = 1
        else:
            above_down_crossings[i] = -1
    band_matrix = np.zeros((G.shape[0], 4), dtype=np.int64)
    band_matrix[:, 0] = name_crossings
    band_matrix[:, 1] = both_sides_crossings[name_crossings]
    band_matrix[:, 2] = above_down_crossings[name_crossings]
    band_matrix[:, 3] = G[G_ind_namelist_to_ind_matrix[name_crossings, 0], 2]

    return band_matrix


# =============================================================================
# @njit(cache=True)
# def genus_zero_band(B):
#     """
#     B:band matrix
#     """
#     if np.any(B[:,1])==1:
#         return False
#     (B_name_list, B_ind_namelist_to_ind_matrix, B_ind_matrix_to_ind_namelist) = index_names_trans(B)
#     above_index_namelist = np.array([],dtype=np.int64)
#     below_index_namelist = np.array([],dtype=np.int64)
#     for i in range(len(B_name_list)):
#         if B[B_ind_namelist_to_ind_matrix[i,0],2] == 1:
#             above_index_namelist = np.append(above_index_namelist,[i])
#         else:
#             below_index_namelist = np.append(below_index_namelist,[i])
#     for i in range(len(above_index_namelist)-1):
#         for j in range(i,len(above_index_namelist)):
#             [a,b] = B_ind_namelist_to_ind_matrix[above_index_namelist[i],:]
#             [c,d] = B_ind_namelist_to_ind_matrix[above_index_namelist[j],:]
#             if (c-a)*(c-b)*(d-a)*(d-b)<0:
#                 return False
#     for i in range(len(below_index_namelist)-1):
#         for j in range(i,len(below_index_namelist)):
#             [a,b] = B_ind_namelist_to_ind_matrix[below_index_namelist[i],:]
#             [c,d] = B_ind_namelist_to_ind_matrix[below_index_namelist[j],:]
#             if (c-a)*(c-b)*(d-a)*(d-b)<0:
#                 return False
#     return True
# =============================================================================


@njit(cache=True)
def genus_zero_band(B):
    """
    B:band matrix
    """
    if np.any(B[:, 1]) == 1:
        return False
    (B_name_list, B_ind_namelist_to_ind_matrix, B_ind_matrix_to_ind_namelist) = index_names_trans(B)
    len_B_name_list = len(B_name_list)
    above_index_namelist = np.array(
        [i for i in range(len_B_name_list) if B[B_ind_namelist_to_ind_matrix[i, 0], 2] == 1]
    )
    below_index_namelist = np.array(
        [i for i in range(len_B_name_list) if B[B_ind_namelist_to_ind_matrix[i, 0], 2] == -1]
    )
    for i in range(len(above_index_namelist) - 1):
        for j in range(i, len(above_index_namelist)):
            a = B_ind_namelist_to_ind_matrix[above_index_namelist[i], 0]
            b = B_ind_namelist_to_ind_matrix[above_index_namelist[i], 1]
            # [a,b] = B_ind_namelist_to_ind_matrix[above_index_namelist[i],:]
            c = B_ind_namelist_to_ind_matrix[above_index_namelist[j], 0]
            d = B_ind_namelist_to_ind_matrix[above_index_namelist[j], 1]
            # [c,d] = B_ind_namelist_to_ind_matrix[above_index_namelist[j],:]
            if (c - a) * (c - b) * (d - a) * (d - b) < 0:
                return False
    for i in range(len(below_index_namelist) - 1):
        for j in range(i, len(below_index_namelist)):
            a = B_ind_namelist_to_ind_matrix[below_index_namelist[i], 0]
            b = B_ind_namelist_to_ind_matrix[below_index_namelist[i], 1]
            # [a,b] = B_ind_namelist_to_ind_matrix[below_index_namelist[i],:]
            c = B_ind_namelist_to_ind_matrix[below_index_namelist[j], 0]
            d = B_ind_namelist_to_ind_matrix[below_index_namelist[j], 1]
            # [c,d] = B_ind_namelist_to_ind_matrix[below_index_namelist[j],:]
            if (c - a) * (c - b) * (d - a) * (d - b) < 0:
                return False
    return True


def real_Gauss_signed(G):
    """
    G: Gauss diagram
    output: True if realisible. Otherwise False
    """
    B = Gauss_to_band(G)
    return genus_zero_band(B)


# =============================================================================
# G = np.array([
#     [1,1,1],
#     [2,-1,1],
#     [3,1,1],
#     [1,-1,1],
#     [2,1,1],
#     [3,-1,1]
# ],int)
# print(real_Gauss_signed(G))
# =============================================================================
