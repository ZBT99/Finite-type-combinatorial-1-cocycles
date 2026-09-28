#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:23:01 2024

@author: butianzhang
"""

import numpy as np

# import cupy as cp
import copy
from itertools import combinations
import math
from numba import njit
# from tqdm import tqdm

###################
#            BASICS


@njit(cache=True)
def index_first_appear(x, arr):
    """
    param: x: an elment in arr
    param: arr: an 1-dim arr
    output:index of x first appears in arr
    """
    for i in range(arr.shape[0]):
        if arr[i] == x:
            return i


@njit(cache=True)
def oppend(G, i):
    """
    :param G: Gauss diagram matrix
    :param i: a row index of G

    Output: another index of the same arrow
    """
    if i < 0 or i >= G.shape[0]:
        print("index required has pass the range of G")
        raise ValueError("No valid index found for oppend.")

    for j in range(G.shape[0]):
        # if G[j,0] == G[i,0] and G[j,1] == -G[i,1]:
        if G[j, 0] == G[i, 0] and j != i:
            return j
    raise ValueError("No valid index found for oppend.")


@njit(cache=True)
def left_right(G, index_along):
    """
    G: Gauss diagram matrix
    index_along: the index in G going along

    Outputs:[next left part index, next right part index]
    """
    row_G = G.shape[0]
    opp = oppend(G, index_along)
    if opp == -1:
        print("error in left_right()")
        raise ValueError("No valid index found for left_right.")

    if G[index_along, 2] == -1 and G[index_along, 1] == 1:
        return ((opp - 1) % row_G, (opp + 1) % row_G)
    elif G[index_along, 2] == -1 and G[index_along, 1] == -1:
        return ((opp + 1) % row_G, (opp - 1) % row_G)
    elif G[index_along, 2] == 1 and G[index_along, 1] == 1:
        return ((opp + 1) % row_G, (opp - 1) % row_G)
    elif G[index_along, 2] == 1 and G[index_along, 1] == -1:
        return ((opp - 1) % row_G, (opp + 1) % row_G)
    else:
        print("Invalid inputs")
        raise ValueError("No valid index found for left_right.")


@njit(cache=True)
def head_foot(G, n):
    A = np.zeros(2, dtype=np.int64)
    s = 0
    for i in range(0, G.shape[0]):
        if G[i, 0] == n:
            A[s] = i
            s += 1
    B = A.copy()
    if G[A[0], 1] == -1 and G[A[1], 1] == 1:
        A[0] = B[1]
        A[1] = B[0]

    #   if G[i,1] == 1:
    #     A[0] = i
    #   else:
    #     A[1] = i
    if np.any(A) == 0:
        print("Sorry, I do not find the the labeled arrow")
    # else:
    return A


@njit(cache=True)
def getname(G):
    row_G = G.shape[0]
    name_list = np.zeros(row_G, dtype=np.int64)
    name_list[0] = G[0, 0]
    if row_G == 1:
        return name_list[0:1]
    s = 1
    for i in range(1, row_G):
        if (G[i, 0] in name_list[0:s]) == 0:  # G[i,0] is new
            name_list[s] = G[i, 0]
            s = s + 1
    return name_list[0:s]


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


@njit(cache=True)
def normalize(G):
    """
    G: Gauss diagram
    return[G_normed,product of all (-1) in crossing of G]
    Remind: if a sign of arrow in G has 0, it does not matter.
    """
    N = G.copy()
    sign = 1
    name_list = np.zeros((G.shape[0]) // 2, dtype=np.int64)
    count_name = 0
    for i in range(G.shape[0]):
        appeared = False
        for j in range(count_name):
            if name_list[j] == G[i, 0]:
                appeared = True
                N[i, 0] = j + 1
                break
        if not appeared:
            name_list[count_name] = G[i, 0]
            N[i, 0] = count_name + 1
            if G[i, 2] == -1:
                sign = -sign
            count_name += 1
    return (N, sign)


@njit(cache=True)
def Gauss_matrix_to_PD(G):
    """
    Input: Gauss matrix with signed 3 cols
    ouput: PD code matrix of G
    Index of rows gives index of arcs, but we plus 1 for all of the names of arcs
    """

    (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(G)
    PD = np.zeros((name_list.shape[0], 4), dtype=np.int64)
    for i in range(name_list.shape[0]):
        if G[ind_namelist_to_ind_matrix[i, 0], 2] == 1:  # positive
            if G[ind_namelist_to_ind_matrix[i, 0], 1] == 1:  # first head
                PD[i, :] = [
                    ind_namelist_to_ind_matrix[i, 1],
                    ind_namelist_to_ind_matrix[i, 0] + 1,
                    ind_namelist_to_ind_matrix[i, 1] + 1,
                    ind_namelist_to_ind_matrix[i, 0],
                ]
            else:  # first foot
                PD[i, :] = [
                    ind_namelist_to_ind_matrix[i, 0],
                    ind_namelist_to_ind_matrix[i, 1] + 1,
                    ind_namelist_to_ind_matrix[i, 0] + 1,
                    ind_namelist_to_ind_matrix[i, 1],
                ]
        else:  # negative
            if G[ind_namelist_to_ind_matrix[i, 0], 1] == 1:  # first head
                PD[i, :] = [
                    ind_namelist_to_ind_matrix[i, 1],
                    ind_namelist_to_ind_matrix[i, 0],
                    ind_namelist_to_ind_matrix[i, 1] + 1,
                    ind_namelist_to_ind_matrix[i, 0] + 1,
                ]
            else:  # first foot
                PD[i, :] = [
                    ind_namelist_to_ind_matrix[i, 0],
                    ind_namelist_to_ind_matrix[i, 1],
                    ind_namelist_to_ind_matrix[i, 0] + 1,
                    ind_namelist_to_ind_matrix[i, 1] + 1,
                ]
    PD %= G.shape[0]
    PD += 1
    return PD


@njit(cache=True)
def examiner(G):
    """
    G: Guass diagram formula
    This function examins whethe G has errors.
    Output: True if no problems are found.
            False if some preblem has been found.
    """
    row_G = G.shape[0]
    col = G.shape[1]

    if row_G % 2 != 0:
        print("Gauss diagram has odd number of rows")
        return False
    if col != 3:
        print("Gauss diagram should have 3 colomns")
        return False

    name_list = getname(G)
    N_names = name_list.shape[0]
    name_index = -np.ones((N_names, 2), dtype=np.int64)
    name_count = np.zeros(N_names, dtype=np.int64)
    for i in range(row_G):
        for j in range(N_names):
            if G[i, 0] == name_list[j]:
                name_index[j, name_count[j] % 2] = i
                name_count[j] = name_count[j] + 1

    for i in range(N_names):
        if name_count[i] != 2:
            print("Gauss diagram matrix is wrong")
            return False
        if G[name_index[i, 0], 1] + G[name_index[i, 1], 1] != 0:
            print("Gauss diagram matrix is wrong")
            return False
        if G[name_index[i, 0], 2] != G[name_index[i, 1], 2]:
            print("Gauss diagram matrix is wrong")
            return False

    return True


@njit(cache=True)
def exam_match_Gauss_matrix(G, A, exam):
    if exam:
        if not examiner(G):
            print("the first diagram is wrong")
            return
        ##0 rows need to be treated in A.
        ## WARNING: This is one of the reason why the name '0' is forbidden for any arrow.
        row_B = A.shape[0]
        for i in range(0, int(A.shape[0])):
            if A[i, 0] == 0:
                row_B = i
                break
        B = A[0:row_B, :]
        if not examiner(B):
            print("the second diagram is wrong")
            return
    else:  # do not exam but treat A
        row_B = A.shape[0]
        for i in range(0, int(A.shape[0])):
            if A[i, 0] == 0:
                row_B = i
                break
        B = A[0:row_B, :]
        if not examiner(B):
            print("the second diagram is wrong")
            return
    return B


@njit(cache=True)
def unsigned_in_combination(G, B, ind_namelist_to_ind_matrix, c):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)

    (D, sign) = normalize(G[index, 0:3])
    D[:, 2] = 0
    if np.any(D - B) == 0:  # D==B
        return sign
    else:
        return 0


@njit(cache=True)
def signed_in_combination(G, B, ind_namelist_to_ind_matrix, c):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)

    (D, _) = normalize(G[index, 0:3])
    if np.any(D - B) == 0:  # D==B
        return 1
    else:
        return 0


def match_Gauss_matrix(G, A, exam=True):
    B = exam_match_Gauss_matrix(G, A, exam)

    ar = B.shape[0] // 2
    cr = G.shape[0] // 2
    if cr < ar:
        return 0
    # now cr >= ar
    result = 0
    (name_list, ind_namelist_to_ind_matrix, _) = index_names_trans(G)
    if B[0, 2] == 0:  # unsigned
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = unsigned_in_combination(G, B, ind_namelist_to_ind_matrix, c)
            result += contr
    if B[0, 2] != 0:  # signed
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = signed_in_combination(G, B, ind_namelist_to_ind_matrix, c)
            result += contr
    return result


@njit(cache=True)
def unsigned_in_combination_conditional(
    G, B, ind_namelist_to_ind_matrix, c, first_foot, last_foot_lowerbounds, last_foot_upperbounds
):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)
    if index[0] == first_foot and last_foot_lowerbounds <= index[-1] <= last_foot_upperbounds:
        (D, sign) = normalize(G[index, 0:3])
        D[:, 2] = 0
        if np.any(D - B) == 0:  # D==B
            return sign
        else:
            return 0
    else:
        return 0


@njit(cache=True)
def signed_in_combination_conditional(
    G, B, ind_namelist_to_ind_matrix, c, first_foot, last_foot_lowerbounds, last_foot_upperbounds
):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)
    if index[0] == first_foot and last_foot_lowerbounds <= index[-1] <= last_foot_upperbounds:
        (D, _) = normalize(G[index, 0:3])
        if np.any(D - B) == 0:  # D==B
            return 1
        else:
            return 0
    else:
        return 0


def match_Gauss_matrix_conditional(
    G, A, first_foot, last_foot_lowerbounds=0, last_foot_upperbounds=np.inf, exam=True
):
    """
    This function is desiged for matching A in G with some conditions.
    G: Gauss diagram matrix
    A: subdiagram to be found in G
    first_foot: index of first foot of A in G
    last_foot_lowerbounds: index of last_foot in G >= last_foot_lowerbounds
    last_foot_upperbounds: index of last_foot in G <= last_foot_upperbounds
    """
    B = exam_match_Gauss_matrix(G, A, exam)

    ar = B.shape[0] // 2
    cr = G.shape[0] // 2
    if cr < ar:
        return 0
    # now cr >= ar
    result = 0
    (name_list, ind_namelist_to_ind_matrix, _) = index_names_trans(G)
    if B[0, 2] == 0:  # unsigned
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = unsigned_in_combination_conditional(
                G,
                B,
                ind_namelist_to_ind_matrix,
                c,
                first_foot,
                last_foot_lowerbounds,
                last_foot_upperbounds,
            )
            result += contr

    if B[0, 2] != 0:  # signed
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = signed_in_combination_conditional(
                G,
                B,
                ind_namelist_to_ind_matrix,
                c,
                first_foot,
                last_foot_lowerbounds,
                last_foot_upperbounds,
            )
            result += contr
    return result


@njit(cache=True)
def dec_to_2(x):
    """
    x: a positive integer
    output: 1-d array of numbers 0, 1 representing x in base-2.
    """
    if x == 0:
        return np.array([0], dtype=np.int64)

    # Use a Python list to collect base-2 digits
    temp = x
    temp_r = []

    # Convert to base-2
    while temp != 0:
        temp_r.append(temp % 2)
        temp = temp // 2

    # Convert the result to a NumPy array
    return np.array(temp_r, dtype=np.int64)


@njit(cache=True)
def dec_to_3(x):
    """
    x: a positive integer
    output: 1-d array of numbers 0, 1, 2 representing x in base-3.
    """
    if x == 0:
        return np.array([0], dtype=np.int64)

    # Use a Python list to collect base-3 digits
    temp = x
    temp_r = []

    # Convert to base-3
    while temp != 0:
        temp_r.append(temp % 3)
        temp = temp // 3

    # Convert the result to a NumPy array
    return np.array(temp_r, dtype=np.int64)


def numb_dia(N):
    return math.factorial(2 * N) // (math.factorial(N) * (2**N))


def diagram_underlying(N, connected=False):
    A = np.zeros([2 * N, 3, numb_dia(N)], dtype=np.int64)
    count = 0
    if N == 0:
        return (A, 1)
    if N == 1:
        A[:, :, 0] = np.array([[1, 0, 0], [1, 0, 0]], int)
        return (A, 1)
    (dia_pre, numb_pre) = diagram_underlying(N - 1, connected)
    for i in range(numb_pre):
        [name_list, ind_namelist_to_ind_matrix, ind_mat_to_ind_name_list] = index_names_trans(
            dia_pre[:, :, i]
        )
        appearance = np.zeros(name_list.shape[0])
        for j in range(2 * (N - 1)):
            appearance[ind_mat_to_ind_name_list[j]] += 1
            if np.all(appearance) == 1:
                index_pre_last = j
                break
        # print('index_pre_last=',index_pre_last)

        for c in combinations(np.arange(index_pre_last + 1, 2 * N - 1 + 1), 2):
            d = np.zeros(2, dtype=np.int64)
            d[0] = c[0]
            d[1] = c[1] - 1
            # print('d=',d)
            if not connected:
                arr = np.array([N, 0, 0], int)
                A[:, :, count] = np.insert(dia_pre[:, :, i], d, [arr, arr], axis=0)
                count += 1
            if connected:
                if d[0] > np.max(ind_namelist_to_ind_matrix[:, 0]) and d[0] <= np.max(
                    ind_namelist_to_ind_matrix[:, 1]
                ):
                    arr = np.array([N, 0, 0], int)
                    A[:, :, count] = np.insert(dia_pre[:, :, i], d, [arr, arr], axis=0)
                    count += 1
    return (A[:, :, 0:count], count)


@njit(cache=True)
def main_diagram_generate_param(N, A, numb_chord):
    numb_directed = numb_chord * 2**N
    numb_directed_signed = numb_directed * 2**N
    B = np.zeros((2 * N, 3, numb_directed), dtype=np.int64)
    C = np.zeros((2 * N, 3, numb_directed_signed), dtype=np.int64)
    # count_chord = 0
    count_directed = 0
    count_directed_signed = 0
    for i in range(numb_chord):
        (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(
            A[:, :, i]
        )
        for j in range(2**N):
            ind_directed0 = dec_to_2(j)
            ind_directed = np.zeros(N, dtype=np.int64)
            ind_directed[0 : ind_directed0.shape[0]] = ind_directed0
            B[:, :, count_directed] = A[:, :, i]
            for ind in range(N):
                ind_directed[ind] = 2 * ind_directed[ind] - 1
                # print('ind_namelist_to_ind_matrix[ind,:]=',ind_namelist_to_ind_matrix[ind,:])
                B[ind_namelist_to_ind_matrix[ind, 0], 1, count_directed] = ind_directed[ind]
                B[ind_namelist_to_ind_matrix[ind, 1], 1, count_directed] = -ind_directed[ind]
            # print('B[:,:,count_directed]')
            # print(B[:,:,count_directed])
            # draw(B[:,:,count_directed],count_directed)
            for k in range(2**N):
                ind_directed_signed0 = dec_to_2(k)
                ind_directed_signed = np.zeros(N, dtype=np.int64)
                ind_directed_signed[0 : ind_directed_signed0.shape[0]] = ind_directed_signed0
                C[:, :, count_directed_signed] = B[:, :, count_directed]
                for ind in range(N):
                    ind_directed_signed[ind] = 2 * ind_directed_signed[ind] - 1
                    C[ind_namelist_to_ind_matrix[ind, 0], 2, count_directed_signed] = (
                        ind_directed_signed[ind]
                    )
                    C[ind_namelist_to_ind_matrix[ind, 1], 2, count_directed_signed] = (
                        ind_directed_signed[ind]
                    )

                count_directed_signed += 1
            count_directed += 1
    return (A, B, C, numb_chord, count_directed, count_directed_signed)


def diagram_generate_param(N, connected=False):
    (A, numb_chord) = diagram_underlying(N, connected)
    return main_diagram_generate_param(N, A, numb_chord)


###################
#          DRAWINGS


@njit(cache=True)
def signstr(x):
    if x > 0:
        return "+"
    if x == 0:
        return ""
    if x < 0:
        return "\u2212"  # ← Unicode minus sign


# sign position correction
@njit(cache=True)
def sp_correction(x):
    if x == 0:
        return 0  # 1
    else:
        return 0.16  # 1.15


@njit(cache=True)
def draw_prepare(Y):
    X = Y.copy()
    for i in range(0, int(X.shape[0])):
        if X[i, 0] == 0:
            X = X[0:i, :]
            break
    N = X.shape[0] // 2
    (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(X)
    head_foot_under_index_namelist = np.zeros((N, 2), dtype=np.int64)
    for i in range(N):
        ind = ind_namelist_to_ind_matrix[i]
        if X[ind[0], 1] == -1 and X[ind[1], 1] == 1:
            head_foot_under_index_namelist[i] = np.array([ind[1], ind[0]], dtype=np.int64)
        elif X[ind[0], 1] == 1 and X[ind[1], 1] == -1:
            head_foot_under_index_namelist[i] = np.array([ind[0], ind[1]], dtype=np.int64)
        else:
            head_foot_under_index_namelist[i] = ind

    x = np.zeros(2 * N)
    y = np.zeros(2 * N)
    theta = np.zeros(2 * N)
    for i in range(2 * N):
        theta[i] = -np.pi / 2 + np.pi / (2 * N) + i * np.pi / N
        x[i] = np.cos(theta[i])
        y[i] = np.sin(theta[i])
    return (head_foot_under_index_namelist, X, theta, x, y, N)


def draw(
    Y,
    tit="",
    color1_ind=np.array([], int),
    color2_ind=np.array([], int),
    color3_ind=np.array([], int),
):
    """
    Y:the Gauss diagram matrix to be draw
    WARNING: We require that names in Y must be strictly positive.
    """
    from matplotlib import pyplot as plt

    a = 0.67
    b = 0.36

    ##a,b for drawing foot with string on the circle....

    (head_foot_under_index_namelist, X, theta, x, y, N) = draw_prepare(Y)

    plt.figure(figsize=(4, 4))
    plt.xlim((-1.5, 1.5))
    plt.ylim((-1.5, 1.5))

    for i in range(N):
        ind_head = head_foot_under_index_namelist[i, 0]
        ind_foot = head_foot_under_index_namelist[i, 1]
        dx = -x[ind_head] + x[ind_foot]
        dy = -y[ind_head] + y[ind_foot]
        dr = np.sqrt(dx**2 + dy**2)

        if X[ind_head, 1] == 0 and X[ind_foot, 1] == 0:
            if (ind_head in color1_ind) == 1 or (ind_foot in color1_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-", color="red", connectionstyle="arc3"),
                )
            elif (ind_head in color2_ind) == 1 or (ind_foot in color2_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-", color="orange", connectionstyle="arc3"),
                )
            elif (ind_head in color3_ind) == 1 or (ind_foot in color3_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-", color="blue", connectionstyle="arc3"),
                )
            else:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-", connectionstyle="arc3"),
                )

        elif X[ind_head, 1] == 1 and X[ind_foot, 1] == -1:
            if (ind_head in color1_ind) == 1 or (ind_foot in color1_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-|>", color="red", connectionstyle="arc3,rad=0"),
                )

            elif (ind_head in color2_ind) == 1 or (ind_foot in color2_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-|>", color="orange", connectionstyle="arc3,rad=0"),
                )
            elif (ind_head in color3_ind) == 1 or (ind_foot in color3_ind) == 1:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-|>", color="blue", connectionstyle="arc3,rad=0"),
                )
            else:
                plt.annotate(
                    signstr(X[ind_foot, 2]),
                    xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
                    xytext=(
                        x[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dx / dr),
                        y[ind_foot]
                        + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                        * sp_correction(X[ind_foot, 2])
                        * (dy / dr),
                    ),  ##arrow foot coordinates
                    arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3,rad=0"),
                )

    from matplotlib import pyplot as plt

    circle = plt.Circle((0, 0), radius=1, color="black", fill=False)
    plt.gcf().gca().add_artist(circle)
    plt.plot(0, -1, "o", color="b")  # basepoint
    plt.title(tit)
    plt.show()
    return


def draw_lr_cfg(ax, G, R3_r_l_type):
    """
    G: a 4 colomns cfg matrix
    R3_r_l_type: lie in {0,1,2,3,4,5} for [l,h,m],[h,m,l],[m,l,h], [l,m,h],[m,h,l],[h,l,m].
    """
    from matplotlib import pyplot as plt

    circle = plt.Circle((0, 0), 1, color="black", fill=False)
    ax.add_artist(circle)
    ax.plot([0], [-1], "bo", markersize=5)
    R3_pt = np.zeros((3, 2), dtype=np.float64)  # [l, m, h]
    phi_0 = -np.pi / 6  # first piece of R3
    phi_1 = np.pi / 2  # second piece of R3
    phi_2 = 7 * np.pi / 6  # third piece of R3
    if R3_r_l_type == 0:
        R3_pt[0] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
    elif R3_r_l_type == 1:
        R3_pt[0] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
    elif R3_r_l_type == 2:
        R3_pt[0] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
    elif R3_r_l_type == 3:
        R3_pt[0] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
    elif R3_r_l_type == 4:
        R3_pt[0] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
    elif R3_r_l_type == 5:
        R3_pt[0] = np.array([np.cos(phi_1), np.sin(phi_1)], dtype=np.float64)
        R3_pt[1] = np.array([np.cos(phi_2), np.sin(phi_2)], dtype=np.float64)
        R3_pt[2] = np.array([np.cos(phi_0), np.sin(phi_0)], dtype=np.float64)
    rho = 0.95
    ax.arrow(
        R3_pt[0, 0],
        R3_pt[0, 1],
        (R3_pt[1, 0] - R3_pt[0, 0]) * rho,
        (R3_pt[1, 1] - R3_pt[0, 1]) * rho,
        head_width=0.05,
        head_length=0.1,
        fc="blue",
        ec="black",
    )
    ax.arrow(
        R3_pt[0, 0],
        R3_pt[0, 1],
        (R3_pt[2, 0] - R3_pt[0, 0]) * rho,
        (R3_pt[2, 1] - R3_pt[0, 1]) * rho,
        head_width=0.05,
        head_length=0.1,
        fc="blue",
        ec="black",
    )
    ax.arrow(
        R3_pt[1, 0],
        R3_pt[1, 1],
        (R3_pt[2, 0] - R3_pt[1, 0]) * rho,
        (R3_pt[2, 1] - R3_pt[1, 1]) * rho,
        head_width=0.05,
        head_length=0.1,
        fc="blue",
        ec="black",
    )

    (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(G)
    len_name_list = len(name_list)
    head_foot_under_index_namelist = np.zeros((len_name_list, 2), dtype=np.int64)
    for i in range(len_name_list):
        ind = ind_namelist_to_ind_matrix[i]
        if G[ind[0], 1] == -1:
            head_foot_under_index_namelist[i] = np.array([ind[1], ind[0]])
        else:
            head_foot_under_index_namelist[i] = ind.copy()

    count = np.zeros(4, dtype=np.int64)
    for i in range(G.shape[0]):
        count[G[i, 3]] += 1
    x = np.zeros(G.shape[0], dtype=np.float64)
    y = np.zeros(G.shape[0], dtype=np.float64)
    theta = np.zeros(G.shape[0], dtype=np.float64)
    for i in range(G.shape[0]):
        sum_0 = count[0]
        sum_1 = sum_0 + count[1]
        sum_2 = sum_1 + count[2]
        # sum_3 = sum_2 + count[3]

        if i < sum_0:
            # print('piece 0')
            theta[i] = -np.pi / 2 + np.pi / (3 * (count[0] + 1)) * (i + 1)
            x[i] = np.cos(theta[i])
            y[i] = np.sin(theta[i])
        elif i < sum_1:
            # print('piece 1')
            theta[i] = -np.pi / 6 + 2 * np.pi / (3 * (count[1] + 1)) * (i - sum_0 + 1)
            x[i] = np.cos(theta[i])
            y[i] = np.sin(theta[i])
        elif i < sum_2:
            # print('piece 2')
            theta[i] = np.pi / 2 + 2 * np.pi / (3 * (count[2] + 1)) * (i - sum_1 + 1)
            x[i] = np.cos(theta[i])
            y[i] = np.sin(theta[i])
        else:
            # print('piece 3')
            theta[i] = 7 * np.pi / 6 + np.pi / (3 * (count[3] + 1)) * (i - sum_2 + 1)
            x[i] = np.cos(theta[i])
            y[i] = np.sin(theta[i])
    # print(theta / np.pi * 180)
    # print('x=',x)
    # print('y=',y)
    from matplotlib import pyplot as plt

    a = 0.67
    b = 1

    for i in range(len_name_list):
        # print('i=',i )
        ind_head = head_foot_under_index_namelist[i, 0]
        ind_foot = head_foot_under_index_namelist[i, 1]
        dx = -x[ind_head] + x[ind_foot]
        dy = -y[ind_head] + y[ind_foot]
        dr = np.sqrt(dx**2 + dy**2)
        # print(x[ind_head],y[ind_head])
        ax.annotate(
            signstr(G[ind_foot, 2]),
            xy=(x[ind_head], y[ind_head]),  ##arrow head coordinates
            xytext=(
                x[ind_foot]
                + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                * sp_correction(G[ind_foot, 2])
                * (dx / dr),
                y[ind_foot]
                + (a + b * np.cos(theta[ind_foot] + 3 * np.pi / 4))
                * sp_correction(G[ind_foot, 2])
                * (dy / dr),
                # x[ind_foot]+ (dx/dr),
                # y[ind_foot]+ (dy/dr)
            ),  ##arrow foot coordinates
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3"),
        )

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_aspect("equal")
    ax.set_title("R3_r_l_type = {}".format(R3_r_l_type))
    ax.set_xticks([])  # 隐藏 x 轴的刻度
    ax.set_yticks([])  # 隐藏 y 轴的刻度
    ax.set_xticklabels([])  # 隐藏 x 轴的刻度标签
    ax.set_yticklabels([])  # 隐藏 y 轴的刻度标签


###################
#            KNOT OPERATIONS


@njit(cache=True)
def reverse_knot(G):
    G0 = G[::-1, :].copy()
    G0[:, 1] = -G0[:, 1]

    G1 = normalize(G0)[0]
    return G1


@njit(cache=True)
def mirror(G):
    """
    :param :G a knot diagram matrix
    Output: the mirrow image of the knot
    """
    G0 = -G
    G0[:, 0] = G[:, 0]
    return G0


@njit(cache=True)
def infty_avancer(G, n):
    """
    param: G: Gauss diagram matrix of a knot
    param: n: move infty n times counter-clockwisely

    Output: G_result

    """
    row_G = G.shape[0]
    G_result = G.copy()
    for i in range(row_G):
        G_result[i, :] = G[(i + n) % row_G, :]
    return normalize(G_result)[0]


@njit(cache=True)
def sum_knot(G0, G1):
    """
    :param G0: Gauss diagram matrix of a knot
    :param G1: Gauss diagram matrix of another knot
    Output: the long knot first going along G0 and then G1

    Warning: sum_knot can make mistakes for GDF with 0 names!!!
    """
    max_name = np.max(getname(G0))
    # row1 = G1.shape[0]
    G2 = G1.copy()
    G2[:, 0] = G1[:, 0] + max_name
    # for i in range(row1):
    #     G2[i,0] = G1[i,0] + max_name
    return np.append(G0, G2, axis=0)


curl = np.zeros([2, 3, 4], int)
curl[:, :, 0] = np.array([[1, 1, 1], [1, -1, 1]], int)
curl[:, :, 1] = np.array([[1, -1, 1], [1, 1, 1]], int)
curl[:, :, 2] = np.array([[1, 1, -1], [1, -1, -1]], int)
curl[:, :, 3] = np.array([[1, -1, -1], [1, 1, -1]], int)


@njit(cache=True)
def addcurl(G, c_type, n):
    """
    :param G: Gauss diagram formula
    :param c_type: type of curl should be 0,1,2,3
    :param n:add n curls, should be a natural number
    """
    G0 = G.copy()
    for i in range(n):
        G0 = sum_knot(G0, curl[:, :, c_type])
    return G0


@njit(cache=True)
def whitney(G):
    row = G.shape[0]
    result = 0
    names = getname(G)
    for i in range(int(row / 2)):
        HF = head_foot(G, names[i])
        if HF[0] - HF[1] < 0:  # head first
            if G[HF[0], 2] > 0:  # sign
                result += 1
            if G[HF[0], 2] < 0:
                result -= 1
        if HF[0] - HF[1] > 0:  # foot first
            if G[HF[0], 2] > 0:  # sign
                result -= 1
            if G[HF[0], 2] < 0:
                result += 1
    return result


@njit(cache=True)
def writhe(G):
    return np.sum(G[:, 2]) // 2


###################
#      REIDEMEISTER MOVE III

R3_table = np.zeros([2, 3, 3, 32], int)
# The R3_table is in the following order:
# Global type: left(0-15)/right(16-31)
# local type: 1-8 (0-7 in index). Local type = up_int[(index+1)/2]
# It might be better to consider the ordering number in 16-base number system.

#################
################# Left types

##########
R3_table[:, :, 0, 0] = np.array([[1, -1, 1], [2, -1, 1]], int)
R3_table[:, :, 1, 0] = np.array([[3, 1, 1], [1, 1, 1]], int)
R3_table[:, :, 2, 0] = np.array([[3, -1, 1], [2, 1, 1]], int)

###########
R3_table[:, :, 0, 2] = np.array([[1, -1, -1], [2, -1, 1]], int)
R3_table[:, :, 1, 2] = np.array([[2, 1, 1], [3, 1, -1]], int)
R3_table[:, :, 2, 2] = np.array([[3, -1, -1], [1, 1, -1]], int)

###########
R3_table[:, :, 0, 4] = np.array([[1, -1, 1], [2, -1, -1]], int)
R3_table[:, :, 1, 4] = np.array([[3, 1, -1], [2, 1, -1]], int)
R3_table[:, :, 2, 4] = np.array([[1, 1, 1], [3, -1, -1]], int)

###########
R3_table[:, :, 0, 6] = np.array([[1, -1, -1], [2, -1, -1]], int)
R3_table[:, :, 1, 6] = np.array([[1, 1, -1], [3, 1, 1]], int)
R3_table[:, :, 2, 6] = np.array([[2, 1, -1], [3, -1, 1]], int)
####
####
for i in range(4):
    temp = np.zeros([6, 3], int)

    temp[0, :] = R3_table[1, :, 0, 2 * i]
    temp[1, :] = R3_table[0, :, 0, 2 * i]

    temp[2, :] = R3_table[1, :, 1, 2 * i]
    temp[3, :] = R3_table[0, :, 1, 2 * i]

    temp[4, :] = R3_table[1, :, 2, 2 * i]
    temp[5, :] = R3_table[0, :, 2, 2 * i]

    temp = normalize(temp)[0]

    R3_table[:, :, 0, 2 * i + 1] = temp[0:2, :]
    R3_table[:, :, 1, 2 * i + 1] = temp[2:4, :]
    R3_table[:, :, 2, 2 * i + 1] = temp[4:6, :]
####
####
for i in range(2, 6):  # for case 6 and 7
    for j in range(3):
        R3_table[:, :, j, i + 8] = R3_table[:, :, j, i]
        for k in range(2):
            R3_table[k, 2, j, i + 8] = -R3_table[k, 2, j, i]

for i in range(0, 2):  # for case 8
    for j in range(3):
        R3_table[:, :, j, i + 14] = R3_table[:, :, j, i]
        for k in range(2):
            R3_table[k, 2, j, i + 14] = -R3_table[k, 2, j, i]

for i in range(6, 8):  # for case 5
    for j in range(3):
        R3_table[:, :, j, i + 2] = R3_table[:, :, j, i]
        for k in range(2):
            R3_table[k, 2, j, i + 2] = -R3_table[k, 2, j, i]

##########
##########
########## For the RIGHT types:

for i in range(0, 8):
    temp0 = np.zeros([6, 3], int)
    temp1 = np.zeros([6, 3], int)

    temp0[0:2, :] = R3_table[:, :, 0, 2 * i + 1]
    temp0[2:4, :] = R3_table[:, :, 2, 2 * i + 1]
    temp0[4:6, :] = R3_table[:, :, 1, 2 * i + 1]

    temp0 = normalize(temp0)[0]

    temp1[0:2, :] = R3_table[:, :, 0, 2 * i]
    temp1[2:4, :] = R3_table[:, :, 2, 2 * i]
    temp1[4:6, :] = R3_table[:, :, 1, 2 * i]

    temp1 = normalize(temp1)[0]

    R3_table[:, :, 0, 2 * i + 16] = temp0[0:2]
    R3_table[:, :, 1, 2 * i + 16] = temp0[2:4]
    R3_table[:, :, 2, 2 * i + 16] = temp0[4:6]

    R3_table[:, :, 0, 2 * i + 17] = temp1[0:2]
    R3_table[:, :, 1, 2 * i + 17] = temp1[2:4]
    R3_table[:, :, 2, 2 * i + 17] = temp1[4:6]

R3_table_matrix = np.zeros([6, 3, 32], int)
for i in range(32):
    R3_table_matrix[0:2, :, i] = R3_table[:, :, 0, i]
    R3_table_matrix[2:4, :, i] = R3_table[:, :, 1, i]
    R3_table_matrix[4:6, :, i] = R3_table[:, :, 2, i]


R3_dhml = np.zeros([3, 2, 32], int)
# row_0: indice of d in order [head,foot]
# row_1: indice of hm
# row_2: indice of ml
for i in range(32):
    if i == 0 or i == 14:
        R3_dhml[:, :, i] = np.array([[3, 0], [2, 4], [5, 1]], int)
    elif i == 1 or i == 15:
        R3_dhml[:, :, i] = np.array([[2, 1], [3, 5], [4, 0]], int)

    elif i == 3 or i == 11:
        R3_dhml[:, :, i] = np.array([[3, 0], [2, 5], [4, 1]], int)
    elif i == 2 or i == 10:
        R3_dhml[:, :, i] = np.array([[2, 1], [3, 4], [5, 0]], int)

    elif i == 5 or i == 13:
        R3_dhml[:, :, i] = np.array([[2, 0], [3, 4], [5, 1]], int)
    elif i == 4 or i == 12:
        R3_dhml[:, :, i] = np.array([[3, 1], [2, 5], [4, 0]], int)

    elif i == 6 or i == 8:
        R3_dhml[:, :, i] = np.array([[2, 0], [3, 5], [4, 1]], int)
    elif i == 7 or i == 9:
        R3_dhml[:, :, i] = np.array([[3, 1], [2, 4], [5, 0]], int)

    elif i == 17 or i == 31:
        R3_dhml[:, :, i] = np.array([[5, 0], [4, 2], [3, 1]], int)
    elif i == 16 or i == 30:
        R3_dhml[:, :, i] = np.array([[4, 1], [5, 3], [2, 0]], int)

    elif i == 18 or i == 26:
        R3_dhml[:, :, i] = np.array([[5, 0], [4, 3], [2, 1]], int)
    elif i == 19 or i == 27:
        R3_dhml[:, :, i] = np.array([[4, 1], [5, 2], [3, 0]], int)

    elif i == 20 or i == 28:
        R3_dhml[:, :, i] = np.array([[4, 0], [5, 2], [3, 1]], int)
    elif i == 21 or i == 29:
        R3_dhml[:, :, i] = np.array([[5, 1], [4, 3], [2, 0]], int)

    elif i == 23 or i == 25:
        R3_dhml[:, :, i] = np.array([[4, 0], [5, 3], [2, 1]], int)
    elif i == 22 or i == 24:
        R3_dhml[:, :, i] = np.array([[5, 1], [4, 2], [3, 0]], int)


@njit(cache=True)
def R3(G, arr):
    """
    G:Gauss diagram matrix
    arr = [i,j,k] being numpy or list
    i:the first index of the two consecutive crossings in the first piece
    j:the first index of the two consecutive crossings in the second piece
    k:the first index of the two consecutive crossings in the third piece
    WARNING: i<j<k

    Outputs: [G_after_R3,R3_type_index, tracked index, dhml]
    """
    i = arr[0]
    j = arr[1]
    k = arr[2]
    row_G = G.shape[0]
    if i < j - 1 < k - 2 < row_G - 3:
        row_G = row_G
    else:
        print("Index for R3 is not valide")
        return

    R3_ind = np.array([i, i + 1, j, j + 1, k, k + 1], dtype=np.int64)

    ####examinator for R3 moves,8 possibilities
    extract = np.zeros((6, 3), dtype=np.int64)
    temp_G = G[:, 0:3]
    if G[i, 1] == -1 and G[i + 1, 1] == -1:
        extract[0:2, :] = temp_G[i : i + 2, :]
        extract[2:4, :] = temp_G[j : j + 2, :]
        extract[4:6, :] = temp_G[k : k + 2, :]
        rot = 0
    elif G[j, 1] == -1 and G[j + 1, 1] == -1:
        extract[0:2, :] = temp_G[j : j + 2, :]
        extract[2:4, :] = temp_G[k : k + 2, :]
        extract[4:6, :] = temp_G[i : i + 2, :]
        rot = 1
    elif G[k, 1] == -1 and G[k + 1, 1] == -1:
        extract[0:2, :] = temp_G[k : k + 2, :]
        extract[2:4, :] = temp_G[i : i + 2, :]
        extract[4:6, :] = temp_G[j : j + 2, :]
        rot = 2
    else:
        print("Index required cannot do R3 ----1")
        return
    # print('extract',extract)
    if not examiner(extract):
        print("Index required cannot do R3 ----2")
        return

    extract = normalize(extract)[0]
    # print('after normalization', extract)
    for s in range(32):
        # any addition; all multiplication
        if not np.any(extract - R3_table_matrix[:, :, s]):  # addition result = 0
            R3_typeindex = s
            break
        else:
            if s == 31:
                print("Index required cannot do R3 ----3")
                return
    # print('R3_typeindex=',R3_typeindex)
    # print('rot=',rot)
    # print('d=',R3_dhml[0,:,R3_typeindex])
    # print('hm=',R3_dhml[1,:,R3_typeindex])
    # print('ml=',R3_dhml[2,:,R3_typeindex])
    dhml = np.zeros((3, 2), dtype=np.int64)
    for a in range(3):
        for b in range(2):
            dhml[a, b] = R3_ind[(R3_dhml[a, b, R3_typeindex] + 2 * rot) % 6]
    if R3_typeindex % 2 == 0:
        pass
    else:
        pass
    # print('d=',dhml[0,:])
    # print('hm=',dhml[1,:])
    # print('ml=',dhml[2,:])
    ####examination finished

    # draw(G,'to be R3ed',np.array([i,i+1,j,j+1,k,k+1],int))
    ##########Contribution part(calculations)

    ################Deformation part
    G0 = G.copy()

    G0[i + 1, :] = G[i, :]
    G0[i, :] = G[i + 1, :]

    G0[j + 1, :] = G[j, :]
    G0[j, :] = G[j + 1, :]

    G0[k + 1, :] = G[k, :]
    G0[k, :] = G[k + 1, :]
    G0 = normalize(G0)[0]
    index_before_to_later = np.arange(0, row_G)

    index_before_to_later[i] = i + 1
    index_before_to_later[i + 1] = i
    index_before_to_later[j] = j + 1
    index_before_to_later[j + 1] = j
    index_before_to_later[k] = k + 1
    index_before_to_later[k + 1] = k
    # draw(G0,'after R3',np.array([index_before_to_later[i],index_before_to_later[i+1],index_before_to_later[j],index_before_to_later[j+1],index_before_to_later[k],index_before_to_later[k+1]],int))

    # dhml_new = np.zeros([3,2],int)
    # for a in range(3):
    # for b in range(2):
    # dhml_new[a,b] = index_before_to_later[dhml[a,b]]

    return (G0, R3_typeindex, index_before_to_later, dhml)


##################
####              CONFIGURATIONS
class Configuration:
    def __init__(self, Gauss_matrix, R3_type_global):
        """
        Gauss_matrix:[-,4] the matrix of the configuration
        R3_type_global is the global matrix of R3
        """
        self.Gauss_matrix = Gauss_matrix
        self.R3_type_global = R3_type_global

        if np.any(Gauss_matrix[:, 2]) == 0:
            self.signed = False
        else:
            self.signed = True

    def exam(self):
        if self.Gauss_matrix.shape[0] == 0:
            return True

        if self.Gauss_matrix.shape != (self.Gauss_matrix.shape[0], 4):
            print("Gauss_matrix of configuration should have 4 colomns ")
            return False
        if not examiner(self.Gauss_matrix[:, 0:3]):
            print("Gauss_matrix of configuration is False")
            return False
        if not examiner(self.R3_type_global):
            print("Gauss_matrix of R3_type_global is False")
            return False
        label = True
        if (
            self.Gauss_matrix[0, 3] >= 0
            and self.Gauss_matrix[0, 3] <= 3
            and self.Gauss_matrix[self.Gauss_matrix.shape[0] - 1, 3] >= 0
            and self.Gauss_matrix[self.Gauss_matrix.shape[0] - 1, 3] <= 3
        ):
            for i in range(self.Gauss_matrix.shape[0] - 1):
                if self.Gauss_matrix[i + 1, 3] - self.Gauss_matrix[i, 3] < 0:
                    print("piece label of Gauss_matrix is False")
                    label = False
        return label

    def total_matrix(self):
        if self.Gauss_matrix.shape[0] == 0:
            return self.R3_type_global

        arr = self.Gauss_matrix[:, 3]
        temp_Gauss_matrix = self.Gauss_matrix[:, 0:3]
        #         print(temp_Gauss_matrix.shape)
        for i in range(temp_Gauss_matrix.shape[0]):
            temp_Gauss_matrix[i, 0] = 3 + temp_Gauss_matrix[i, 0]
        #         print(2*arr)
        result = np.insert(self.R3_type_global, 2 * arr, temp_Gauss_matrix, axis=0)
        #         print(result)
        return result

    def R3_index(self):
        """
        This function returns the R3_index to do R3 in total_matrix()
        """


@njit(cache=True)
def get_split(G, R3_index):
    G_split = np.zeros((G.shape[0] - 6, 4), dtype=np.int64)
    G_split[0 : R3_index[0], 0:3] = G[0 : R3_index[0], :]
    G_split[R3_index[0] : R3_index[1] - 2, 0:3] = G[R3_index[0] + 2 : R3_index[1], :]
    G_split[R3_index[1] - 2 : R3_index[2] - 4, 0:3] = G[R3_index[1] + 2 : R3_index[2], :]
    G_split[R3_index[2] - 4 : G.shape[0] - 6, 0:3] = G[R3_index[2] + 2 : G.shape[0], :]

    G_split[R3_index[0] : R3_index[1] - 2, 3] = 1
    G_split[R3_index[1] - 2 : R3_index[2] - 4, 3] = 2
    G_split[R3_index[2] - 4 : G.shape[0] - 6, 3] = 3
    return G_split


@njit(cache=True)
def cfg_unsigned_in_combination(G_split, cfg_Gauss_matrix_normed, ind_namelist_to_ind_matrix, c):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)

    (D, sign) = normalize(G_split[index, 0:4])
    D[:, 2] = 0
    if np.any(D - cfg_Gauss_matrix_normed) == 0:  # D==cfg_Gauss_matrix_normed
        return sign
    else:
        return 0


@njit(cache=True)
def cfg_signed_in_combination(G_split, cfg_Gauss_matrix_normed, ind_namelist_to_ind_matrix, c):
    c_array = np.array(c, dtype=np.int64)
    ind_temp = ind_namelist_to_ind_matrix[c_array, :]
    ind_temp = ind_temp.flatten()
    index = np.sort(ind_temp)

    (D, _) = normalize(G_split[index, 0:4])
    if np.any(D - cfg_Gauss_matrix_normed) == 0:  # D==cfg_Gauss_matrix_normed
        return 1
    else:
        return 0


def R3_configuration(G, R3_index, cfg):
    """
    G: Gauss diagram
    R3_index: [i,j,k] to perform R3 in G
    cfg: a Configuration class
    """
    if (
        np.any(
            normalize(
                G[
                    [
                        R3_index[0],
                        R3_index[0] + 1,
                        R3_index[1],
                        R3_index[1] + 1,
                        R3_index[2],
                        R3_index[2] + 1,
                    ],
                    :,
                ]
            )[0]
            - normalize(cfg.R3_type_global)[0]
        )
        != 0
    ):
        return 0
    #     print('R3 is the same')
    G_split = get_split(G, R3_index)

    ar = cfg.Gauss_matrix.shape[0] // 2
    cr = G_split.shape[0] // 2
    if cr < ar:
        return 0
    if ar == 0:
        return 1

    # R3_temp = R3(G,R3_index[0],R3_index[1],R3_index[2])# [G0, R3_typeindex, index_before_to_later, dhml]
    (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(G_split)
    cfg_Gauss_matrix_normed = normalize(cfg.Gauss_matrix)[0]

    result = 0
    if not cfg.signed:  # unsigned configuration
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = cfg_unsigned_in_combination(
                G_split, cfg_Gauss_matrix_normed, ind_namelist_to_ind_matrix, c
            )
            result += contr
    else:  # signed configuration
        for c in combinations(np.arange(name_list.shape[0]), ar):
            contr = cfg_signed_in_combination(
                G_split, cfg_Gauss_matrix_normed, ind_namelist_to_ind_matrix, c
            )
            result += contr

    return result


class Formula:
    def __init__(self, seq_cfg, seq_coeff):
        """
        Gauss_matrix: the matrix of the configuration
        R3_type_global is the global matrix of R3
        """
        self.length = len(seq_coeff)
        self.seq_cfg = seq_cfg
        self.seq_coeff = seq_coeff


def R3_formula(G, R3_index, fml):
    """
    G: Gauss diagram
    R3_index: [i,j,k] to perform R3 in G
    fml: a Formula class
    """
    result = 0
    for i in range(fml.length):
        result += R3_configuration(G, R3_index, fml.seq_cfg[i]) * fml.seq_coeff[i]

    return result


# ARNAUD_MORTIER
seq_cfg = []
seq_coeff = []
matrix_0 = np.array([[1, 1, 0, 0], [1, -1, 0, 1]], int)
matrix_1 = np.array([[1, -1, 0, 0], [1, 1, 0, 1]], int)
matrix_2 = np.array([[1, 1, 0, 1], [1, -1, 0, 3]], int)

for i in range(16, 32):
    seq_cfg += [Configuration(matrix_0, infty_avancer(R3_table_matrix[:, :, i], 2))]
    seq_coeff += [(-1) ** i]
# for i in range(len(seq_cfg)):
#     draw(seq_cfg[i].total_matrix(),[i,seq_coeff[i]],[0,3])

for i in range(16):
    seq_cfg += [Configuration(matrix_1, infty_avancer(R3_table_matrix[:, :, i], 4))]
    seq_coeff += [(-1) ** i]

for i in range(16):
    seq_cfg += [Configuration(matrix_2, infty_avancer(R3_table_matrix[:, :, i], 4))]
    seq_coeff += [(-1) ** i]

seq_coeff = np.array(seq_coeff)
Arnaud_Mortier = Formula(seq_cfg, seq_coeff)


def generate_cfg_matrix_unsigned(numb_arrows):
    """ """

    (chord, directed, directed_signed, numb_chord, count_directed, count_directed_signed) = (
        diagram_generate_param(numb_arrows)
    )
    x = np.zeros(4, int)  # x0+x1+x2+x3 = 2*numb_arrow
    list_cfg_matrix = []
    for i in range(count_directed):
        for c in combinations(np.arange(1, 2 * numb_arrows + 4), 3):
            x[0] = c[0] - 1
            x[1] = c[1] - c[0] - 1
            x[2] = c[2] - c[1] - 1
            x[3] = 2 * numb_arrows + 3 - c[2]
            temp_matrix = np.zeros([directed.shape[0], 4], int)
            temp_matrix[:, 0:3] = copy.deepcopy(directed[:, :, i])
            temp_matrix[0 : x[0], 3] = 0
            temp_matrix[x[0] : x[0] + x[1], 3] = 1
            temp_matrix[x[0] + x[1] : x[0] + x[1] + x[2], 3] = 2
            temp_matrix[x[0] + x[1] + x[2] : x[0] + x[1] + x[2] + x[3], 3] = 3
            list_cfg_matrix += [temp_matrix]
    return np.array(list_cfg_matrix, dtype=np.int64)


def generate_list_right_left_formula(numb_arrows):
    """
    param: numb_arrows:number of arrows outside R3
    ouput: list_formula
    """
    list_formula = []
    if numb_arrows == 0:
        for r in range(3):
            for t in range(2):
                seq_cfg = []
                seq_coeff = np.array([])
                for i in range(t * 16, (t + 1) * 16):  # left/right type
                    seq_cfg += [
                        Configuration(
                            np.zeros([0, 4], int), infty_avancer(R3_table_matrix[:, :, i], 2 * r)
                        )
                    ]
                    seq_coeff = np.append(seq_coeff, [(-1) ** i])

                list_formula += [Formula(seq_cfg, seq_coeff)]
        return list_formula

    list_cfg_matrix = generate_cfg_matrix_unsigned(numb_arrows)
    for s in range(len(list_cfg_matrix)):
        for r in range(3):
            for t in range(2):
                # t ==0 left type
                # t ==1 right type
                cfg_matrix = list_cfg_matrix[s]
                seq_cfg = []
                seq_coeff = np.array([])
                for i in range(t * 16, (t + 1) * 16):  # left/right type
                    seq_cfg += [
                        Configuration(cfg_matrix, infty_avancer(R3_table_matrix[:, :, i], 2 * r))
                    ]
                    seq_coeff = np.append(seq_coeff, [(-1) ** i])
                list_formula += [Formula(seq_cfg, seq_coeff)]
    return list_formula


@njit(cache=True)
def R3_identify_r_l(G, R3_index, R3_matrix_cfg):
    sum_R3_piece_1 = R3_matrix_cfg[0, 1] + R3_matrix_cfg[1, 1]
    sum_R3_piece_2 = R3_matrix_cfg[2, 1] + R3_matrix_cfg[3, 1]
    sum_R3_piece_3 = R3_matrix_cfg[4, 1] + R3_matrix_cfg[5, 1]

    sum_G_piece_1 = G[R3_index[0], 1] + G[R3_index[0] + 1, 1]
    sum_G_piece_2 = G[R3_index[1], 1] + G[R3_index[1] + 1, 1]
    sum_G_piece_3 = G[R3_index[2], 1] + G[R3_index[2] + 1, 1]

    if (
        sum_R3_piece_1 != sum_G_piece_1
        or sum_R3_piece_2 != sum_G_piece_2
        or sum_R3_piece_3 != sum_G_piece_3
    ):
        return False
    else:
        return True


@njit(cache=True)
def get_R3_sign(R3_matrix):
    """
    R3_matrix[6,3]
    output: the sign of the R3
    """
    if R3_matrix[2, 0] == R3_matrix[0, 0]:
        if R3_matrix[4, 0] == R3_matrix[1, 0]:
            sign = 1
        else:  # R3_matrix[5,0] ==  R3_matrix[1,0]:
            sign = -1
    elif R3_matrix[3, 0] == R3_matrix[0, 0]:
        if R3_matrix[5, 0] == R3_matrix[1, 0]:
            sign = 1
        else:  # R3_matrix[4,0] ==  R3_matrix[1,0]:
            sign = -1
    elif R3_matrix[4, 0] == R3_matrix[0, 0]:
        if R3_matrix[3, 0] == R3_matrix[1, 0]:
            sign = 1
        else:  # R3_matrix[2,0] ==  R3_matrix[1,0]:
            sign = -1
    else:  # R3_matrix[5,0] ==  R3_matrix[0,0]:
        if R3_matrix[2, 0] == R3_matrix[1, 0]:
            sign = 1
        else:  # R3_matrix[3,0] ==  R3_matrix[1,0]:
            sign = -1
    return sign


@njit(cache=True)
def r_l_R3_type(R3_matrix):
    """
    R3_matrix: [6,3]
    output: the index of the order : [l,h,m],[h,m,l],[m,l,h], [l,m,h],[m,h,l],[h,l,m].
    Warning: This function does not check whether R3_matrix is right for R3.
    """
    a = R3_matrix[0, 1] + R3_matrix[1, 1]
    b = R3_matrix[2, 1] + R3_matrix[3, 1]
    if a == -2:
        if b == 2:
            index = 0
        else:  # b == 0:
            index = 3
    elif a == 0:
        if b == -2:
            index = 2
        else:  # b == 2:
            index = 4
    else:  # a == 2
        if b == 0:
            index = 1
        else:  # b == -2:
            index = 5
    return index


@njit(cache=True)
def get_R3_pairs_type(R3_matrix):
    a = R3_matrix[0, 1] + R3_matrix[1, 1]
    b = R3_matrix[2, 1] + R3_matrix[3, 1]
    alter_numb_h = 0  # alter number of h
    if a == -2:
        if b == 2:
            global_type = 0  # [l,h,m]
            if R3_matrix[1, 0] == R3_matrix[2, 0] or R3_matrix[5, 0] == R3_matrix[2, 0]:
                alter_numb_h += 1
            if R3_matrix[0, 0] == R3_matrix[3, 0] or R3_matrix[4, 0] == R3_matrix[3, 0]:
                alter_numb_h += 1
        else:  # b == 0:
            global_type = 3  # [l,m,h]
            if R3_matrix[1, 0] == R3_matrix[4, 0] or R3_matrix[3, 0] == R3_matrix[4, 0]:
                alter_numb_h += 1
            if R3_matrix[0, 0] == R3_matrix[5, 0] or R3_matrix[2, 0] == R3_matrix[5, 0]:
                alter_numb_h += 1
    elif a == 0:
        if b == -2:
            global_type = 2  # [m,l,h]
            if R3_matrix[1, 0] == R3_matrix[4, 0] or R3_matrix[3, 0] == R3_matrix[4, 0]:
                alter_numb_h += 1
            if R3_matrix[0, 0] == R3_matrix[5, 0] or R3_matrix[2, 0] == R3_matrix[5, 0]:
                alter_numb_h += 1
        else:  # b == 2:
            global_type = 4  # [m,h,l]
            if R3_matrix[1, 0] == R3_matrix[2, 0] or R3_matrix[5, 0] == R3_matrix[2, 0]:
                alter_numb_h += 1
            if R3_matrix[0, 0] == R3_matrix[3, 0] or R3_matrix[4, 0] == R3_matrix[3, 0]:
                alter_numb_h += 1
    else:  # a == 2
        if b == 0:
            global_type = 1  # [h,m,l]
            if R3_matrix[3, 0] == R3_matrix[0, 0] or R3_matrix[5, 0] == R3_matrix[0, 0]:
                alter_numb_h += 1
            if R3_matrix[2, 0] == R3_matrix[1, 0] or R3_matrix[4, 0] == R3_matrix[1, 0]:
                alter_numb_h += 1
        else:  # b == -2:
            global_type = 5  # [h,l,m]
            if R3_matrix[3, 0] == R3_matrix[0, 0] or R3_matrix[5, 0] == R3_matrix[0, 0]:
                alter_numb_h += 1
            if R3_matrix[2, 0] == R3_matrix[1, 0] or R3_matrix[4, 0] == R3_matrix[1, 0]:
                alter_numb_h += 1
    sum_sign = np.sum(R3_matrix[:, 2])
    if alter_numb_h == 0:
        if sum_sign == 2:
            local_type = 0
        else:  # sum_sign == -2:
            local_type = 1
    elif alter_numb_h == 1:
        if sum_sign == 2:
            local_type = 2
        elif sum_sign == -2:
            local_type = 3
        elif sum_sign == 6:
            local_type = 4
        else:  # sum_sign == -6:
            local_type = 5
    else:  # alter_numb_h == 2:
        if sum_sign == 2:
            local_type = 6
        else:  # sum_sign == -2:
            local_type = 7
    R3_pairs_type = global_type * 8 + local_type
    return R3_pairs_type


# =============================================================================
# def R3_right_left_formula(G,R3_index,r_l_fml):
#     """
#     G: Gauss diagram
#     R3_index: [i,j,k] to perform R3 in G
#     fml: a right or left type Formula class generated by
#             generate_list_right_left_formula()
#     """
#     cfg = r_l_fml.seq_cfg[0] # select one cfg
#     R3_matrix_cfg = cfg.R3_type_global
#
#     be_or_not = R3_identify_r_l(G,R3_index,R3_matrix_cfg)
#     if be_or_not == False:
#         return 0
#
#     ####  R3 has matched
#     #### now find the R3_sign
#     R3_sign = get_R3_sign(G,R3_index)
#
#     ####R3 part is ok
#     ####now do the match part for the rest
#
#     G_split = get_split(G,R3_index)
#     ar = cfg.Gauss_matrix.shape[0] // 2
#     cr = G_split.shape[0] // 2
#     if cr < ar :
#         return 0
#     if ar == 0:
#         return 1
#
#     (name_list, ind_namelist_to_ind_matrix, ind_matrix_to_ind_namelist) = index_names_trans(G_split)
#     cfg_Gauss_matrix_normed = normalize(cfg.Gauss_matrix)[0]
#
#     result = 0
#     if cfg.signed == False: #unsigned configuration
#         for c in combinations(np.arange(name_list.shape[0]),ar):
#             contr = cfg_unsigned_in_combination(G_split,cfg_Gauss_matrix_normed,ind_namelist_to_ind_matrix,c)
#             result += contr
#     else:#signed configuration
#         for c in combinations(np.arange(name_list.shape[0]),ar):
#             contr = cfg_signed_in_combination(G_split,cfg_Gauss_matrix_normed,ind_namelist_to_ind_matrix,c)
#             result += contr
#
#     return result*R3_sign
# =============================================================================
