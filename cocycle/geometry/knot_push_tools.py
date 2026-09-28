"""Knot pushing helpers recovered from the original local project."""

import numpy as np
from numba import njit
from .basic_tools_jit import getname, left_right, oppend


@njit(cache=True)
def R2_forpush(G, i):
    """
    G: Guass diagram matrix of a long knot
    i: the first index of the consecutive two to be pushed
    Outputs:[G_result, index_tracked in G_result]
    """
    # G0 = G.copy()
    row0 = G.shape[0]
    ind_tr = np.arange(0, row0)

    names = getname(G)
    name_max = np.max(names)
    a = oppend(G, i)
    b = oppend(G, i + 1)

    codirection = G[i, 1] * G[i + 1, 1] * G[i, 2] * G[i + 1, 2]
    if (
        codirection == 1
    ):  # codirection == 1 if the two strands have the same direction; otherwise -1.
        if left_right(G, a)[1] == i + 1:
            if G[i + 1, 1] == 1:  # Going BELOW !!!
                # print('同向，在右边，下面走')
                arr = np.array(
                    [
                        [name_max + 1, 1, 1],
                        [name_max + 2, 1, -1],
                        [name_max + 1, -1, 1],
                        [name_max + 2, -1, -1],
                    ],
                    dtype=np.int64,
                )

            else:  # going above
                # print('同向，在右边，上面走')
                arr = np.array(
                    [
                        [name_max + 1, -1, -1],
                        [name_max + 2, -1, 1],
                        [name_max + 1, 1, -1],
                        [name_max + 2, 1, 1],
                    ],
                    dtype=np.int64,
                )

        else:  # left_right(G0,a)[0] == i+1
            if G[i + 1, 1] == 1:  # Going BELOW !!!
                # print('同向，在左边，下面走')
                arr = np.array(
                    [
                        [name_max + 1, 1, -1],
                        [name_max + 2, 1, 1],
                        [name_max + 1, -1, -1],
                        [name_max + 2, -1, 1],
                    ],
                    dtype=np.int64,
                )

            else:  # going above
                # print('同向，在左边，上面走')
                arr = np.array(
                    [
                        [name_max + 1, -1, 1],
                        [name_max + 2, -1, -1],
                        [name_max + 1, 1, 1],
                        [name_max + 2, 1, -1],
                    ],
                    dtype=np.int64,
                )
        #         G0 = np.insert(G0,[a,a,b,b],arr,axis = 0)
        #         print('case_1')
        if a < b:
            # for j in range(row0):
            #     if j>= a and j < b:
            #         ind_tr[j] += 2
            #     elif j >= b:
            #         ind_tr[j] += 4
            ind_tr[a:b] += 2
            ind_tr[b:row0] += 4
            #             G0 = np.insert(G0,[a,a,b,b],arr,axis = 0)
            G_new = np.zeros((G.shape[0] + 4, 3), dtype=np.int64)
            G_new[0:a, :] = G[0:a, :]
            G_new[a : a + 2, :] = arr[0:2, :]
            G_new[a + 2 : b + 2, :] = G[a:b, :]
            G_new[b + 2 : b + 4, :] = arr[2:4, :]
            G_new[b + 4 : G_new.shape[0], :] = G[b : G.shape[0], :]

        else:  # a > b
            # for j in range(row0):
            #     if j>= b and j < a:
            #         ind_tr[j] += 2
            #     elif j >= a:
            #         ind_tr[j] += 4
            ind_tr[b:a] += 2
            ind_tr[a:row0] += 4
            #             G0 = np.insert(G0,[a,a,b,b],arr,axis = 0)
            G_new = np.zeros((G.shape[0] + 4, 3), dtype=np.int64)
            G_new[0:b, :] = G[0:b, :]
            G_new[b : b + 2, :] = arr[2:4, :]
            G_new[b + 2 : a + 2, :] = G[b:a, :]
            G_new[a + 2 : a + 4, :] = arr[0:2, :]
            G_new[a + 4 : G_new.shape[0], :] = G[a : G.shape[0], :]

    else:  # codirection == -1
        if left_right(G, a)[1] == i + 1:
            if G[i + 1, 1] == 1:  # Going BELOW !!!
                # print('异向，在右边，下面走')
                arr = np.array(
                    [
                        [name_max + 1, 1, -1],
                        [name_max + 2, 1, 1],
                        [name_max + 2, -1, 1],
                        [name_max + 1, -1, -1],
                    ],
                    dtype=np.int64,
                )
            else:  # Going above
                # print('异向，在右边，上面走')
                arr = np.array(
                    [
                        [name_max + 1, -1, 1],
                        [name_max + 2, -1, -1],
                        [name_max + 2, 1, -1],
                        [name_max + 1, 1, 1],
                    ],
                    dtype=np.int64,
                )
        else:  # left_right(G0,a)[0] == i+1
            if G[i + 1, 1] == 1:  # Going BELOW !!!
                #                 print('异向，在左边，下面走')
                arr = np.array(
                    [
                        [name_max + 1, 1, 1],
                        [name_max + 2, 1, -1],
                        [name_max + 2, -1, -1],
                        [name_max + 1, -1, 1],
                    ],
                    dtype=np.int64,
                )
            else:  # Going above
                # print('异向，在左边，上面走')
                arr = np.array(
                    [
                        [name_max + 1, -1, -1],
                        [name_max + 2, -1, 1],
                        [name_max + 2, 1, 1],
                        [name_max + 1, 1, -1],
                    ],
                    dtype=np.int64,
                )
        if b != row0 - 1:
            if a != b + 1:
                #                 G0 = np.insert(G0,[a,a,b+1,b+1],arr,axis = 0)
                #                 #draw(G0,"a={},b+1={}".format(a,b+1))
                #                 print('case_2')
                G_new = np.zeros((G.shape[0] + 4, 3), dtype=np.int64)
                if a <= b + 1:
                    G_new[0:a, :] = G[0:a, :]
                    G_new[a : a + 2, :] = arr[0:2, :]
                    G_new[a + 2 : b + 1 + 2, :] = G[a : b + 1, :]
                    G_new[b + 1 + 2 : b + 1 + 4, :] = arr[2:4, :]
                    G_new[b + 1 + 4 : G_new.shape[0], :] = G[b + 1 : G.shape[0], :]
                else:  # a > b+1
                    G_new[0 : b + 1, :] = G[0 : b + 1, :]
                    G_new[b + 1 : b + 1 + 2, :] = arr[2:4, :]
                    G_new[b + 1 + 2 : a + 2, :] = G[b + 1 : a, :]
                    G_new[a + 2 : a + 4, :] = arr[0:2, :]
                    G_new[a + 4 : G_new.shape[0], :] = G[a : G.shape[0], :]

            else:
                arr_1 = arr.copy()
                arr_1[0, :] = arr[2, :]
                arr_1[1, :] = arr[3, :]
                arr_1[2, :] = arr[0, :]
                arr_1[3, :] = arr[1, :]
                #                 G0 = np.insert(G0,[a,a,b+1,b+1],arr_1,axis = 0)
                #                 #draw(G0,"a={},b+1={}".format(a,b+1))
                #                 print('case_3')
                G_new = np.zeros((G.shape[0] + 4, 3), dtype=np.int64)
                if a <= b + 1:
                    G_new[0:a, :] = G[0:a, :]
                    G_new[a : a + 2, :] = arr_1[0:2, :]
                    G_new[a + 2 : b + 1 + 2, :] = G[a : b + 1, :]
                    G_new[b + 1 + 2 : b + 1 + 4, :] = arr_1[2:4, :]
                    G_new[b + 1 + 4 : G_new.shape[0], :] = G[b + 1 : G.shape[0], :]
                else:  # a > b+1:
                    G_new[0 : b + 1, :] = G[0 : b + 1, :]
                    G_new[b + 1 : b + 1 + 2, :] = arr_1[2:4, :]
                    G_new[b + 1 + 2 : a + 2, :] = G[b + 1 : a, :]
                    G_new[a + 2 : a + 4, :] = arr_1[0:2, :]
                    G_new[a + 4 : G_new.shape[0], :] = G[a : G.shape[0], :]
        #                 else: # a == b+1
        #                     G_new[0 : a,:] = G[0 : a,:]
        #                     G_new[a : a+4,:] = arr_1
        #                     G_new[a+4 : G_new.shape[0],:] = G[a : G.shape[0], :]

        else:
            #             G0 = np.insert(G0,[a,a],arr[0:2,:],axis = 0)
            #             G0 = np.append(G0,arr[2:4,:],axis = 0)
            #             print('case_4')
            G_new = np.zeros((G.shape[0] + 4, 3), dtype=np.int64)
            G_new[0:a, :] = G[0:a, :]
            G_new[a : a + 2, :] = arr[0:2, :]
            G_new[a + 2 : G.shape[0] + 2, :] = G[a : G.shape[0], :]
            G_new[G.shape[0] + 2 : G.shape[0] + 4, :] = arr[2:4, :]

        if a < b:
            # for j in range(row0):
            #     if j>= a and j <= b:
            #         ind_tr[j] += 2
            #     elif j > b:
            #         ind_tr[j] += 4
            ind_tr[a : b + 1] += 2
            ind_tr[b + 1 : row0] += 4

        else:  # a > b
            # for j in range(row0):
            #     if j > b and j < a:
            #         ind_tr[j] += 2
            #     elif j >= a:
            #         ind_tr[j] += 4
            ind_tr[b + 1 : a] += 2
            ind_tr[a:row0] += 4
    return (G_new, ind_tr)


@njit(cache=True)
def R2_push_pre(G, ind, direction=-1):
    """
    G:Gauss matrix
    ind:the index consecutive to be pushed
    direction=-1:push ind to ind-1,
    direction =1: push ind to ind+1
    """
    if direction == -1:
        if ind > 0 and ind < G.shape[0]:
            #             print('direc=-1')
            #             print('ind-1=',ind-1)
            (G_R2, index_before_to_later) = R2_forpush(G, ind - 1)
        else:
            raise ValueError("Index cannot do push backward.")
    elif direction == 1:
        #         print('direc=1')
        #         print('ind=',ind)
        if ind >= 0 and ind < G.shape[0] - 1:
            (G_R2, index_before_to_later) = R2_forpush(G, ind)
        else:
            raise ValueError("Index cannot do push forward.")
    else:
        raise ValueError("direction must be -1 or 1")
    return (G_R2, index_before_to_later)


@njit(cache=True)
def circledetect(G, starts, ends):
    """
    :param G: Gauss diagram having deleted last circle
    :param starts: an 1-d index array for G giving every beginning point
    :param ends: an index array for G giving all possible ends
    Outputs:
    [matrix with each row a circle, the index array giving the meaningful length of each row, new_starts, new_ends]
    """
    # Notice that starts in the same TURN will not get collide with each other on their orbits.

    row = G.shape[0]
    ns = starts.shape[0]  # number of starting points
    # ends = -np.ones(ns,dtype = np.int64)
    # for i in range(ns):
    # ends[i] = (starts[i]+1) % row # + adapt to push in anti-direction
    # print('ends=',ends)
    C = -np.ones((ns, row), dtype=np.int64)
    nc = np.zeros(ns, dtype=np.int64)  # the length of meaningful part of each row in C

    for i in range(ns):
        step = starts[i]
        # print('step=', step)
        s = 0
        C[i, s] = step
        step = oppend(G, step)
        # print('step=', step)
        s = s + 1
        C[i, s] = step

        while (((step - 1) % row) in ends) == 0:  # when step is not in starts
            step = (step - 1) % row  # -1 means we push in the anti-direction
            # print('step=', step)
            s = s + 1
            C[i, s] = step
            step = oppend(G, step)
            # print('step=', step)
            s = s + 1
            C[i, s] = step
            # print('next while condition: belong(%d )' % ((step-1) % row), ends)
    ####
    nc_total = 0
    for i in range(ns):
        num_c = 0
        for j in range(row):
            if C[i, j] >= 0:
                num_c = num_c + 1
                nc_total += 1
        nc[i] = num_c

    together = -np.ones(
        nc_total, dtype=np.int64
    )  # making all rows in C into one row(just meaningful part)
    count = 0
    for i in range(ns):
        for j in range(nc[i]):
            together[count] = C[i, j]
            count += 1
    # print('together')
    # print(together)

    new_starts = -np.ones(row, dtype=np.int64)
    count = 0
    for i in range(nc_total // 2):
        if ((together[2 * i] - 1) in together) == 0:  # odd terms -1 not in together
            new_starts[count] = together[2 * i] - 1
            count = count + 1
    ind_cut = 0
    for i in range(row):
        if new_starts[i] >= 0:
            ind_cut = ind_cut + 1
    new_starts = new_starts[0:ind_cut]

    new_ends = -np.ones(nc_total // 2, dtype=np.int64)
    for i in range(nc_total // 2):
        new_ends[i] = together[2 * i + 1]

    return (C, nc, new_starts, new_ends)


@njit(cache=True)
def seifertcircle(G):
    """
    :param G: Gauss diagram
    Outputs:[Circle matrix[each starts,:,each turn], number of turns in total,
    decipher[each turn,meaning for parts for C[:,:,turn]] inherite from circledetect,
    dofd: 1d array. On dofd[i] is the meaningful lenghth of decipher[i,:]
    """
    row_G = G.shape[0]
    C = -np.ones((row_G, row_G, row_G // 2), dtype=np.int64)
    decipher = -np.ones((row_G // 2, row_G), dtype=np.int64)

    turn = 0
    # print(np.array([row_G],dtype = np.int64))
    temp = circledetect(
        G, np.array([row_G - 1], dtype=np.int64), np.array([row_G - 1], dtype=np.int64)
    )
    # print(temp[0].shape[0],temp[0].shape[1])
    # print(temp[1])
    C[0 : temp[0].shape[0], 0 : temp[0].shape[1], turn] = temp[0]
    decipher[turn, 0 : temp[1].shape[0]] = temp[1]

    while temp[2].shape[0] != 0:
        turn = turn + 1
        temp = circledetect(G, temp[2], temp[3])
        # print('turn = ', turn)
        # print(temp[1])
        C[0 : temp[0].shape[0], 0 : temp[0].shape[1], turn] = temp[0]
        decipher[turn, 0 : temp[1].shape[0]] = temp[1]

    turn = turn + 1

    dofd = -np.ones(turn, dtype=np.int64)
    for i in range(turn):
        for j in range(row_G):
            if decipher[i, j] < 0:
                dofd[i] = j
                break
    return (C, turn, decipher, dofd)
