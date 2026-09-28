from .geometry.basic_tools_jit import *
from .geometry.knot_push_tools import seifertcircle, R2_push_pre
from . import matching as btc
@njit
def R3_deform(G, R3_index):
    G_new = G.copy()

    i = R3_index[0]
    j = R3_index[1]
    k = R3_index[2]

    G_new[i] = G[i+1]
    G_new[i+1] = G[i]
    G_new[j] = G[j+1]
    G_new[j+1] = G[j]
    G_new[k] = G[k+1]
    G_new[k+1] = G[k]

    return G_new

# @njit   
def push_rl_cfgs(G, ind, arr_cfgs, arr_R3_pairs_types, direction = -1):
    """
    G: Gauss matrix of a long knot
    nG: a class cabled_knot
    ind:the index to be pushed in G_single(on the path,so here is begining of the consecutive 2)
    direction=-1: push backward
    # direction = 1: push forward
    """
    # draw(G,'single G before push in push_rl_cfgs', np.array([ind], int))
    (G_R2,index_before_R2_to_later) = R2_push_pre(G,ind,direction)

    a = oppend(G_R2,index_before_R2_to_later[ind-1])
    b = oppend(G_R2,index_before_R2_to_later[ind])
    if G[ind-1,1]*G[ind,1]*G[ind-1,2]*G[ind,2] == 1:# codirection ==1, two strands have the same direction
        arr_ind = np.sort(np.array([a-1,index_before_R2_to_later[ind-1],b-1],dtype = np.int64))
    else:
        arr_ind = np.sort(np.array([a-1,index_before_R2_to_later[ind-1],b],dtype = np.int64))
    (
        G_new, 
        index_before_R3_to_later, 
        contribution
        ) = btc.R3_arr_rl_cfgs(G_R2, arr_ind,  arr_cfgs, arr_R3_pairs_types)
    # draw(G_new, 'G_new, contribution = {}'.format(contribution), np.append(arr_ind, arr_ind + 1))

    index_before_to_later = index_before_R3_to_later[index_before_R2_to_later]
    return (G_new, index_before_to_later, contribution )

# @njit
def push_along_rl_cfgs(G, ind, path, arr_cfgs, arr_R3_pairs_types, direction = -1):
    """
    G: Gauss matrix
    ind: the index to be pushed along the path a
    path: path[0] need to ind-1(if direction == -1)/ind+1(if direction == 1)
    path need have even number of terms i.e., [ind0,oppend(ind0),...] the path from serfertcicle...
    arr_index_cfg: an array of index in [0,2520)
    """
    contribution = np.zeros(len(arr_R3_pairs_types), dtype = np.int64)
    G_new = G.copy()
    if path[0] != ind+direction:
        raise TypeError("ind not pushable along path")
    index_before_to_later = np.arange(G.shape[0])
    path_evolve = path.copy()
    len_path = len(path) // 2
    for i in range(len_path):        
        (
            G_new, 
            index_before_push_to_later, 
            contri
            ) = push_rl_cfgs(G_new, path_evolve[2*i]-direction, arr_cfgs, arr_R3_pairs_types, direction = -1)
        contribution += contri
        index_before_to_later = index_before_push_to_later[index_before_to_later]
        path_evolve = index_before_push_to_later[path_evolve] 
    
    return (G_new, index_before_to_later, contribution )

# @njit
def push_through_rl_cfgs(K, G, arr_cfgs, arr_R3_pairs_types):
    """
    K:the 1st long knot 
    G:the 2nd long knot
    arr_index_cfg: an array of index in [0,2520)
    we view K as a small point pushing it from left through the big G
    Output: contribution 
    """
    K1 = K
    G1 = G
    S1 = sum_knot(K1,G1)
    S = S1.copy()
    
    row_K = K1.shape[0]
    row_G = G1.shape[0]

    circles = seifertcircle(K1)
    turn_total = circles[1]
    length_decipher = circles[3]
    decipher = circles[2]
    ind_tr = np.arange(0,row_K)
#     print('turn_total',turn_total)
#     print('length_decipher',length_decipher)
#     print('decipher',decipher)
    len_arr_index_cfg = len(arr_R3_pairs_types)
    contribution = np.zeros(len_arr_index_cfg ,dtype = np.int64)

    for i in range(row_G):
        # print('i=',i)
        for j in range(turn_total):
            # print('j=',j)
            matrix_path = circles[0][0:length_decipher[j],:,j]
            #print('matrix_path',matrix_path)
            newmatrix = matrix_path.copy()
            for k in range(length_decipher[j]):
                # print('k=',k)
                length_path = decipher[j,k]
                #print('length_path=',length_path)
                for a in range(length_decipher[j]):
                    newmatrix[a, 0 : decipher[j,a]] = ind_tr[matrix_path[a, 0 : decipher[j, a]]]
                # (G_new, index_before_to_later, contribution ) = push_along_cfgs(G, ind, path, arr_index_cfg, direction = -1)
                (S, index_before_to_later, contri ) = push_along_rl_cfgs(
                    S,
                    newmatrix[k,0:length_path][0]+1,
                    newmatrix[k,0:length_path],
                    arr_cfgs, 
                    arr_R3_pairs_types,
                    direction =-1
                )
                ind_tr = index_before_to_later[ind_tr]
                #########Calculation part
                
                contribution += contri
                #########
#         print('contribution=',contribution)
        # the following part is faster than using njit. 
        if i < row_G-1:   
            S_temp = S1.copy()
            S_temp[0:(i+1),:] = S1[row_K:(row_K+i+1),:]
            S_temp[(i+1):(i+1+row_K),:] = S1[0:row_K,:]
            S_temp[(i+1+row_K):(row_K+row_G),:] = S1[(row_K+i+1):(row_K+row_G),:]
            ind_tr = np.arange(0,row_K)+i+1
            S = S_temp

    return contribution

# @njit
def exchange_loop_rl_cfgs(K, G, arr_cfgs, arr_R3_pairs_types):
    """
    arr_index_cfg: an array of index in [0,2520)
    Output: contribution
    """
    contribution_0 = push_through_rl_cfgs(K, G, arr_cfgs, arr_R3_pairs_types)
    contribution_1 = push_through_rl_cfgs(G, K, arr_cfgs, arr_R3_pairs_types)
    contribution = contribution_0 + contribution_1

    return contribution

@njit
def curl_prepare(G):
    """
    param: G: Gauss diagram matrix
    param: ind_tr: index in G to be tracked
    output:[G_result, index_before_to_after]
    G_result is the matrix with 1 curls added for FH loop.

    """
    row_G = G.shape[0]
    name_max = np.max(G[:,0])
    index_before_to_after = np.arange(G.shape[0])
    # ind = ind_namelist_to_ind_matrix[ind_matrix_to_ind_namelist[row_G-1],0]
    ind = oppend(G,row_G-1)

    if G[ind,1] == 1:
        arr = np.array([
            [name_max+1,-1,G[ind,2]],
            [name_max+1,1,G[ind,2]],
        ],dtype = np.int64)
    elif G[ind,1] == -1:
        arr = np.array([
            [name_max+1,1,G[ind,2]],
            [name_max+1,-1,G[ind,2]],
        ],dtype = np.int64)

    G0 = np.zeros((row_G+2,G.shape[1]),dtype = np.int64)
    G0[0 : ind] = G[0 : ind]
    G0[ind : ind + 2] = arr
    G0[ind + 2 : G0.shape[0]] = G[ind : row_G]

    index_before_to_after[ind : row_G] += 2
    return (G0, index_before_to_after)


# @njit
def Fox_Hatcher_rl_cfgs(K, arr_cfgs, arr_R3_pairs_types):
    """
    K:Gauss matrix 
    arr_index_cfg: an array of index in [0,2520)
    Output:contribution
    """
    row_K = K.shape[0]
    K_temp = K.copy()
    contribution = np.zeros(len(arr_R3_pairs_types), dtype = np.int64)
    for i in range(row_K):
        # print('i=',i)
        circles = seifertcircle(K_temp)
        turn_total = circles[1]
        length_decipher = circles[3]
        decipher = circles[2]
        (K_evolve, ind_before_curl_to_after) = curl_prepare(K_temp)
        ind_tr = ind_before_curl_to_after
        # print('turn_total = ', turn_total)
        for j in range(turn_total):
            # print('j=',j)
            matrix_path = circles[0][0:length_decipher[j],:,j]
            newmatrix = matrix_path.copy()
#             print('length_decipher[j]=',length_decipher[j])
            for k in range(length_decipher[j]):
                # print('k=',k)
                length_path = decipher[j,k]
#                 print('length_path',length_path)
                for a in range(length_decipher[j]):
                    newmatrix[a, 0 : decipher[j,a]] = ind_tr[matrix_path[a, 0 : decipher[j, a]]]
                if j == 0:#fisrt move after having added a curl 
                    if length_path != 2:
                        (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,2:length_path][0]+1, 
                            newmatrix[k,2:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                        
                        ind_tr = index_before_to_later[ind_tr]
                        contribution += contri
                        # print('contri = ', contri)
                else:# not the first move
                    (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,0:length_path][0]+1, 
                            newmatrix[k,0:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                    ind_tr = index_before_to_later[ind_tr]
                    contribution += contri
                    # print('contri = ', contri)
        # draw(K_evolve, 'i = {}'.format(i), ind_tr)
        if i < row_K-1:
            K_temp = infty_avancer(K_temp,-1)
    
    return contribution

def half_Fox_Hatcher_rl_cfgs(K, arr_cfgs, arr_R3_pairs_types):
    """
    K:Gauss matrix 
    arr_index_cfg: an array of index in [0,2520)
    Output:contribution
    """
    row_K = K.shape[0]
    K_temp = K.copy()
    contribution = np.zeros(len(arr_R3_pairs_types), dtype = np.int64)
    # draw(K_temp,'begin')
    for i in range(row_K // 2):
        # print('i=',i)
        circles = seifertcircle(K_temp)
        turn_total = circles[1]
        length_decipher = circles[3]
        decipher = circles[2]
        (K_evolve, ind_before_curl_to_after) = curl_prepare(K_temp)
        ind_tr = ind_before_curl_to_after
        # print('turn_total = ', turn_total)
        for j in range(turn_total):
            # print('j=',j)
            matrix_path = circles[0][0:length_decipher[j],:,j]
            newmatrix = matrix_path.copy()
#             print('length_decipher[j]=',length_decipher[j])
            for k in range(length_decipher[j]):
                # print('k=',k)
                length_path = decipher[j,k]
#                 print('length_path',length_path)
                for a in range(length_decipher[j]):
                    newmatrix[a, 0 : decipher[j,a]] = ind_tr[matrix_path[a, 0 : decipher[j, a]]]
                if j == 0:#fisrt move after having added a curl 
                    if length_path != 2:
                        (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,2:length_path][0]+1, 
                            newmatrix[k,2:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                        
                        ind_tr = index_before_to_later[ind_tr]
                        contribution += contri
                        # print('contri = ', contri)
                else:# not the first move
                    (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,0:length_path][0]+1, 
                            newmatrix[k,0:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                    ind_tr = index_before_to_later[ind_tr]
                    contribution += contri
                    # print('contri = ', contri)
        # draw(K_evolve, 'i = {}'.format(i), ind_tr)
        if i < row_K-1:
            K_temp = infty_avancer(K_temp,-1)
    # draw(K_temp,'half')
    return contribution

def root_rolling_rl_cfgs(K, n, arr_cfgs, arr_R3_pairs_types):
    """
    K:Gauss matrix 
    n: a positive integer, to do n-th root of rolling
    arr_cfgs: an array of cfgs
    arr_R3_pairs_types: arr_R3_pairs_types

    Output:arr of contribution
    """
    row_K = K.shape[0]
    K_temp = K.copy()
    contribution = np.zeros(len(arr_R3_pairs_types), dtype = np.int64)
    # draw(K_temp,'begin')
    for i in range(row_K // n):
        # print('i=',i)
        circles = seifertcircle(K_temp)
        turn_total = circles[1]
        length_decipher = circles[3]
        decipher = circles[2]
        (K_evolve, ind_before_curl_to_after) = curl_prepare(K_temp)
        ind_tr = ind_before_curl_to_after
        # print('turn_total = ', turn_total)
        for j in range(turn_total):
            # print('j=',j)
            matrix_path = circles[0][0:length_decipher[j],:,j]
            newmatrix = matrix_path.copy()
#             print('length_decipher[j]=',length_decipher[j])
            for k in range(length_decipher[j]):
                # print('k=',k)
                length_path = decipher[j,k]
#                 print('length_path',length_path)
                for a in range(length_decipher[j]):
                    newmatrix[a, 0 : decipher[j,a]] = ind_tr[matrix_path[a, 0 : decipher[j, a]]]
                if j == 0:#fisrt move after having added a curl 
                    if length_path != 2:
                        (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,2:length_path][0]+1, 
                            newmatrix[k,2:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                        
                        ind_tr = index_before_to_later[ind_tr]
                        contribution += contri
                        # print('contri = ', contri)
                else:# not the first move
                    (K_evolve,index_before_to_later,contri ) = push_along_rl_cfgs(
                            K_evolve, 
                            newmatrix[k,0:length_path][0]+1, 
                            newmatrix[k,0:length_path], 
                            arr_cfgs, 
                            arr_R3_pairs_types, 
                            direction = -1
                            )
                    ind_tr = index_before_to_later[ind_tr]
                    contribution += contri
                    # print('contri = ', contri)
        # draw(K_evolve, 'i = {}'.format(i), ind_tr)
        if i < row_K-1:
            K_temp = infty_avancer(K_temp,-1)
    # draw(K_temp,'half')
    return contribution

