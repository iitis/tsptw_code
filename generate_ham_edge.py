import os
from itertools import *
import dimod
from instance_edge import *
from collections import defaultdict


# Generates hamiltonian tour contsraints
def ham_tour(instance):
    Map, elist = instance.Map, instance.elist
    G, G_0 = instance.G, instance.G_0
    n = instance.n

    Q = defaultdict(lambda: 0)

    # Term 1
    for u in G.neighbors(0):
        Q[Map.get(((0, u), 1)), Map.get(((0, u), 1))] -= 1
        for v in G.neighbors(0):
            if u < v:
                Q[Map.get(((0, u), 1)), Map.get(((0, v), 1))] += 2

    # Term 2
    for t, e in product(range(2, n), elist):
        Q[Map.get((e, t)), Map.get((e, t))] -= 1
        for f in elist:
            if e < f:
                Q[Map.get((e, t)), Map.get((f, t))] += 2
    # Term 3
    for u in G.successors(0):
        Q[Map.get(((u, 0), n)), Map.get(((u, 0), n))] -= 1
        for v in G.successors(0):
            if u < v:
                Q[Map.get(((u, 0), n)), Map.get(((v, 0), n))] += 2
    # Term 4
    for u in range(1, n):
        for v, t in product(G_0.neighbors(u), range(2, n)):
            Q[Map.get(((u, v), t)), Map.get(((u, v), t))] -= 1
            for w in G_0.neighbors(u):
                if v < w:
                    Q[Map.get(((u, v), t)), Map.get(((u, w), t))] += 2
                for d in range(t, n):
                    if v == w and t < d:
                        Q[Map.get(((u, v), t)), Map.get(((u, v), d))] += 2
                    elif t < d:
                        Q[Map.get(((u, v), t)), Map.get(((u, w), d))] += 2
    for u in G.successors(0):
        Q[Map.get(((u, 0), n)), Map.get(((u, 0), n))] -= 1
        for v, t in product(G_0.neighbors(u), range(2, n)):
            Q[Map.get(((u, 0), n)), Map.get(((u, v), t))] += 2

    # Term 5
    for v in G.neighbors(0):
        for w in G_0.neighbors(v):
            Q[Map.get(((0, v), 1)), Map.get(((v, w), 2))] -= 1

    # Term 6
    for t, v in product(range(2, n - 1), range(1, n)):
        for w, u in product(G_0.neighbors(v), G_0.successors(v)):
            if w != u:
                Q[Map.get(((u, v), t)), Map.get(((v, w), t + 1))] -= 1

    # Term 7
    for v in G.successors(0):
        for u in G_0.successors(v):
            Q[Map.get(((u, v), n - 1)), Map.get(((v, 0), n))] -= 1

    return Q


# Generates the cost of the tour
def cost(instance):
    Map, elist = instance.Map, instance.elist
    n, C, G = instance.n, instance.C, instance.G

    Q = defaultdict(lambda: 0)

    for v in G.neighbors(0):
        Q[Map.get(((0, v), 1)), Map.get(((0, v), 1))] += C[0][v]
    for e, t in product(elist, range(2, n)):
        Q[Map.get((e, t)), Map.get((e, t))] += C[e[0]][e[1]]
    for v in G.successors(0):
        Q[Map.get(((v, 0), n)), Map.get(((v, 0), n))] += C[v][0]

    return Q


# Add latest finish time constraints for each i
def latest_cons(instance, i):
    G, C, elist, latest = instance.G, instance.C, instance.elist, instance.latest
    Map, Map_s, Map_w = instance.Map, instance.Map_s, instance.Map_w
    D, W = instance.D_arr[i - 1], instance.W_arr[i - 1]
    W_arr = instance.W_arr

    Q = defaultdict(lambda: 0)

    if i == 1:
        # a term
        for v in G.neighbors(0):
            Q[Map.get(((0, v), 1)), Map.get(((0, v), 1))] += (C[0][v] - latest[v]) ** 2
            for u in G.neighbors(0):
                if u > v:
                    Q[Map.get(((0, v), 1)), Map.get(((0, u), 1))] += 2 * (C[0][v] - latest[v]) * (C[0][u] - latest[u])
        # b term
        for d in range(D):
            Q[Map_s.get((d, 1)), Map_s.get((d, 1))] += (2 ** d) * (2 ** d)
            for s in range(d + 1, D):
                Q[Map_s.get((d, 1)), Map_s.get((s, 1))] += 2 * (2 ** d) * (2 ** s)
            # 2ab
            for v in G.neighbors(0):
                Q[Map.get(((0, v), 1)), Map_s.get((d, 1))] += 2 * (2 ** d) * (C[0][v] - latest[v])
    else:
        # a term
        term_count = defaultdict(lambda: 0)
        for t in (range(1, i)):
            for k in range(W_arr[t-1]):
                Q[Map_w.get((k, t)), Map_w.get((k, t))] += (2 ** k) * (2 ** k)
                term_count['a'] += 1
                for g in (range(t, i)):
                    for l in range(W_arr[g-1]):
                        if g == t and k < l:
                            Q[Map_w.get((k, t)), Map_w.get((l, t))] += 2 * (2 ** k) * (2 ** l)
                            term_count['a'] += 1
                        elif g > t:
                            Q[Map_w.get((k, t)), Map_w.get((l, g))] += 2 * (2 ** k) * (2 ** l)
                            term_count['a'] += 1
                # ab
                for v in G.neighbors(0):
                    Q[Map_w.get((k, t)), Map.get(((0, v), 1))] += 2 * (2 ** k) * (C[0][v])
                    term_count['ab'] += 1
                # ac
                for tau, e in product(range(2, i + 1), elist):
                    Q[Map_w.get((k, t)), Map.get((e, tau))] += 2 * (2 ** k) * (C[e[0]][e[1]])
                    term_count['ac'] += 1
                # ad
                for e in elist:
                    Q[Map_w.get((k, t)), Map.get((e, i))] += -2 * (2 ** k) * (latest[e[1]])
                    term_count['ad'] += 1
                # ae
                for d in range(D):
                    Q[Map_w.get((k, t)), Map_s.get((d, i))] += 2 * (2 ** k) * (2 ** d)
                    term_count['ae'] += 1
        # b term
        for v in G.neighbors(0):
            Q[Map.get(((0, v), 1)), Map.get(((0, v), 1))] += (C[0][v]) ** 2  # a^2
            term_count['b'] += 1
            for u in G.neighbors(0):
                if u > v:
                    Q[Map.get(((0, v), 1)), Map.get(((0, u), 1))] += 2 * (C[0][v]) * (C[0][u])
                    term_count['b'] += 1
            # bc
            for t, e in product(range(2, i + 1), elist):
                Q[Map.get(((0, v), 1)), Map.get((e, t))] += 2 * (C[0][v]) * (C[e[0]][e[1]])
                term_count['bc'] += 1
            # bd
            for e in elist:
                Q[Map.get(((0, v), 1)), Map.get((e, i))] += -2 * (C[0][v]) * (latest[e[1]])
                term_count['bd'] += 1
            # be
            for d in range(D):
                Q[Map.get(((0, v), 1)), Map_s.get((d, i))] += 2 * (C[0][v]) * (2 ** d)
                term_count['be'] += 1
        # c term
        for t, e in product(range(2, i + 1), elist):
            Q[Map.get((e, t)), Map.get((e, t))] += (C[e[0]][e[1]]) ** 2
            term_count['c'] += 1
            for g, f in product(range(t, i + 1), elist):
                if g == t and e < f:
                    Q[Map.get((e, t)), Map.get((f, t))] += 2 * (C[e[0]][e[1]]) * (C[f[0]][f[1]])
                    term_count['c'] += 1
                elif g > t:
                    Q[Map.get((e, t)), Map.get((f, g))] += 2 * (C[e[0]][e[1]]) * (C[f[0]][f[1]])
                    term_count['c'] += 1
            # cd
            for ee in elist:
                Q[Map.get((e, t)), Map.get((ee, i))] += -2 * (C[e[0]][e[1]]) * (latest[ee[1]])
                term_count['cd'] += 1
            # ce
            for d in range(D):
                Q[Map.get((e, t)), Map_s.get((d, i))] += 2 * (C[e[0]][e[1]]) * (2 ** d)
                term_count['ce'] += 1
        # d term
        for e in elist:
            Q[Map.get((e, i)), Map.get((e, i))] += (latest[e[1]]) ** 2
            term_count['d'] += 1
            for f in elist:
                if e < f:
                    Q[Map.get((e, i)), Map.get((f, i))] += 2* (latest[e[1]]) * latest[f[1]]
                    term_count['d'] += 1
            # de
            for d in range(D):
                Q[Map.get((e, i)), Map_s.get((d, i))] += -2 * latest[e[1]] * (2 ** d)
                term_count['de'] += 1
        # e term
        for d in range(D):
            Q[Map_s.get((d, i)), Map_s.get((d, i))] += (2 ** d) * (2 ** d)
            term_count['e'] += 1
            for s in range(d + 1, D):
                Q[Map_s.get((d, i)), Map_s.get((s, i))] += 2 * (2 ** d) * (2 ** s)
                term_count['e'] += 1
    # for key in term_count.keys():
    #    print(key,term_count[key])
    return Q


# Add earliest finish time constraints for each i
def earliest_cons(instance, i):
    G, C, elist, earliest = instance.G, instance.C, instance.elist, instance.earliest
    Map, Map_u, Map_w = instance.Map, instance.Map_u, instance.Map_w
    L, W = instance.L, instance.W_arr[i - 1]
    W_arr = instance.W_arr
    Q = defaultdict(lambda: 0)

    if i == 1:
        # a
        for k in range(W):
            Q[Map_w.get((k, 1)), Map_w.get((k, 1))] += (2 ** k) * (2 ** k)
            for l in range(k + 1, W):
                Q[Map_w.get((k, 1)), Map_w.get((l, 1))] += 2 * (2 ** k) * (2 ** l)
            # ab
            for v in G.neighbors(0):
                Q[Map.get(((0, v), 1)), Map_w.get((k, 1))] += 2 * (C[0][v] - earliest[v]) * (2 ** k)
            # ac
            for d in range(L):
                Q[Map_w.get((k, 1)), Map_u.get((d, 1))] += -2 * (2 ** k) * (2 ** d)
        # b
        for v in G.neighbors(0):
            Q[Map.get(((0, v), 1)), Map.get(((0, v), 1))] += (C[0][v] - earliest[v]) ** 2
            for u in G.neighbors(0):
                if u > v:
                    Q[Map.get(((0, v), 1)), Map.get(((0, u), 1))] += 2 * (C[0][v] - earliest[v]) * (
                                C[0][u] - earliest[u])
            for d in range(L):
                Q[Map.get(((0, v), 1)), Map_u.get((d, 1))] += -2 * (2 ** d) * (C[0][v] - earliest[v])

        # c
        for d in range(L):
            Q[Map_u.get((d, 1)), Map_u.get((d, 1))] += (2 ** d) * (2 ** d)
            for s in range(d + 1, L):
                Q[Map_u.get((d, 1)), Map_u.get((s, 1))] += 2 * (2 ** d) * (2 ** s)

    else:
        # a term
        term_count = defaultdict(lambda: 0)
        for t in (range(1, i + 1)):
            for k in range(W_arr[t-1]):
                Q[Map_w.get((k, t)), Map_w.get((k, t))] +=  (2 ** k) * (2 ** k)
                term_count['a'] += 1
                for g in range(t, i + 1):
                    for l in range(W_arr[g-1]):
                        if g == t and k < l:
                            Q[Map_w.get((k, t)), Map_w.get((l, g))] += 2 * (2 ** k) * (2 ** l)
                            term_count['a'] += 1
                        elif g > t:
                            Q[Map_w.get((k, t)), Map_w.get((l, g))] += 2 * (2 ** k) * (2 ** l)
                            term_count['a'] += 1

                # ab
                for v in G.neighbors(0):
                    Q[Map_w.get((k, t)), Map.get(((0, v), 1))] += 2 * (2 ** k) * (C[0][v])
                    term_count['ab'] += 1
                # ac
                for tau, e in product(range(2, i + 1), elist):
                    Q[Map_w.get((k, t)), Map.get((e, tau))] += 2 * (2 ** k) * (C[e[0]][e[1]])
                    term_count['ac'] += 1
                # ad
                for e in elist:
                    Q[Map_w.get((k, t)), Map.get((e, i))] += -2 * (2 ** k) * (earliest[e[1]])
                    term_count['ad'] += 1
                # ae
                for d in range(L):
                    Q[Map_w.get((k, t)), Map_u.get((d, i))] += -2 * (2 ** k) * (2 ** d)
                    term_count['ae'] += 1
        # b: xc term
        for v in G.neighbors(0):
            Q[Map.get(((0, v), 1)), Map.get(((0, v), 1))] += (C[0][v]) ** 2  # a^2
            term_count['b'] += 1
            for u in G.neighbors(0):
                if u > v:
                    Q[Map.get(((0, v), 1)), Map.get(((0, u), 1))] += 2 * (C[0][v]) * (C[0][u])
                    term_count['b'] += 1

            # bc
            for t, e in product(range(2, i + 1), elist):
                Q[Map.get(((0, v), 1)), Map.get((e, t))] += 2 * (C[0][v]) * (C[e[0]][e[1]])
                term_count['bc'] += 1
            # bd
            for e in elist:
                Q[Map.get(((0, v), 1)), Map.get((e, i))] += -2 * (C[0][v]) * (earliest[e[1]])
                term_count['bd'] += 1
            # be
            for d in range(L):
                Q[Map.get(((0, v), 1)), Map_u.get((d, i))] += -2 * (C[0][v]) * (2 ** d)
                term_count['be'] += 1
        # c term
        for t, e in product(range(2, i + 1), elist):
            Q[Map.get((e, t)), Map.get((e, t))] += (C[e[0]][e[1]]) ** 2
            term_count['c'] += 1
            for g, f in product(range(t, i + 1), elist):
                if g == t and e < f:
                    Q[Map.get((e, t)), Map.get((f, t))] += 2 * (C[e[0]][e[1]]) * (C[f[0]][f[1]])
                    term_count['c'] += 1
                elif g > t:
                    Q[Map.get((e, t)), Map.get((f, g))] += 2 * (C[e[0]][e[1]]) * (C[f[0]][f[1]])
                    term_count['c'] += 1
            # cd
            for ee in elist:
                Q[Map.get((e, t)), Map.get((ee, i))] += -2 * (C[e[0]][e[1]]) * (earliest[ee[1]])
                term_count['cd'] += 1
            # ce
            for d in range(L):
                Q[Map.get((e, t)), Map_u.get((d, i))] += -2 * (C[e[0]][e[1]]) * (2 ** d)
                term_count['ce'] += 1
        # d term
        for e in elist:
            Q[Map.get((e, i)), Map.get((e, i))] += (earliest[e[1]]) ** 2
            term_count['d'] += 1
            for f in elist:
                if e < f:
                    Q[Map.get((e, i)), Map.get((f, i))] += 2*(earliest[e[1]]) * earliest[f[1]]
                    term_count['d'] += 1
            # de
            for d in range(L):
                Q[Map.get((e, i)), Map_u.get((d, i))] += 2 * earliest[e[1]] * (2 ** d)
                term_count['de'] += 1
        # e s_term
        for d in range(L):
            Q[Map_u.get((d, i)), Map_u.get((d, i))] += (2 ** d) * (2 ** d)
            term_count['e'] += 1
            for s in range(d + 1, L):
                Q[Map_u.get((d, i)), Map_u.get((s, i))] += 2 * (2 ** d) * (2 ** s)
                term_count['e'] += 1
        # for key in term_count.keys():
        #    print(key,term_count[key])
    return Q


# Generate Time Windows contraints
def time_window(instance):
    n = instance.n
    H_t, H = defaultdict(lambda: 0), defaultdict(lambda: 0)
    for i in range(1,n):
        H_t_i_1 = latest_cons(instance, i)
        H_t_i_2 = earliest_cons(instance, i)

        for k in H_t_i_1.keys():
            if k in H_t_i_2.keys():
                H_t[k] +=  H_t_i_1.get(k, 0) + H_t_i_2.get(k, 0)
            else:
                H_t[k] += H_t_i_1.get(k, 0)

        for k in H_t_i_2.keys():
            if k not in H_t_i_1.keys():
                H_t[k] += H_t_i_2.get(k, 0)

    return H_t


# Removes variables corresponding to invalid edges
def remove_zeros(Q, zeros):
    newQ = {}
    if zeros != []:
        for key, value in Q.items():
            if all(f"({z}," not in str(key) and f", {z})" not in str(key) for z in zeros):
                newQ[key] = value
    else:
        newQ = Q
    return newQ


# Combines the constraints and returns binary quadratic model
def generate_bqm_edge(instance, P1,P2):
    H_c = ham_tour(instance)
    H_m = cost(instance)
    H_t = time_window(instance)

    H = {}
    for k in chain(H_c.keys(), H_m.keys(), H_t.keys()):
        H[k] = P1 * H_c.get(k, 0) +  H_m.get(k, 0) + P2 * H_t.get(k, 0)

    H = remove_zeros(H, instance.invalid)
    n = instance.n

    return dimod.BinaryQuadraticModel.from_qubo(H, offset = (4 + (n-2) + (n-1) + (n-1-2))*P1 )


# Returns ising formulatiom
def ising_edge(instance, p1, p2):

    P1 = p1*instance.C.max()
    P2 = (1/p2)*instance.C.max()

    bqm = generate_bqm_edge(instance, P1, P2)
    return Instance_edge.convert_and_scale(bqm)

