import numpy as np
import dimod

import tsp_helpers as tsp

'''
Solution format [((0,1),1), ((1,2),2),... ]
Route format (0,2,...,n-1)
'''

### The first four functions are used for route to array conversion ###

def route_to_array_edge(route,ins):

    Map, Map_s, Map_w, Map_u = ins.Map, ins.Map_s, ins.Map_w, ins.Map_u
    earliest, latest, C = ins.earliest, ins.latest, ins.C

    arr = np.zeros(len(Map) + len(Map_s) + len(Map_w) + len(Map_u))

    for u, v, t in zip(route, route[1:], range(1, len(route))):
        arr[Map[((u, v), t)]] = 1
    arr[Map[((route[-1], route[0]), t + 1)]] = 1
    s, w, u = soln_to_spin(route, earliest, latest, C)

    for pairs in s:
        arr[Map_s[pairs]] = 1
    for pairw in w:
        arr[Map_w[pairw]] = 1
    for pairu in u:
        arr[Map_u[pairu]] = 1

    return arr


# Converts binary list of variables to unary
def binary_to_unary(list, n):
    ulist = np.zeros(n)
    for l in list:
        ulist[l[1]] += 2 ** l[0]
    return ulist


# Generates spin assignment from solution
def soln_to_spin(soln, earliest, latest, C):
    c = 0
    slack, waiting, slackl = [], [], []

    for i in range(len(soln) - 1):
        c += C[soln[i]][soln[i + 1]]
        if c > latest[soln[i + 1]]:
            return False
        else:
            slack.append(latest[soln[i + 1]] - c)
        if c < earliest[soln[i + 1]]:
            waiting.append(earliest[soln[i + 1]] - c)
            slackl.append(0)
            c = earliest[soln[i + 1]]
        else:
            waiting.append(0)
            slackl.append(c - earliest[soln[i + 1]])
    return var_to_blist(slack), var_to_blist(waiting), var_to_blist(slackl)

# Converts variables into binary list
def var_to_blist(slack):
    s = []
    for i in range(len(slack)):
        slack = np.array(slack, dtype=int)
        b = np.binary_repr(slack[i])
        for j in range(len(b)):
            if b[j] != '0':
                s.append((len(b) - j - 1, i + 1))
    return s


#Converts array to integer vars
def array_to_vars_edge(x, Map_s, Map_w, Map_u, n):

    k = (n - 1) ** 3 - 2 * (n - 1) ** 2 + 3 * (n - 1)

    s_list = [list(Map_s)[i] for i in range(len((Map_s))) if x[k + i] == 1]
    sl = binary_to_unary(s_list,n)
    k += len(Map_s)
    w_list = [list(Map_w)[i] for i in range(len((Map_w))) if x[k + i] == 1]
    wl = binary_to_unary(w_list,n)
    k += len(Map_w)
    u_list = [list(Map_u)[i] for i in range(len((Map_u))) if x[k + i] == 1]
    ul = binary_to_unary(u_list,n)

    return sl,wl,ul


### The functions defined below are for checking whether the solution is valid and optimal
def check_solution_edge(soln, n, C, earliest, latest):

    cost = calc_cost_edge(C, soln)
    valid = check_valid_soln_edge(n, soln)
    route = soln_to_route_edge(soln)
    window = tsp.check_window(route, earliest, latest, C)
    return valid, route, window, cost


# Check validity of solution
def check_valid_soln_edge(n, soln):
    if n != len(soln):
        return False
    check = np.zeros(n)
    for i in range(len(soln) - 1):
        if soln[i][0][1] != soln[i + 1][0][0]:
            return False
        check[soln[i][0][0]] += 1
    check[soln[-1][0][0]] += 1
    if soln[0][0][0] != 0 or soln[-1][0][1] != 0:
        return False
    if not (check == np.ones(n)).all():
        return False
    return True


# Converts array to edge soln format
def array_to_soln_edge(x, Map, invalid):

    for z in invalid:
        x = np.insert(x, z , -1)

    return [list(Map)[i] for i in range(len(Map)) if x[i] == 1]


# Converts route format to soln
def route_to_soln_edge(route):
    soln = [((u, v), t) for u, v, t in zip(route, route[1:], range(1, len(route)))]
    soln.append(((route[-1], route[0]), len(route)))
    return soln


# Converts solution to route
def soln_to_route_edge(soln):
    return [soln[i][0][0] for i in range(len(soln))]


# Evaluates samples in the sampleset
def evaluate_sampleset_edge(ins, sampleset):

    results = []
    for datum in sampleset.data(fields=['sample', 'energy']):
        spin_sol = dimod.sampleset.as_samples(datum.sample)[0][0]
        energy = datum.energy
        x = tsp.spin_to_binary(spin_sol)

        soln = array_to_soln_edge(x, ins.Map, ins.invalid)
        cost = calc_cost_edge(ins.C, soln)
        valid = check_valid_soln_edge(ins.n, soln)
        windows = tsp.check_window(soln_to_route_edge(soln), ins.earliest, ins.latest, ins.C)

        results.append([valid, windows, cost, energy])

    return results

# Calculate cost of a solution
def calc_cost_edge(C, soln):
    return sum(C[e[0][0]][e[0][1]] for e in soln)