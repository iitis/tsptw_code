import dimod
from tsp_helpers import *
from math import floor, log

'''
Solution format [(0,1), (1,2),... ]
Route format (0,2,...,n-1)
'''


### Functions for route to array conversion ###

def route_to_array_ilp(route, ins):
    C, n = ins.C, ins.n
    earliest, latest = ins.earliest, ins.latest
    invalid = ins.invalid

    # First segment of the array representing the route
    route_array = soln_to_route_array_ilp(n, route, invalid)

    # Service and waiting times
    service, wait = service_wait(n, C, route, earliest)

    # Conversion to binary array
    s_arr = intvars_to_binary(service, C, latest, n, earliest, "service")
    w_arr = intvars_to_binary(wait, C, latest, n, earliest, "wait")

    # Conversion of slack variables to binary arrays
    arr24 = intvars_to_binary(eq24slack(n, service, wait, C, route_array, earliest), C, latest, n, earliest, "eq24")
    arr25 = intvars_to_binary(eq25slack(n, service, wait, C, route_array, latest, earliest), C, latest, n, earliest, "eq25")
    arr27 = eq2728_to_binary(eq27slack(n, C, service, latest, route, wait, earliest, invalid), C, earliest, latest, n,
                             invalid)
    arr28 = eq2728_to_binary(eq28slack(n, C, service, latest, route, wait, earliest, invalid), C, earliest, latest, n,
                             invalid)

    # Combining the arrays
    return route_array + s_arr + w_arr + arr24 + arr25 + arr27 + arr28


def array_to_soln_ilp(x, n, invalid):
    x = add_vars(x, invalid, n)
    soln = x[0:(n - 1) * n]
    return soln

# Converts route to ILP soln format
def route_to_soln_ilp(route):

    soln = [(route[i], route[i + 1]) for i in range(len(route) - 1)]
    soln.append((route[len(route) - 1], route[0]))

    return soln


# Creates part of the array corresponding to route
def soln_to_route_array_ilp(n, route, invalid):
    array = []

    soln = route_to_soln_ilp(route)
    for i, j in permutations(range(0, n), 2):
        if (i, j) in soln:
            array += [1]
        elif [i, j] not in invalid:
            array += [0]
    return array


# Given the route, generates service and waiting times
def service_wait(n, C, route, earliest):

    service, wait = [0] * n, [0] * n
    duration = 0
    soln = route_to_soln_ilp(route)

    for (i, j) in soln:
        duration += C[i][j]
        if duration < earliest[j]:
            service[j] = earliest[j]
            wait[j] = earliest[j] - duration
            duration = earliest[j]
        else:
            service[j] = duration
            wait[j] = 0
    service[0] = 0
    for i in range(len(service)):
        service[i] -= earliest[i]

    return service, wait

# Conversion of integer vars to binary array
def intvars_to_binary(array, C, latest, n, earliest, var):

    bservice_size = [0] * n
    for i in range(1, n):
        if var == "service":
            bservice_size[i] = floor(log(latest[i] - earliest[i], 2)) + 1
        elif var == "wait":
            if earliest[i] - C[0][i] > 0:
                bservice_size[i] = floor(log(earliest[i] - C[0][i], 2)) + 1
            else:
                bservice_size[i] = 0
        elif var == "eq24":
            bservice_size[i] = floor(log(latest[i], 2)) + 1
        elif var == "eq25":
            bservice_size[i] = floor(log(latest[i] - C[0][i], 2)) + 1

    service_arr = [0] * sum(bservice_size)

    k = 0
    for i in range(1, n):
        if var == "service":
            bi = mapint(bmap(latest[i] - earliest[i]), (array[i]))
        elif var == "eq24":
            bi = mapint(bmap(latest[i]), (array[i]))
        elif var == "eq25":
            bi = mapint(bmap(latest[i] - C[0][i]), (array[i]))
        elif var == "wait":
            if bservice_size[i] != 0:
                bi = mapint(bmap(earliest[i] - C[0][i]), (array[i]))
            else:
                continue

        for b in bi:
            service_arr[k] = int(b)
            k += 1

    return service_arr


# Binary conversion
def bmap(n):
    map = []
    k = 0
    while n > 2 ** k:
        map.append(2 ** k)
        n -= 2 ** k
        k += 1
    map.append(n)
    return map


# Map binary number to n variables
def mapint(map, n):

    if n == 0:
        return [0] * len(map)
    arr = [0] * len(map)
    k = len(map) - 1
    while n > 0:
        if n >= map[k]:
            n -= map[k]
            arr[k] = 1
        k -= 1
    return arr


# Generates integer slack variables
def eq24slack(n, service, wait, C, array, earliest):

    int24 = [0] * n
    for v in range(0, n):
        int24[v] = service[v] - wait[v] - C[0][v] * array[v - 1] + earliest[v]

    return int24


# Generates integer slack variables
def eq25slack(n, service, wait, C, array, latest, earliest):

    int25 = [0] * n
    for v in range(0, n):
        int25[v] = latest[v] - service[v] + wait[v] - (latest[v] - C[0][v]) * array[v - 1] - earliest[v]

    return int25

# Generates integer slack variables
def eq27slack(n, C, service, latest, route, wait, earliest, invalid):

    int27 = {}

    for u, v in permutations(range(1, n), 2):
        if [u, v] not in invalid:
            x = 0
            for i in range(1, len(route) - 1):
                if u == route[i] and v == route[i + 1]:
                    x = 1
            int27[(u, v)] = latest[u] - C[0][v] - service[u] + service[v] - wait[v] - (
                        latest[u] - C[0][v] + C[u][v]) * x - earliest[u] + earliest[v]

    return int27


# Generates integer slack variables
def eq28slack(n, C, service, latest, route, wait, earliest, invalid):

    int28 = {}

    for u, v in permutations(range(1, n), 2):
        if [u, v] not in invalid:
            x = 0
            for i in range(1, len(route) - 1):
                if u == route[i] and v == route[i + 1]:
                    x = 1
            int28[(u, v)] = latest[v] - earliest[u] - service[v] + wait[v] + service[u] - (
                        latest[v] - earliest[u] - C[u][v]) * x + earliest[u] - earliest[v]
    return int28


# Conversion of slack vars to binary arrs
def eq2728_to_binary(array, C, earliest, latest, n, invalid):

    sizes = {}

    for u, v in permutations(range(1, n), 2):
        if [u, v] not in invalid:
            sizes[(u, v)] = floor(log(-earliest[u] + latest[v] + latest[u] - C[0][v], 2)) + 1

    eq_arr = [0] * sum(sizes.values())

    k = 0
    for u, v in permutations(range(1, n), 2):
        if [u, v] not in invalid:
            bi = mapint(bmap(-earliest[u] + latest[v] + latest[u] - C[0][v]), int(array[(u, v)]))
            for b in bi:
                eq_arr[k] = int(b)
                k += 1

    return eq_arr


##############################################################################

### The functions defined below are for checking whether the solution is valid and optimal
def check_solution_ilp(soln, n, C, earliest, latest):

    cost = calc_cost_ilp(soln, n, C)
    valid = check_valid_soln_ilp(soln, n)
    if valid:
        order = []
        route = find_order(soln, order, 0, n)
        route = [0] + route[:-1]
        window = check_window(route, earliest, latest, C)
    else:
        route = []
        window = False

    return valid, route, window, cost

def find_order(soln, order, k, n):
    for j in range(k * (n - 1), k * (n - 1) + n - 1):
        if soln[j] == 1:
            if j % (n - 1) < k:
                order.append(j % (n - 1))
            else:
                order.append(j % (n - 1) + 1)
        if len(order) == n:
            return order
    return find_order(soln, order, order[-1], n)


# Checks whether an ilp solution is valid
def check_valid_soln_ilp(soln, n):
    arr = [0] * n
    k = 0
    for i, j in permutations(range(0, n), 2):

        if soln[k] == 1:
            arr[i] += 1
        k += 1
    return arr == [1] * n


# Adds missing variables to the array
def add_vars(x, invalid, n):
    k = 0
    for i, j in permutations(range(n), 2):
        if [i, j] in invalid:
            x = np.insert(x, k, -1)
        k += 1
    return x

####################################################################################

# Calculate cost of a solution
def calc_cost_ilp(soln, n, C):

    k = 0
    cost = 0
    for i, j in permutations(range(0, n), 2):
        if soln[k] == 1:
            cost += C[i][j]
        k += 1
    return cost


# Evaluates samples in the sampleset
def evaluate_sampleset_ilp(ins, sampleset):

    results = []
    for datum in sampleset.data(fields=['sample', 'energy']):
        x = dimod.sampleset.as_samples(datum.sample)[0][0]
        energy = datum.energy

        n = ins.n
        x = add_vars(x, ins.invalid, n)
        route_array = x[0:(n - 1) * n]

        valid, route, window, cost = check_solution_ilp(route_array, n, ins.C, ins.earliest, ins.latest)
        results.append([valid, window, cost, energy])

    return results


'''
def cplex(n, C, earliest, latest):
    qp = QuadraticProgram()

    add_variables2(qp, n, earliest, latest, C)

    ham_tour(qp, n)
    cost(qp, n, C)
    time_window(qp, n, C, earliest, latest)

    return qp


def check_cplex():
    dir = "instances"
    for instance in os.listdir(dir):
        instance = instance[:-4]
        ins = Instance(f"{instance}", "instances")
        print()
        print("------------", instance, "-----------")
        print(ins.C, ins.earliest, ins.latest)

        print("Optimal solution is:", brute_force_tsptw(ins.C, ins.latest, ins.earliest))

        problem = cplex(ins.n, ins.C, ins.earliest, ins.latest)
        optimizer = CplexOptimizer() if CplexOptimizer.is_cplex_installed() else None
        if optimizer: result = optimizer.solve(problem)
        print(result)

        valid, route, window, cost = check_solution(result, ins.n, ins.C, ins.earliest, ins.latest)
        print("Route is ", route, ". Time windows is", window, ". Cost is ", cost)
'''
