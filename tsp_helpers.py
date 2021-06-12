from itertools import permutations
import numpy as np
import edge_helpers as eh

'''
Route format (0,2,...,n-1)
'''

# Conversion between binary and spin variables
def spin_to_binary(spin):
    return [1 if s == 1 else 0 for s in spin]


# Conversion between binary and spin variables
def binary_to_spin(binary):
    return [1 if b == 1 else -1 for b in binary]


# Check if time windows is obeyed
def check_window(route, earliest, latest, C):
    c = 0
    for i in range(len(route) - 1):
        c += C[route[i]][route[i + 1]]
        if c > latest[route[i + 1]]:
            return False
        elif c < earliest[route[i + 1]]:
            c = earliest[route[i + 1]]
    return True


# Checks if the tsp is metric
def is_metric_tsp(C,n):
    for i,j,k in permutations(range(n),3):
        if C[i][j]+C[j][k] <= C[i][k]:
            return False
    return True


# Find best TSP solution
def brute_force_tsp(C):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    costs = [eh.calc_cost_edge(C, eh.route_to_soln_edge(route)) for route in perm]
    return perm[np.argmin(costs)], min(costs)


# Find all routes
def brute_force_tsptw_all(C, latest, earliest):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    all_costs = [(eh.calc_cost_edge(C, eh.route_to_soln_edge(route)),route,check_window(route, earliest, latest, C)) for route in perm]
    return all_costs


# Find best TSPTW solution
def brute_force_tsptw(C, latest, earliest):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    costs = [1000 if not check_window(route, earliest, latest, C) else eh.calc_cost_edge(C, eh.route_to_soln_edge(route)) for route in perm]
    return perm[np.argmin(costs)], min(costs)
