from itertools import permutations
import numpy as np

'''
Route format (0,2,...,n-1)
Solution format [((0,1),1), ((1,2),2),... ]
'''

# Converts solution to route
def soln_to_route(soln):
    return [soln[i][0][0] for i in range(len(soln))]

# Converts route format to soln
def route_to_soln(route):
    soln = [((u, v), t) for u, v, t in zip(route, route[1:], range(1, len(route)))]
    soln.append(((route[-1], route[0]), len(route)))
    return soln

# Calculate cost of a solution
def calc_cost(C, soln):
    return sum(C[e[0][0]][e[0][1]] for e in soln)

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


def is_metric_tsp(C,n):
    for i,j,k in permutations(range(n),3):
        if C[i][j]+C[j][k] <= C[i][k]:
            return False
    return True

# Find best TSP solution
def brute_force_tsp(C):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    costs = [calc_cost(C, route_to_soln(route)) for route in perm]
    return perm[np.argmin(costs)], min(costs)

# Find all routes
def brute_force_tsptw_all(C, latest, earliest):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    all_costs = [(calc_cost(C, route_to_soln(route)),route,check_window(route, earliest, latest, C)) for route in perm]
    return all_costs

# Find best TSPTW solution
def brute_force_tsptw(C, latest, earliest):
    perm = [[0] + list(route) for route in permutations(range(1, len(C)))]
    costs = [1000 if not check_window(route, earliest, latest, C) else calc_cost(C, route_to_soln(route)) for route in perm]
    return perm[np.argmin(costs)], min(costs)
