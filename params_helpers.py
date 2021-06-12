
# Finds the optmial cost
def find_opt_cost(data):
    min_cost = data[data.A == -1]['cost'].values[0]
    return min_cost


# Finds the energy corresponding to optimal route
def find_opt_energy(data):
    opt_energy = data[data.A == -1]['energy'].values[0]
    return opt_energy


# Returns true if optimal route is observed
def opt_found(data):
    min_cost = find_opt_cost(data)
    data_d = data[(data.cost == min_cost) & (data.valid == True) & (data.windows == True) & (data.A != -1)]
    return len(data_d) > 0


# Calculates the probability of samples encoding optimal route
def probability(data):
    min_cost = find_opt_cost(data)
    data_d = data[(data.cost == min_cost) & (data.valid == True) & (data.windows == True) & (data.A != -1)]
    return len(data_d) / (len(data) - 1)


# Calculates the probability of samples encoding optimal route, returns 0 if the lowest energy sample is not optimal
def min_prob(data):
    min_cost = find_opt_cost(data)
    if data.valid.iloc[0] == True & data.windows.iloc[0] == True and data.cost.iloc[0] == min_cost:
        return probability(data)
    else:
        return 0
