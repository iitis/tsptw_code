from collections import defaultdict

import dimod
from qiskit.optimization import QuadraticProgram
from itertools import *
from qiskit.optimization.converters import LinearEqualityToPenalty, IntegerToBinary, InequalityToEquality
from instance_ilp import Instance_ilp


# Variables for the route and cost hamiltonian
def add_variables1(qp, n, invalid):
    for i, j in permutations(range(0, n), 2):
        if ([i, j] not in invalid):
            qp.binary_var(f"x{i}{j}")


# Variables for the time windows hamiltonian
def add_variables2(qp, n, earliest, latest, C, invalid):
    for i, j in permutations(range(0, n), 2):
        if ([i, j] not in invalid):
            qp.binary_var(f"x{i}{j}")

    for i in range(1, n):
        qp.integer_var(earliest[i], latest[i], f"s{i}")

    for i in range(1, n):
        if earliest[i] - C[0][i] > 0:
            qp.integer_var(0, earliest[i] - C[0][i], f"w{i}")


# Route constraint1
def eq17(qp, n, invalid):
    for i in range(0, n):
        eq17 = {}
        for j in range(0, n):
            if i != j and ([i, j] not in invalid):
                eq17[f"x{i}{j}"] = 1
        qp.linear_constraint(linear=eq17, sense='E', rhs=1, name=f"eq17{i}")


# Route constraint2
def eq18(qp, n, invalid):
    for j in range(0, n):
        eq18 = {}
        for i in range(0, n):
            if i != j and ([i, j] not in invalid):
                eq18[f"x{i}{j}"] = 1
        qp.linear_constraint(linear=eq18, sense='E', rhs=1, name=f"eq18{j}")


# Time windows1
def eq24(qp, C, n, earliest):
    for i in range(1, n):
        if earliest[i] - C[0][i] > 0:
            qp.linear_constraint(linear={f"w{i}": -1, f"s{i}": 1, f"x{0}{i}": -C[0][i]}, sense='GE', rhs=0,
                                 name=f"eq24{i}")
        else:
            qp.linear_constraint(linear={f"s{i}": 1, f"x{0}{i}": -C[0][i]}, sense='GE', rhs=0, name=f"eq24{i}")


# Time windows2
def eq25(qp, C, n, latest, earliest):
    for i in range(1, n):
        if earliest[i] - C[0][i] > 0:
            qp.linear_constraint(linear={f"w{i}": -1, f"s{i}": 1, f"x{0}{i}": latest[i] - C[0][i]}, sense='LE',
                                 rhs=latest[i],
                                 name=f"eq25{i}")
        else:
            qp.linear_constraint(linear={f"s{i}": 1, f"x{0}{i}": latest[i] - C[0][i]}, sense='LE',
                                 rhs=latest[i],
                                 name=f"eq25{i}")


# Time windows3
def eq27(qp, C, earliest, latest, n, invalid):
    for i, j in permutations(range(1, n), 2):
        if ([i, j] not in invalid):

            cons = defaultdict(lambda: 0)
            cons[f"s{i}"] = 1
            if earliest[j] - C[0][j] > 0:
                cons[f"w{j}"] = 1
            cons[f"s{j}"] = -1

            cons[f"x{i}{j}"] = latest[i] - C[0][j] + C[i][j]
            qp.linear_constraint(linear=cons, sense='LE', rhs=latest[i] - C[0][j],
                                 name=f"eq27{i}{j}")


# Time windows4
def eq28(qp, C, earliest, latest, n, invalid):
    for i, j in permutations(range(1, n), 2):
        if ([i, j] not in invalid):
            cons = defaultdict(lambda: 0)
            cons[f"s{i}"] = -1
            if earliest[j] - C[0][j] > 0:
                cons[f"w{j}"] = -1
            cons[f"s{j}"] = 1

            cons[f"x{i}{j}"] = latest[j] - earliest[i] - C[i][j]
            qp.linear_constraint(linear=cons, sense='LE', rhs=latest[j] - earliest[i],
                                 name=f"eq28{i}{j}")


# Generates route constraints
def ham_tour(qp, n, invalid):
    eq17(qp, n, invalid)
    eq18(qp, n, invalid)


# Generates objective function
def cost(qp, n, C, invalid):
    T1 = defaultdict(lambda: 0)
    for i, j in permutations(range(n), 2):
        if [i, j] not in invalid:
            T1[f"x{i}{j}"] =  C[i][j]
    qp.minimize(linear=T1)


# Generates time windows
def time_window(qp, n, C, earliest, latest, invalid):
    eq24(qp, C, n, earliest)
    eq25(qp, C, n, latest, earliest)
    eq27(qp, C, earliest, latest, n, invalid)
    eq28(qp, C, earliest, latest, n, invalid)


# Hamiltonian for route and cost
def ham1(n, C, P1, invalid):
    qp = QuadraticProgram()
    add_variables1(qp, n, invalid)
    ham_tour(qp, n, invalid)
    cost(qp, n,  C, invalid)

    lineq2penalty = LinearEqualityToPenalty(penalty=P1)
    qubo1 = lineq2penalty.convert(qp)

    return qubo1


# Hamiltonian for time windows
def ham2(n, C, earliest, latest, E, invalid):
    qp = QuadraticProgram()
    add_variables2(qp, n, earliest, latest, C, invalid)
    time_window(qp, n, C, earliest, latest, invalid)

    ineq2eq = InequalityToEquality()
    qp = ineq2eq.convert(qp)

    int2bin = IntegerToBinary()
    qp = int2bin.convert(qp)

    lineq2penalty = LinearEqualityToPenalty(penalty=E)
    qub = lineq2penalty.convert(qp)

    return qub


# Combines hamiltonians and returns binary quadratic model
def generate_bqm_ilp(ins,P1,P2):
    n = ins.n
    invalid = ins.invalid
    earliest , latest = ins.earliest, ins.latest
    C = ins.C

    qubo1 = ham1(n, C, P1, invalid)
    qubo2 = ham2(n, C, earliest, latest, P2, invalid)

    L = {}
    for k in chain(qubo1.objective.linear.to_dict().keys(), qubo2.objective.linear.to_dict().keys()):
        L[k] = qubo1.objective.linear.to_dict().get(k, 0) + qubo2.objective.linear.to_dict().get(k, 0)

    Q = {}
    for k in chain(qubo1.objective.quadratic.to_dict().keys(), qubo2.objective.quadratic.to_dict().keys()):
        Q[k] = qubo1.objective.quadratic.to_dict().get(k, 0) + qubo2.objective.quadratic.to_dict().get(k, 0)

    of = qubo1.objective.constant + qubo2.objective.constant
    q, of = dimod.AdjVectorBQM(L, Q, of, vartype=dimod.BINARY).to_qubo()

    return dimod.BQM.from_qubo(q, of)


# Return ising formulation
def ising_ilp(ins,p1,p2):
    P1 = p1 * ins.C.max()
    P2 = (1/p2) * ins.C.max()

    bqm = generate_bqm_ilp(ins,P1,P2)
    h, J, of = Instance_ilp.convert_and_scale(bqm)

    return h, J, of





