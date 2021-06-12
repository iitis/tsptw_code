import argparse
import os
from tsp_helpers import *

# Generates the time windows where the window range is selected randomly between l_t and u_t
def generate_tw(n, l_t, u_t, C):
    w = np.zeros((n, 2))
    for i in range(1, n):
        w[i][0] = np.random.randint(C[0][i], l_t)
    for i in range(1, n):
        w[i][1] = np.random.randint(l_t+1, u_t)
    return w

# Generates the cost matrix where the costs are between l and u
def generate_cost(n,l,u,sym):
    metric = False
    while not metric:
        if sym:
            c = np.random.randint(l, u, size=(n, n))
            C = np.tril(c) + np.tril(c, -1).T
        else:
            C = np.random.randint(l, u, size=(n, n))
            np.fill_diagonal(C, 0)
        metric = is_metric_tsp(C, n)
    return C

# Generates a tsptw instance where the underlying graph is complete
def complete(dir,n_instances, sym, n, l, u, l_t, u_t, txt):
    i = 0
    tsym = 1 if sym else 0
    while i < n_instances:
        C = generate_cost(n,l,u,sym)
        w = generate_tw(n, l_t, u_t, C)

        _, cost = brute_force_tsptw(C, w.T[1],w.T[0])
        _ , costt = brute_force_tsp(C)

        if(i < n_instances/2) and (cost == 1000 or costt != cost):
            continue
        elif(i >= n_instances/2) and (cost == 1000 or costt == cost):
            continue
        if txt:
            file_name = f"{dir}_txt/complete_{tsym}_{n}_{i}"
            f = open(file_name, "w")
            print(n, file=f)
            for row in C:
                print(*row, file=f)
            for row in w:
                print(*row, file=f)
        file_name = f"{dir}/complete_{tsym}_{n}_{i}"
        np.savez(file_name, n=n, weight=C, window=w)
        i += 1

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=str,
                        help="Name of the folder to store instances.")
    parser.add_argument("no", type=int,
                        help="Number of instances to be generated.")
    parser.add_argument("n", type=int,
                        help="Number of cities including the depot.")
    parser.add_argument("-sym", action="store_true", default=False,
                        help="Whether the graph is symmetric or not. Default is False.")
    parser.add_argument("-l", type=int, default=0,
                        help="Lower bound for entries in the cost matrix. Default is 0.")
    parser.add_argument("-u", type=int, default=10,
                        help="Upper bound for entries in the cost matrix. Default is 10.")
    parser.add_argument("-et", type=int, default=10,
                        help="Upper bound for the earliest start time. Default is 20.")
    parser.add_argument("-lt", type=int, default=40,
                        help="Upper bound for the latest start time. Default is 40.")
    parser.add_argument("-txt", action="store_true", default=False,
                        help="Whether to save the file in txt format as well. Default is False.")
    parser.set_defaults(pre=True)
    args = parser.parse_args()

    try:
        os.mkdir(args.dir)
    except OSError:
        print("Creation of the directory failed")
    else:
        print("Successfully created the directory")

    if args.txt:
        try:
            os.mkdir(f"{args.dir}_txt")
        except OSError:
            print("Creation of the directory failed")
        else:
            print("Successfully created the directory")

    complete(args.dir, args.no, args.sym, args.n, args.l, args.u, args.et, args.lt, args.txt)



