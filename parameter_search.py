import argparse
from itertools import product
from experiment import *
import pandas as pd
from params_helpers import *


# Search parameters for ILP formulation
def search_ilp(insdir, out, lp1, up1, lp2, up2):
    try:
        os.mkdir(out)
    except OSError:
        print("Creation of the directory failed or directory already exists")

    for instance_name in os.listdir(insdir):
        instance_name = instance_name[:-4]

        print("-------------------File ", instance_name, "is processed------------------")
        try:
            os.mkdir(f"{out}\{instance_name}")
        except OSError:
            print("Creation of the directory failed")

        ins = Instance_ilp(instance_name, insdir)
        earliest, latest = ins.earliest, ins.latest
        C = ins.C

        # Find optimal route by brute force
        opt_route, opt_cost = brute_force_tsptw(C, latest, earliest)

        # Convert route to array
        x = route_to_array_ilp(opt_route, ins)

        index = 0

        # Parameter search
        for p1, p2 in product(np.linspace(lp1, up1, num=(up1 - lp1 + 1)), np.linspace(lp2, up2, num=(up2 - lp2 + 1))):
            # Get ising formulation
            h, J, of = ising_ilp(ins, p1, p2)

            # Run simulated annealing
            sampleset = anneal(h, J, BETA_RANGE, NUM_READS, NUM_SWEEPS, BETA_SCHEDULE_TYPE)

            # Evaluate each sample
            Results = evaluate_sampleset_ilp(ins, sampleset)

            # Prepare pandas dataframe for the results
            Results.append((False, False, 100, 0, 0, 0, 0))

            data = pd.DataFrame(Results)
            data.columns = ['valid', 'windows', 'cost', 'energy', 'A', 'B', 'E']
            data['A'] = p1
            data['B'] = 1
            data['E'] = p2

            # Energy of the optimal route
            energy = dimod.ising_energy(binary_to_spin(x), h, J)
            data.loc[len(data)] = [True, True, opt_cost, energy, -1, 0, 0]

            # Store the results
            data.to_pickle(f"{out}\{instance_name}\{instance_name}_{index}")
            index += 1


# Search parameters for edge-based formulation
def search_edge(insdir, out, lp1, up1, lp2, up2):

    for instance_name in os.listdir(insdir):
        instance_name = instance_name[:-4]

        print("-------------------File ", instance_name, "is processed------------------")
        try:
            os.mkdir(f"{out}\{instance_name}")
        except OSError:
            print(f"Creation of the directory {out}\{instance_name} failed")
            quit()

        ins = Instance_edge(instance_name, insdir)
        C = ins.C
        earliest, latest = ins.earliest, ins.latest

        # Find optimal route by brute force
        opt_route, opt_cost = brute_force_tsptw(C, latest, earliest)
        x = route_to_array_edge(opt_route, ins)

        index = 0
        # Parameter search
        for p1, p2 in product(np.linspace(lp1, up1, num=(up1 - lp1 + 1)), np.linspace(lp2, up2, num=(up2 - lp2 + 1))):

            # Get ising formulation
            h, J, of = ising_edge(ins, p1, p2)

            # Run simulated annealing
            sampleset = anneal(h, J, BETA_RANGE, NUM_READS, NUM_SWEEPS, BETA_SCHEDULE_TYPE)

            # Evaluate each sample
            Results = evaluate_sampleset_edge(ins, sampleset)

            # Prepare pandas dataframe for the results
            Results.append((False, False, 100, 0, 0, 0, 0))

            data = pd.DataFrame(Results)
            data.columns = ['valid', 'windows', 'cost', 'energy', 'A', 'B', 'E']
            data['A'] = p1
            data['B'] = 1
            data['E'] = p2

            # Energy of the optimal route
            energy = dimod.ising_energy(binary_to_spin(x), h, J)
            data.loc[len(data)] = [True, True, opt_cost, energy, -1, 0, 0]

            # Store the results

            data.to_pickle(f"{out}/{instance_name}/{instance_name}_{index}")
            index += 1


def probs(out):
    for folder in os.listdir(out):
        if "summary" not in folder:
            columns = ['A', 'E', 'p', 'mp']
            results = pd.DataFrame(columns=columns)

            for f in os.listdir(f"{out}/{folder}"):

                data = pd.read_pickle(f"{out}/{folder}/{f}")
                data.reset_index(drop=True, inplace=True)
                data = data[data.cost!=100]

                p = probability(data)
                mp = min_prob(data)

                results = results.append(dict(zip(results.columns, [data.A[0], data.E[0],p,mp])), ignore_index=True)
                write_npz_file(f"{out}/summary", folder, results=results)

def find_params(out):
    dir = f"{out}/summary"
    columns = ['A', 'E', 'p', 'mp', 'instance', 'n']
    df = pd.DataFrame(columns=columns)
    for instance in os.listdir(dir):
        data = load_npz_file(dir, instance[:-4], "results")

        n=0
        if "_0_3" in instance:
            n = 3
        elif "_0_4" in instance:
            n = 4
        elif "_0_5" in instance:
            n = 5

        ns = [[f"{n}"]] * len(data)
        ins = [[f"{instance[:-4]}"]] * len(data)

        data = load_npz_file(dir, instance[:-4], "results")
        data = np.append(data, ins, axis=1)
        data = np.append(data, ns, axis=1)

        df = df.append(pd.DataFrame(data, columns=df.columns), ignore_index=True)

        dict = {}

    for folder in os.listdir(dir):
        soln = (df[df.instance == folder[:-4]].sort_values('mp', ascending=False).iloc[0])
        dict[folder[:-4]] = (soln.A, soln.E)
    write_npz_file(out, "params", dict = dict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("ins", type=str,
                        help="Name of the folder which stores the instances.")
    parser.add_argument("out", type=str, default="instances",
                        help="Name of the folder to store the results. Default is results.")
    parser.add_argument("model", type=str, choices=["e", "i"],
                        help="Formulation to be used. Type e for edge and i for ILP.")
    parser.add_argument("-lp1", type=int, default=1,
                        help="Lower bound for the parameter p1. Default is 1.")
    parser.add_argument("-up1", type=int, default=5,
                        help="Upper bound for the parameter p1. Default is 5.")
    parser.add_argument("-lp2", type=int, default=1,
                        help="Lower bound for the parameter p2. Default is 1.")
    parser.add_argument("-up2", type=int, default=100,
                        help="Upper bound for the parameter p2. Default is 100.")
    parser.add_argument("-reads", type=int, default=100,
                        help="Number of reads for simulated annealing. Default is 100.")
    parser.add_argument("-steps", type=int, default=1000,
                        help="Number of steps for simulated annealing. Default is 1000.")
    parser.add_argument("-lb", type=int, default=5,
                        help="Lower bound for beta range. Default is 5.")
    parser.add_argument("-ub", type=int, default=100,
                        help="Upper bound for beta range. Default is 100.")
    parser.add_argument("-schedule", type=str, default="geometric",
                        help="Schedule type. Default is geometric")
    parser.set_defaults(pre=True)
    args = parser.parse_args()

    try:
        os.mkdir(args.out)
    except OSError:
        if os.path.isdir(args.out):
            print(f"Directory {args.out} already exists. Try removing the directory.")
        else:
            print(f"Creation of the directory {args.out} failed")
        quit()

    try:
        os.mkdir(f"{args.out}/summary")
    except OSError:
        print(f"Creation of the directory {args.out}/summary failed")
        quit()


    BETA_RANGE = (args.lb, args.ub)
    NUM_READS = args.reads
    NUM_SWEEPS = args.steps
    BETA_SCHEDULE_TYPE = args.schedule


    if args.model == "e":
        search_edge(args.ins, args.out, args.lp1, args.up1, args.lp2, args.up2)
    elif args.model == "i":
        search_ilp(args.ins, args.out, args.lp1, args.up1, args.lp2, args.up2)

    probs(args.out)
    find_params(args.out)
