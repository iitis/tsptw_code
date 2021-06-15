import argparse
import pickle
import neal
from dwave.system import EmbeddingComposite, DWaveSampler, LeapHybridSampler
from dimod import *
import pandas as pd
from edge_helpers import *
from generate_ham_edge import *
from generate_ham_ilp import *
from ilp_helpers import *



def anneal(h,J,BETA_RANGE,NUM_READS,NUM_SWEEPS,BETA_SCHEDULE_TYPE):
    s = neal.SimulatedAnnealingSampler()

    sampleset = s.sample_ising(h, J, beta_range=BETA_RANGE, num_sweeps=NUM_SWEEPS, num_reads=NUM_READS,
                               beta_schedule_type=BETA_SCHEDULE_TYPE)

    return sampleset

def float_range(mini,maxi):
    # Define the function with default arguments
    def float_range_checker(arg):
        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < mini or f > maxi:
            raise argparse.ArgumentTypeError("must be in range [" + str(mini) + " .. " + str(maxi)+"]")
        return f

    return float_range_checker

def exp(ins, model, alg, p1, p2):
    if model == "e":
        h, J, of = ising_edge(ins, p1, p2)
    elif model == "i":
        h, J, of = ising_ilp(ins, p1, p2)

    if alg == 'sa':
        sampleset = anneal(h,J,BETA_RANGE,NUM_READS,NUM_SWEEPS,BETA_SCHEDULE_TYPE)
    elif alg == 'qa':
        sampler = EmbeddingComposite(DWaveSampler())
        sampleset = sampler.sample_ising(h, J, num_reads=NUM_READS, auto_scale='false', annealing_time=TIME, chain_strength=CHAIN)
    elif alg == 'hyb':
        sampler = LeapHybridSampler()
        sampleset = sampler.sample_ising(h,J)

    opt_route, opt_cost = brute_force_tsptw(ins.C, ins.latest, ins.earliest)
    x = route_to_array_edge(opt_route,ins) if model =="e" else route_to_array_ilp(opt_route, ins)
    energy = dimod.ising_energy(binary_to_spin(x), h, J)

    return sampleset, opt_cost, energy

def store_data(sampleset, ins, p1, p2, model, alg, out, opt_cost, energy):

    results = []
    for datum in sampleset.data(fields=['sample', 'energy']):
        x = dimod.sampleset.as_samples(datum.sample)[0][0]

        if model == "e":
            soln = array_to_soln_edge(x, ins.Map, ins.invalid)
            v, route, w, cost = check_solution_edge(soln, ins.n, ins.C, ins.earliest, ins.latest)
        elif model == "i":
            soln = array_to_soln_ilp(x, ins.n, ins.invalid)

            v, route, w, cost = check_solution_ilp(soln, ins.n, ins.C, ins.earliest, ins.latest)

        f = w and v
        o = cost == opt_cost and f
        results.append((x, route, cost, datum.energy, w, v, f, o))

    df = pd.DataFrame(results)
    df.columns = ['sample', 'route', 'cost', 'energy', 'windows', 'valid', 'feasible', 'optimal']

    if p1-floor(p1)==0:
        p1 = int(p1)
    if p2-floor(p2) == 0:
        p2 = int(p2)
    if alg == "qa" and CHAIN!=2:
        filename = f"{ins.name}_{p1}_{p2}_{alg}_{CHAIN}_{TIME}"
    elif alg == "qa" and CHAIN==2:
        filename = f"{ins.name}_{p1}_{p2}_{alg}_{TIME}"
    else:
        filename = f"{ins.name}_{p1}_{p2}_{alg}"
    df.to_pickle(f"{out}/{filename}")
    sdf = sampleset.to_serializable()
    with open(f"{out}/{filename}_sample", 'wb') as handle:
        pickle.dump(sdf, handle)
    with open(f"{out}/{filename}_energy", 'wb') as handle:
        pickle.dump(energy, handle)
    return df

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("instances", type=str,
                        help="Name of the folder which stores the instances.")
    parser.add_argument("out", type=str,
                        help="Name of the folder to store the results.")
    parser.add_argument("model", type=str, choices=["e", "i"],
                        help="Formulation to be used. Type e for edge and i for ILP.")
    parser.add_argument("alg", type=str, choices=["qa", "sa", "hyb"],
                        help="Type qa for quantum annealing, sa for simulated annealing and hyb for hybrid solvers")
    parser.add_argument("-ins", type=str,
                        help="Name of the instance. If not entered, then the experiment is run for all instances in the folder.")
    parser.add_argument("-dictf", type=str,
                        help="Location of the dictionary storing penalty parameters.")
    parser.add_argument("-dict", type=str,
                        help="Name of the dictionary storing penalty parameters.")
    parser.add_argument("-p1", type=float,
                        help="p1 penalty parameter. This is required if no dict parameter is given.")
    parser.add_argument("-p2", type=float,
                        help="p2 penalty parameter. This is required if no dict parameter is given.")
    parser.add_argument("-time", type=int, default=50,
                        help="Annealing time. Default is 50.")
    parser.add_argument("-reads", type=int, default=100,
                        help="Number of reads for simulated and quantum annealing. Default is 1000.")
    parser.add_argument("-chain", type=float_range(0,2), default=2,
                        help="Chain strength. Default is 2.")
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

    if not args.ins and not args.dict:
        parser.error('-dict is required when the instance name is not specified.')
    elif args.ins and (not args.p2 or not args.p1):
        parser.error('-p1 and -p2 are required')
    elif args.dict and (args.p1 or args.p2) :
        parser.error('You should enter either dict or p1 and p2.')

    BETA_RANGE = (args.lb, args.ub)
    NUM_READS = args.reads
    NUM_SWEEPS = args.steps
    BETA_SCHEDULE_TYPE = args.schedule
    TIME = args.time
    CHAIN  = args.chain

    try:
        os.mkdir(args.out)
    except OSError:
        if os.path.isdir(args.out):
            print(f"Directory {args.out} already exists. Results may be overwritten.")
        else:
            print(f"Creation of the directory {args.out} failed.")
            quit()


    if args.ins:
        ins = Instance_edge(args.ins,args.instances) if args.model == "e" else Instance_ilp(args.ins,args.instances)
        sampleset, opt_cost, energy = exp(ins, args.model, args.alg, args.p1, args.p2)
        store_data(sampleset, ins, args.p1, args.p2, args.model, args.alg, args.out, opt_cost, energy)
    else:
        for instance_name in os.listdir(args.instances):
            instance_name = instance_name[:-4]
            print("-----------------",instance_name, "is processed.-----------------")
            ins = Instance_edge(instance_name, args.instances) if args.model == "e" else Instance_ilp(instance_name,
                                                                                                 args.instances)

            fdict = pd.read_pickle(f"{args.dictf}/{args.dict}")
            p1 = fdict[instance_name][0]
            p2 = fdict[instance_name][1]

            sampleset, opt_cost, energy = exp(ins, args.model, args.alg, p1, p2)
            store_data(sampleset, ins, p1, p2, args.model, args.alg, args.out, opt_cost, energy)



