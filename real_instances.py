import argparse
import pickle
from math import log,floor
import pandas as pd
from instance_ilp import *
from instance_node import *
from utils import *
import os
import matplotlib.pyplot as plt
from matplotlib import rc


def generate_npz(dir_orig,dir,ext):
    for instance_name in os.listdir(dir_orig):
        lines = tuple(open(f"{dir_orig}/{instance_name}", 'r'))
        n = int(lines[0])
        M = []
        for i in range(int(n)):
            M.append(list(map(int, lines[i + 1].split())))

        N = []
        for i in range(int(n)):
            N.append(list(map(int, lines[n + i + 1].split())))
        N = np.array(N)
        k = (len(ext)+1)*-1
        instance_name = instance_name[:k]
        file_name = f"{dir}/{instance_name}"
        np.savez(file_name, n=n, weight=M, window=N)


def get_vars_instance(instance_name,dir,n,model):

    A,B,E = 1,1,1
    if model =="Edge":
        ins = Instance_edge(instance_name, dir)
        n = ins.n - 1
        edge_var = n**3 - 2*n*n + 3*n - len(ins.invalid)
        wait_var = sum(ins.W_arr)
        d_var = sum(ins.D_arr)
        L_var = ins.L*(n)
        varmax = edge_var + wait_var + d_var + L_var
        var = edge_var + wait_var
        slack = d_var + L_var

    elif model == "ILP":
        ins = Instance_ilp(instance_name,dir)
        n = ins.n
        earliest, latest = ins.earliest, ins.latest
        C =  ins.C
        edge_var = 0
        for i, j in permutations(range(0, n), 2):
            if ([i, j] not in ins.invalid):
                edge_var+=1
        service_var, wait_var,eq24var, eq25var,eq27var =0,0,0,0,0

        for i in range(1,ins.n):
            service_var += floor(log(ins.latest[i]-ins.earliest[i] ,2) )+1
            if ins.earliest[i] - ins.C[0][i] > 0:
                wait_var += floor(log(earliest[i] - C[0][i], 2) )+1
            eq24var += floor(log(latest[i],2))+1
            eq25var += floor(log(latest[i]-C[0][i] , 2) )+1
        for i,j in permutations(range(1,ins.n),2):
            if [i,j] not in ins.invalid:
                eq27var+= floor(log(-earliest[i] + latest[j] + latest[i] - C[0][j]  ,2))+1

        varmax = edge_var + service_var + wait_var + eq24var + eq25var + 2*eq27var
        var = edge_var + service_var + wait_var
        slack = eq24var + eq25var + 2*eq27var
    else:
        ins = Instance_node(instance_name, dir)
        n = ins.n - 1
        node_var = n ** 2
        wait_var = sum(ins.W_arr)
        d_var = sum(ins.D_arr)
        L_var = ins.L * (n)
        varmax = node_var + wait_var + d_var + L_var
        var = node_var + wait_var
        slack = d_var + L_var
    return var, slack, varmax

def num_var(dir, out):

    for instance_name in os.listdir(dir):

        instance_name = instance_name[:-4]

        C = load_npz_file(dir, instance_name, 'weight')
        n = len(C)
        qubo_vars,qubo_slack,qubo_all = get_vars_instance( instance_name, dir, n, "Edge")
        hobo_vars,hobo_slack,hobo_all = get_vars_instance( instance_name, dir, n, "Node")
        ilp_vars,ilp_slack,ilp_all = get_vars_instance( instance_name, dir, n, "ILP")

        write_npz_file(out, f"{instance_name}", n=n, qubo_vars = qubo_vars, hobo_vars = hobo_vars,
                       ilp_vars = ilp_vars, qubo_slack = qubo_slack, hobo_slack = hobo_slack,ilp_slack = ilp_slack,
                       qubo_all=qubo_all, hobo_all=hobo_all, ilp_all=ilp_all )


def analysis(dir,plot):
    df = pd.DataFrame([],
                      columns=["instance_name", "size", "qubo_all", "qubo_slack", "hobo_all", "hobo_slack", "ilp_all",
                               "ilp_slack"])

    for instance_name in os.listdir(dir):
        n= int(load_npz_file(dir, instance_name[:-4], "n"))
        qubo_vars = int(load_npz_file(dir, instance_name[:-4], "qubo_vars"))
        hobo_vars = int(load_npz_file(dir, instance_name[:-4], "hobo_vars"))
        ilp_vars = int(load_npz_file(dir, instance_name[:-4], "ilp_vars"))

        qubo_all = int(load_npz_file(dir, instance_name[:-4], "qubo_all"))
        hobo_all = int(load_npz_file(dir, instance_name[:-4], "hobo_all"))
        ilp_all = int(load_npz_file(dir, instance_name[:-4], "ilp_all"))

        qubo_slack = int(load_npz_file(dir, instance_name[:-4], "qubo_slack"))
        hobo_slack = int(load_npz_file(dir, instance_name[:-4], "hobo_slack"))
        ilp_slack = int(load_npz_file(dir, instance_name[:-4], "ilp_slack"))

        qubo_slack = np.round(qubo_slack / qubo_all * 100, 1)
        hobo_slack = np.round(hobo_slack / hobo_all * 100, 1)
        ilp_slack = np.round(ilp_slack / ilp_all * 100, 1)

        df = pd.DataFrame(np.insert(df.values, 0,
                                    values=[instance_name[:-4], n, qubo_all, qubo_slack, hobo_all, hobo_slack,
                                            ilp_all, ilp_slack], axis=0))
        df = pd.DataFrame(df.values,
                          columns=["instance_name", "size", "qubo_all", "qubo_slack", "hobo_all", "hobo_slack",
                                   "ilp_all", "ilp_slack"])

    with open(f"{dir}/analysis", 'wb') as handle:
        pickle.dump(df, handle)

    if plot:
        rc("text", usetex=True)
        rc("font", family="serif")
        plt.rcParams.update({'font.size': 25})

        plt.figure(figsize=(12, 6))
        plt.scatter(df["size"], (df.qubo_all), marker='x', s=100)
        plt.scatter(df["size"], (df.hobo_all), marker='.', s=200)
        plt.scatter(df["size"], (df.ilp_all), marker='+', s=200)

        plt.ylabel("Logical variables")
        plt.xlabel("Number of cities")
        plt.legend(["Edge-based", "ILP", "Node-based"])
        plt.yscale('log')

        try:
            os.mkdir("plots")
        except OSError:
            print("Creation of the directory failed or directory already exists")

        plt.savefig(f"plots/real_instances_vars.pdf", bbox_inches='tight')

if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("insf", type=str,
                        help="Location of the folder which stores the instances.")
    parser.add_argument("ins", type=str,
                        help="Name of the folder which stores the instances.")
    parser.add_argument("ext", type=str,
                        help="Extension of the instances.")
    parser.add_argument("out", type=str, default="vars",
                        help="Name of the folder to store the results. Default is vars.")
    parser.add_argument("-plot", action="store_true", default=False,
                        help="Whether the results will be plotted or not.")

    args = parser.parse_args()

    try:
        os.mkdir(args.out)
    except OSError:
        print("Creation of the directory failed or directory already exists")

    try:
        os.mkdir(f"{args.out}/{args.ins}_npz")
    except OSError:
        print("Creation of the directory failed or directory already exists")

    try:
        os.mkdir(f"{args.out}/{args.ins}_count")
    except OSError:
        print("Creation of the directory failed or directory already exists")

    generate_npz(f"{args.insf}/{args.ins}",f"{args.out}/{args.ins}_npz",args.ext)
    num_var(f"{args.out}/{args.ins}_npz",f"{args.out}/{args.ins}_count")
    analysis(f"{args.out}/{args.ins}_count",args.plot)



