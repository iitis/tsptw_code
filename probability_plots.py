import argparse
import os
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import re


def calc_prob(model,dir,cities):

    min_found, fprob, oprob = defaultdict(lambda: 0), defaultdict(lambda: 0), defaultdict(lambda: 0)

    for instance_dir in os.listdir(dir):
        if "energy" not in instance_dir and "sample" not in instance_dir and "qa" in instance_dir:
            j = [m.start() for m in re.finditer(r"_", instance_dir)][2]
            k = [m.start() for m in re.finditer(r"_", instance_dir)][3]
            n = instance_dir[j-1:k-2]

            if n==cities:

                no = instance_dir[j+1:k]
                data = pd.read_pickle(f"{dir}/{instance_dir}")

                oprob[no] = len(data[data.optimal == True]) / len(data)
                fprob[no] = len(data[data.feasible == True]) / len(data)
                min_found[no] = True if (data['optimal'].iloc[0] == True) else False

    df = pd.DataFrame()

    df["Instances"] = list(oprob.keys()) * 2
    df["Probability"] = list(oprob.values()) + list(fprob.values())
    df["Optimal"] = ["Optimal"] * len(oprob.keys()) + ["Feasible"] * len(oprob.keys())


    return df


def plot_prob(df,cities,model):
    o = df[df.Optimal == "Optimal"]
    f = df[df.Optimal == "Feasible"]
    x = range(len(o))
    rc("text", usetex=True)
    rc("font", family="serif")

    plt.figure(figsize=(8, 6), dpi=80)

    ax = plt.bar(x, f["Probability"], color='#9c59dd')
    for i, bar in enumerate(ax.patches):
        hatch = '\\\\'
        bar.set_hatch(hatch)
        bar.set_edgecolor('#5d02b5')
    ax = plt.bar(x, o["Probability"], color='#e1a84e')
    for i, bar in enumerate(ax.patches):
        hatch = '//'
        bar.set_hatch(hatch)
        bar.set_edgecolor('#a16302')

    plt.xlabel('Instances', fontsize=30)
    plt.ylabel('Probability', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)

    plt.xticks(range(0, 10))
    plt.tick_params(left=True, bottom=True)
    plt.savefig(f"plots/probs_{cities}_{model}.pdf", bbox_inches='tight')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("model", type=str, choices=["e", "i"],
                        help="Formulation to be used. Type e for edge and i for ILP.")
    parser.add_argument("cities", type=str,
                        help="Number of cities.")
    parser.add_argument("-res", type=str, default="",
                        help="Name of the folder which stores the experiment results")


    args = parser.parse_args()

    if not os.path.exists("plots"):
        try:
            os.mkdir("plots")
        except OSError:
            print("Creation of the directory plots failed.")

    if args.res=="":
        if args.model =="e":
            args.res = "data/dwave_results_edge"
        elif args.model =="i":
            args.res = "data/dwave_results_ilp"

    df = calc_prob(args.model,args.res,args.cities)
    plot_prob(df,args.cities,args.model)
