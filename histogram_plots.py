import argparse
import os
from collections import defaultdict
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
import matplotlib.patches as mpatches
import re
def new_data_gen(dir):

    energyo, energyf, energyn = defaultdict(lambda: 0), defaultdict(lambda: 0), defaultdict(lambda: 0)

    for instance_dir in os.listdir(dir):
        if "energy" not in instance_dir and "sample" not in instance_dir:
            data = pd.read_pickle(f"{dir}/{instance_dir}")

            k = [m.start() for m in re.finditer(r"_", instance_dir)][3]
            instance_dir = instance_dir[:k]

            df = pd.DataFrame(data)
            df.loc[((df.feasible == True) & (df.optimal == False)), 'optimal'] = "Feasible"
            df.loc[df.feasible == False, 'optimal'] = "Not feasible"
            df.loc[df.optimal == True, 'optimal'] = "Optimal"

            energyo[instance_dir] = df[df.optimal == "Optimal"].energy
            energyf[instance_dir] = df[df.optimal == "Feasible"].energy
            energyn[instance_dir] = df[df.optimal == "Not feasible"].energy


    return energyo, energyf, energyn


def plt_hist(axis, data, hatch, label, color,ecolor, e):
    binwidth = 0.1
    counts, edges = np.histogram(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
    edges = np.repeat(edges, 2)
    hist = np.hstack((0, np.repeat(counts, 2), 0))

    outline, = axis.plot(edges, hist, linewidth=1.3, color=color)
    axis.fill_between(edges, hist, 0,
                      edgecolor=ecolor, hatch=hatch, label=label,
                      facecolor=color,alpha=0.4)  ## < removes facecolor
    axis.set_ylim(0, None, auto=True)


def plt_hists(ins, model, energyo, energyf, energynf,dir):

    e = np.load(f"{dir}/{ins}_energy", allow_pickle=True)

    k = [m.start() for m in re.finditer(r"_", ins)][3]
    instance = ins[:k]
    h1 = '---'
    d1 = energyo[instance]
    lab1 = 'Optimal'
    h2 = '\\\\\\'
    d2 = energyf[instance]
    lab2 = 'Feasible'
    h3 = '///'
    d3 = energynf[instance]
    lab3 = 'Not feasible'
    nf_color = '#009966'
    nfe_color = '#013624'
    o_color = '#e1a84e'
    oe_color = '#a16200'
    f_color = '#9c59dd'
    fe_color = '#430085'

    rc("text", usetex=True)
    rc("font", family="serif")

    fig, ax = plt.subplots(1)
    fig.set_figheight(5)
    fig.set_figwidth(8)
    plt.xlabel("Energy")
    plt.ylabel("Count")

    if d1.size> 0:
        plt_hist(ax,d1,h1,lab1,o_color,oe_color,e)
    if d2.size> 0:
        plt_hist(ax,d2,h2,lab2,f_color,fe_color,e)
    if d3.size > 0:
        plt_hist(ax,d3,h3,lab3,nf_color,nfe_color,e)

    plt.axvline(e, 0, 1, color='red')
    plt.tick_params(left=True, bottom=True)


    o = mpatches.Patch(facecolor=o_color, alpha=0.5, label='Optimal', hatch=h1)
    f = mpatches.Patch(facecolor=f_color, alpha=0.3, label='Feasible', hatch=h2)
    nf = mpatches.Patch(facecolor=nf_color, alpha=0.3, label='Not feasible', hatch=h3)

    df = [d1, d2, d3]
    legs = [o, f, nf]
    arr = [legs[i] for i in range(len(df)) if df[i].size > 0]

    plt.legend(handles=arr)

    for item in ([ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels() +
                 ax.get_legend().get_texts()):
        item.set_fontsize(22)

    plt.savefig(f"plots/hist_{ins}_{model}.pdf", bbox_inches='tight')

    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("model", type=str, choices=["e", "i"],
                        help="Formulation to be used. Type e for edge and i for ILP.")
    parser.add_argument("-res", type=str, default="",
                        help="Name of the folder which stores the experiment results")
    parser.add_argument("-ins", type=str,
                        help="Name of the instance. If not entered, then the experiment is run for all instances in the folder.")

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

    energyo, energyf, energyn = new_data_gen(args.res)
    if args.ins:
        plt_hists(args.ins, args.model, energyo, energyf, energyn,args.res)
    else:
        for instance_name in os.listdir(args.res):
            if "energy" not in instance_name and "sample" not in instance_name:
                plt_hists(instance_name, args.model, energyo, energyf, energyn,args.res)