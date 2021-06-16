import argparse
import os

from utils import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def penalty_search_plot(out,instance,mode,model):

    results = pd.read_pickle(f"{out}/summary/{instance}")

    results = results.astype({"E": int, "A": int})

    plot_data = results.pivot_table(index='A', columns='E', values="mp")
    sns.set(font_scale=1.5)
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")
    f = sns.heatmap(plot_data, cmap="rocket_r")

    f.axhline(y=0, color='k',
              linewidth=1)

    f.axhline(y=4.99, color='k',
              linewidth=1)

    f.axvline(x=0, color='k',
              linewidth=1)

    f.axvline(x=99.5, color='k',
              linewidth=1)


    plt.xticks(rotation='horizontal')
    for ind, label in enumerate(f.get_xticklabels()):
        if ind % 3 == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)
    plt.tick_params(left=True, bottom=True)

    p1 = "$p_1$"
    p2 = "$p_2$"

    plt.xlabel(p2)
    plt.ylabel(p1)

    f.figure.set_figwidth(5)
    f.figure.set_figheight(3)
    f.figure.savefig(f"plots/params_{instance}_{model}_{mode}.pdf", bbox_inches='tight', dpi=f.figure.dpi)

    plt.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("results", type=str, default="",
                        help="Name of the folder which stores the parameter search results")
    parser.add_argument("model", type=str, choices=["e", "i"],
                        help="Formulation to be used. Type e for edge and i for ILP.")
    parser.add_argument("mode", type=str,
                        help="Type p to plot the probabilities and mp to plot the probabilities where probability is set to 0 if the lowest energy sample does not encode an optimal route")
    parser.add_argument("-ins", type=str,
                        help="Name of the instance.")

    args = parser.parse_args()

    if not os.path.exists("plots"):
        try:
            os.mkdir("plots")
        except OSError:
            print("Creation of the directory plots failed.")

    if args.ins:
        penalty_search_plot(args.results, args.ins, args.mode, args.model)
    else:
        for instance_name in os.listdir(args.results):
            if "summary" not in instance_name and "params" not in instance_name:
                penalty_search_plot(args.results, instance_name, args.mode, args.model)



