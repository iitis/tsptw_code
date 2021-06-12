import argparse
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import rc

def plot_chain_lengths(model,dir):

    logical, chain_length, physical, size = defaultdict(lambda: 0), defaultdict(lambda: 0), defaultdict(
        lambda: 0), defaultdict(lambda: 0)
    for instance_name in os.listdir(dir):
        if "sample" in instance_name:
            data = pd.read_pickle(f"{dir}/{instance_name}")
            l = list(data["info"]['embedding_context']['embedding'].values())
            chain_length[instance_name] = len(max(l, key=len))
            logical[instance_name] = data['num_variables']
            l = list(data["info"]['embedding_context']['embedding'].values())
            physical[instance_name] = len(set.union(*[set(x) for x in l]))
            size[instance_name] = instance_name[11]
    df = pd.DataFrame()
    df['Logical variables'] = list(logical.values())
    df['Physical variables'] = list(physical.values())
    df['Chain length'] = list(chain_length.values())
    df['Number of cities'] = list(size.values())
    df['Model'] = model
    return df



if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("-resi", type=str,default = "data/dwave_results_ilp",
                        help="Name of the folder which stores the experiment results for ILP formulation.")
    parser.add_argument("-rese", type=str,default = "data/dwave_results_edge",
                        help="Name of the folder which stores the experiment results for edge-based formulation.")


    args = parser.parse_args()

    try:
        os.mkdir("plots")
    except OSError:
        print("Creation of the directory failed or directory already exists")


    df1 = plot_chain_lengths("Edge-based",args.rese)
    df2 = plot_chain_lengths("ILP", args.resi )
    df = df1.append(df2)

    sns.set(font_scale=4, style='white')


    rc("text", usetex=True)
    rc("font", family="serif")

    p = sns.color_palette("tab10")[:3]
    g = sns.relplot(x="Logical variables", y="Physical variables", hue="Number of cities",
                    style="Number of cities", legend=False,
                    s=1000, height=10, aspect=1.3, data=df,
                    col="Model", palette=p)



    c3 = mlines.Line2D([], [], color=p[0], marker='o', linestyle='None',
                       markersize=20, label='3 cities')
    c4 = mlines.Line2D([], [], color=p[1], marker='X', linestyle='None',
                       markersize=20, label='4 cities')
    c5 = mlines.Line2D([], [], color=p[2], marker='s', linestyle='None',
                       markersize=20, label='5 cities')

    plt.legend(handles=[c3, c4, c5], bbox_to_anchor=(1.05, 1))

    plt.xlabel("Logical variables")
    plt.ylabel("Physical variables")

    g.axes.flat[0].set_title('Edge-based')
    g.axes.flat[1].set_title('ILP')
    g.axes.flat[0].tick_params(left=True, bottom=True)
    g.axes.flat[1].tick_params(left=True, bottom=True)

    g.savefig(f"plots/variables.pdf", bbox_inches='tight')



