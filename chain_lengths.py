from minorminer import find_embedding
from dwave.cloud import Client
from dwave.embedding import embed_bqm
from dwave.embedding.utils import edgelist_to_adjacency
import pandas as pd
from generate_ham_edge import *
from instance_edge import *


dir = "data/instances"

client = Client.from_config()
solver = client.get_solver(qpu=True)
target_edgelist = solver.edges

s,t,c,n = [], [],[],[]
for instance_name in os.listdir(dir):
    instance_name = instance_name[:-4]

    ins = Instance_edge(instance_name, dir)  # Generate a new instance object

    P1,P2 = 1,1
    print(instance_name)

    bqm = generate_bqm_edge(ins,P1,P2)

    source_edgelist = list(bqm.quadratic) + [(v, v) for v in bqm.linear]

    embedding = find_embedding(source_edgelist, target_edgelist)
    target_adjacency = edgelist_to_adjacency(target_edgelist)
    bqm_embedded = embed_bqm(bqm, embedding, target_adjacency)

    num_source_vars = len(bqm.variables)
    target_vars = set(bqm_embedded.variables)
    num_target_variables = len(bqm_embedded.variables)
    max_chain_length = max(len(target_vars.intersection(chain))
                           for chain in embedding.values())
    s.append(num_source_vars)
    t.append(num_target_variables)
    c.append(max_chain_length)
    n.append(ins.n)


    df = pd.DataFrame()
    df['Logical variables'] = s
    df['Physical variables'] = t
    df['Chain length'] = c
    df['Number of cities'] = n
print(df)

