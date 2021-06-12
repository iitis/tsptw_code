from instance import *


class Instance_edge(Instance):

    def __init__(self, name, dir):
        super().__init__(name,dir)
        self.W_arr = []  # Upper bound for waiting times
        self.L = 0  # Upper bound for slack vars.
        self.D_arr = []  # Upper bound for slack var.

        self.G = nx.from_numpy_matrix(self.C, create_using=nx.DiGraph)  # Graph
        self.G_0 = self.G.subgraph(range(1, self.n))  # Subgraph without node 0
        self.elist = list(self.G_0.edges())  # Edge list without node 0

        self.Map = self.generate_map()  # Map for the variables
        self.Map_w, self.Map_u, self.Map_s = self.generate_map_w()  # Map for the remaining variables
        self.invalid = self.invalid_edges(self.C,self.earliest,self.latest,self.Map,self.n,self.elist)

    def invalid_edges(self,C, earliest,latest,Map,n,elist):
        zeros = []

        for u, v in elist:
            if earliest[u] + C[u][v] > latest[v]:
                if u != 0 and v != 0:
                    z = [Map.get(((u, v), t)) for t in range(2, n)]
                elif u == 0:
                    z = Map.get(((u, v), 1))
                elif v == 0:
                    z = Map.get(((u, v), n))
                zeros += z

        return sorted(zeros)

    # Create map for the waiting times and slack variables
    def generate_map_w(self):
        n, C = self.n, self.C
        earliest, latest = self.earliest, self.latest
        Map_w, Map_s, Map_u= {}, {}, {}

        k = (n - 1) ** 3 - 2 * (n - 1) ** 2 + 3 * (n - 1)

        self.D_arr = [find_bound(C, t, latest) for t in range(1, n)]
        for t in range(1, n):
            for d in range(0, self.D_arr[t - 1]):
                Map_s[(d, t)] = k
                k = k + 1

        self.W_arr = [find_bound(C, t, earliest) for t in range(1, n)]
        for t in range(1, n):
            for w in range(0, self.W_arr[t - 1]):
                Map_w[(w, t)] = k
                k = k + 1

        self.L = max(math.floor(math.log(max(1, max(latest - earliest)), 2)) + 1, 1)
        for t in range(1, n):
            for d in range(self.L):
                Map_u[(d, t)] = k
                k = k + 1

        return Map_w, Map_u, Map_s

    #Generate map for x_uv^i variables
    def generate_map(self):
        n = self.n
        Map = {}
        counter = 0
        for v in self.G.neighbors(0):
            Map[(0, v), 1] = counter
            counter = counter + 1
        for t in range(2, n):
            for e in self.elist:
                Map[e, t] = counter
                counter = counter + 1
        for u in self.G.successors(0):
            Map[(u, 0), n] = counter
            counter = counter + 1
        return Map

    # Save data
    def store_data(self,dir):
        write_npz_file(dir, f"{self.name}_data", Map=self.Map, Map_w=self.Map_w, Map_u=self.Map_u,
                       Map_s=self.Map_s, E=self.earliest,
                       L=self.latest, h=self.h, J=self.J)