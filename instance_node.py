from itertools import product
from instance_edge import *


class Instance_node(Instance):

    def __init__(self, name, dir):
        super().__init__(name, dir)


        self.W_arr = []  # Upper bound for waiting times
        self.L = 0  # Upper bound for slack vars.
        self.D_arr = []  # Upper bound for slack var.
        self.Map = self.generate_map()  # Map for the variables
        self.Map_w, self.Map_u, self.Map_s = Instance_edge.generate_map_w(self)  # Map for the remaining variables

    #Generate map for x_uv^i variables
    def generate_map(self):
        n = self.n
        Map = {}
        counter = 0
        for v,t in product(range(1, n),range(1,n)):
            Map[(v, t)] = counter
            counter = counter + 1
        return Map

