import networkx as nx
from utils import *


class Instance:

    def __init__(self, name, dir):
        self.name = name  # filename
        self.dir = dir  # directory

        self.C = load_npz_file(dir, name, 'weight')  # Cost matrix
        self.n = len(self.C)  # Number of cities

        self.w = load_npz_file(dir, name, 'window')  # Time windows
        self.earliest = self.w.T[0]  # Ready times
        self.latest = self.w.T[1]  # Due times

        self.h = {} #ising model h
        self.J = {} #ising model J
        self.o = 0  # ising model offset

        self.update_earliest()

    def update_earliest(self):
        for i in range(1,self.n):
            if self.C[0][i]> self.earliest[i]:
                self.earliest[i] = self.C[0][i]

    # Store h, J
    def store_hJ(self, index):
        write_npz_file("data_p", f"{self.name}_{index}_data", h=self.h, J = self.J, o = self.o)

    #Save hamiltonian
    def write_hamiltonian(self, model, dir, index=0):

        filename = f"{dir}/{self.name}/{index}" if index != 0 else f"{dir}/{self.name}"
        f = open(filename, "w")

        self.h, self.J, self.o = Instance.convert_and_scale(model)

        for key, value in self.h.items():
            print(key, key, value, file=f)
        for key, value in self.J.items():
            print(key[0], key[1], value, file=f)

    @staticmethod
    def convert_and_scale(model):
        # Convert to spin and scale
        h, J, offset = model.to_ising()
        return Instance.scaleH(h, J,offset)

    @staticmethod
    def scaleH(h,J,offset):
        minh, maxh, minj, maxj = -2, 2, -1, 1
        scaleFactor = max(max(max(h.values())/maxh,0), max(min(h.values())/minh,0), max(max(J.values())/maxj, 0), max(min(J.values())/minj,0))
        scaled_h, scaled_J = {},{}
        for key in h.keys():
            scaled_h[key] = h[key] / scaleFactor
        for key in J.keys():
            scaled_J[key] = J[key] / scaleFactor
        return scaled_h, scaled_J, offset/scaleFactor
