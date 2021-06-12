from itertools import permutations
from instance import *

class Instance_ilp(Instance):

    def __init__(self, name, dir):
        super().__init__(name, dir)
        self.invalid = self.invalid_edges(self.C,self.earliest,self.latest,self.n)

    def invalid_edges(self,C, earliest,latest,n):

        invalid = []

        for i, j in permutations(range(n), 2):
            if j != 0 and earliest[i] + C[i][j] > latest[j]:
                invalid.append([i, j])

        return invalid
