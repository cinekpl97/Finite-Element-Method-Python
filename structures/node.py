from numpy import *


class Node:
    def __init__(self, x, y, node_id, temperature):
        self.node_id = node_id
        self.x = x
        self.y = y
        self.temperature = temperature

    def show_nodes_params(self):
        print("Node ID: {}, Coordinates: {}, {},\t Temperature: {}".format(self.node_id, self.x, self.y, self.temperature))
