from numpy import *


class Element:
    def __init__(self, element_id, node_1, node_2, node_3, node_4, k):
        self.element_id = element_id
        self.nodes = [node_1, node_2, node_3, node_4]
        self.heated_areas = [0, 0, 0, 0]
        self.k = k
        self.element_matrix_h = zeros([4, 4])
        self.element_matrix_c = zeros([4, 4])
        self.element_vector_p = zeros(4)

    def show_element_params(self):
        print("Element ID: {}, Elements' nodes: {}, {}, {}, {}  k={}".format(self.element_id, self.nodes[0].node_id,
                                                                             self.nodes[1].node_id,
                                                                             self.nodes[2].node_id,
                                                                             self.nodes[3].node_id, self.k))

    # def show_element_matrices(self):
    #     print("Matrix H:\n{} {} {} {}\n{} {} {} {}\n{} {} {} {}\n{} {} {} {}\n".format(self.H_matrix[0]))
