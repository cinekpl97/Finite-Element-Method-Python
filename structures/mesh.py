from numpy import *
from settings import *
import matplotlib.pyplot as plt
import seaborn as sb


class Mesh:
    def __init__(self):
        settings = get_settings()
        self.elements = []
        self.nodes = []
        self.nodes_count = settings["nH"] * settings["nL"]
        self.nH = settings["nH"]
        self.nL = settings["nL"]
        self.ambient_temperature = settings["ambient_temperature"]
        self.simulation_step_time = settings["simulation_step_time"]
        self.simulation_time = settings["simulation_time"]
        self.global_matrix_h = zeros((self.nodes_count, self.nodes_count))
        self.global_matrix_c = zeros((self.nodes_count, self.nodes_count))
        self.vector_p = zeros(self.nodes_count)

    def print(self):
        for nodes in self.nodes:
            nodes.show_nodes_params()

            print()

        for elements in self.elements:
            elements.show_element_params()

    def aggregation(self):
        for element in self.elements:
            for i in range(4):
                for j in range(4):
                    node_id_i = element.nodes[i].node_id - 1
                    node_id_j = element.nodes[j].node_id - 1

                    value_from_local_matrix_h = element.element_matrix_h.item((i, j))
                    global_matrix_h_value = self.global_matrix_h.item(
                        (node_id_i, node_id_j)) + value_from_local_matrix_h
                    self.global_matrix_h.itemset((node_id_i, node_id_j), global_matrix_h_value)

                    matrix_c_local_value = element.element_matrix_c.item((i, j))
                    matrix_c_global_value = self.global_matrix_c.item((node_id_i, node_id_j)) + matrix_c_local_value
                    self.global_matrix_c.itemset((node_id_i, node_id_j), matrix_c_global_value)

            # vector needs only one loop, cause its one-dimensional
            for i in range(4):
                node_id_i = element.nodes[i].node_id - 1

                vactor_p_local_value = element.element_vector_p.item(i)
                vector_p_global_value = self.vector_p.item(node_id_i) + vactor_p_local_value
                self.vector_p.itemset(node_id_i, vector_p_global_value)

    def heat_equation_solving(self):
        dt = self.simulation_step_time
        print("\n")
        print(self.global_matrix_c)

        left_equation = self.global_matrix_h + (self.global_matrix_c / dt)

        max_temp_list = []
        min_temp_list = []
        for step in range(0, self.simulation_time, dt):

            # Wektor t0 - temperatury kazdego z wezlow
            t0 = zeros(self.nodes_count)
            for i in range(self.nodes_count):
                t0.itemset(i, self.nodes[i].temperature)

            right_equation = -self.vector_p + (self.global_matrix_c / dt) @ t0

            t1 = linalg.solve(left_equation, right_equation)
            if step == 0:
                self.heat_print(step)

            for i in range(len(t1)):
                self.nodes[i].temperature = t1[i]

            print('Time {}'.format(step + dt))
            print('Maximum temperature: {}'.format(max(t1)))
            print('Minimum temperature: {}'.format(min(t1)))
            max_temp_list.append(max(t1))
            min_temp_list.append(min(t1))
            self.heat_print(step)
            print()

    def heat_print(self, step):
        nodes_heatmap = zeros([self.nH, self.nL])
        index = 0
        for i in range(self.nH):
            for j in range(self.nL):
                nodes_heatmap.itemset((i, j), self.nodes[index].temperature)
                index += 1
        plt.figure(figsize=(self.nH + 4, self.nL + 4))
        sb.heatmap(nodes_heatmap, 0, self.ambient_temperature - 300, annot=True, square=False, fmt="f", cmap="coolwarm",
                   yticklabels='',
                   xticklabels='')
        plt.savefig('chart/time{}.png'.format(step + 50))
