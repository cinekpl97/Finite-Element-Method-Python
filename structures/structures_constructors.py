from shapes.matrix import matrix_h_local, matrix_c, vector_p
from structures.element import Element
from structures.mesh import Mesh
from structures.point import Point
from structures.node import Node
from shapes import matrix
from settings import *


class Structures:
    def __init__(self):
        settings = get_settings()
        self.mesh = Mesh()
        self.H = settings["H"]
        self.L = settings["L"]
        self.nH = settings["nH"]
        self.nL = settings["nL"]
        self.k = settings["k"]
        self.ro = settings["ro"]
        self.number_of_elements = (self.nH - 1) * (self.nL - 1)
        self.base_temperature = settings["base_temperature"]

        # creating nodes
        dH = self.H / (self.nH - 1)
        dL = self.L / (self.nL - 1)
        print(dH)
        print(dL)
        id_temp = 1
        x = 0
        y = 0
        for i in range(self.nL):
            y = 0
            for j in range(self.nH):
                temp_node = Node(x, y, id_temp, self.base_temperature)
                self.mesh.nodes.append(temp_node)
                y += dH
                id_temp += 1
            x += dL

        # creating elements
        id_element = 1
        id1 = -1
        for i in range(self.nL - 1):
            id1 = id1 + 1
            id2 = id1 + self.nH
            id3 = id2 + 1
            id4 = id1 + 1

            for j in range(self.nH - 1):
                element = Element(id_element,
                                  self.mesh.nodes[id1],
                                  self.mesh.nodes[id2],
                                  self.mesh.nodes[id3],
                                  self.mesh.nodes[id4],
                                  self.k)
                self.mesh.elements.append(element)

                if self.mesh.nodes[id1].y == 0 and self.mesh.nodes[id2].y == 0:
                    element.heated_areas[0] = 1
                if self.mesh.nodes[id2].x == self.L and self.mesh.nodes[id3].x == self.L:
                    element.heated_areas[1] = 1
                if self.mesh.nodes[id3].y == self.H and self.mesh.nodes[id4].y == self.H:
                    element.heated_areas[2] = 1
                if self.mesh.nodes[id4].x == 0 and self.mesh.nodes[id1].x == 0:
                    element.heated_areas[3] = 1

                id_element += 1
                id1 += 1
                id2 += 1
                id3 += 1
                id4 += 1

    def mesh_create(self):

        for single_element in self.mesh.elements:
            single_element.element_matrix_h = matrix_h_local(single_element)
            single_element.element_matrix_c = matrix_c(single_element)
            single_element.element_vector_p = vector_p(single_element)
        self.mesh.aggregation()
