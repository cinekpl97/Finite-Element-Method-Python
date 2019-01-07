from structures.point import Point
from settings import *
from shapes import shapes_functions_and_deratives
from math import *
from numpy import *

settings = get_settings()
k = settings["k"]
c = settings["c"]
ro = settings["ro"]
alfa = settings["alfa"]
ambient_temperature = settings["ambient_temperature"]


# | dN  |   |   dx   dy    | |dN|
# |dksi |   |  dksi  deta  | |dx|
# |     |   |              | |  |
# |     |=  |   dx   dy    | |dN|
# |dN   |   |  deta  dksi  | |dy|
# |deta |   |______________| |  |
#                   |
#                   v
def jakobian_matrix(element, integration_point):
    dx_dksi, dx_deta, dy_dksi, dy_deta = 0, 0, 0, 0
    dn_ksi = shapes_functions_and_deratives.dn_dksi(integration_point.eta)
    dn_eta = shapes_functions_and_deratives.dn_deta(integration_point.ksi)

    for i in range(4):
        dx_dksi += dn_ksi[i] * element.nodes[i].x
        dx_deta += dn_eta[i] * element.nodes[i].x
        dy_dksi += dn_ksi[i] * element.nodes[i].y
        dy_deta += dn_eta[i] * element.nodes[i].y
    return matrix([[dx_dksi, dx_deta], [dy_dksi, dy_deta]])


# -----------------------
# |          |           |
# |  p4      |      p3   |
# |          |           |
# |----------------------|
# |          |           |
# |  p1      |      p2   |
# |----------|-----------
#            |
#            v
def create_integration_point_set():
    const_value = 1 / sqrt(3)
    integration_point_set = [Point(-const_value, -const_value), Point(const_value, -const_value),
                             Point(const_value, const_value),
                             Point(-const_value, const_value)]
    return integration_point_set


def matrix_h(element):
    integration_point_set = create_integration_point_set()
    dn_dx_dn_dy_and_transponated_multiply_k = []
    for i in range(4):
        integration_point = integration_point_set[i]
        dn_ksi = shapes_functions_and_deratives.dn_dksi(integration_point.eta)
        dn_eta = shapes_functions_and_deratives.dn_deta(integration_point.ksi)

        jakobian = jakobian_matrix(element, integration_point_set[i])

        det_jakobian = linalg.det(jakobian)
        inv_jakobian = jakobian.I

        dn_dx = zeros([4])
        dn_dy = zeros([4])
        for j in range(4):
            dn_dx.itemset(j, inv_jakobian.item(0) * dn_ksi[j] + inv_jakobian.item(1) * dn_eta[j])
            dn_dy.itemset(j, inv_jakobian.item(2) * dn_ksi[j] + inv_jakobian.item(3) * dn_eta[j])

        dn_dx_t = dn_dx.reshape(4, 1)
        dn_dy_t = dn_dy.reshape(4, 1)

        dn_dx_and_dn_dx_t_multiply = dn_dx * dn_dx_t * det_jakobian
        dn_dy_and_dn_dy_t_multiply = dn_dy * dn_dy_t * det_jakobian
        dn_dx_dn_dy_and_transponated_multiply_k.append((dn_dx_and_dn_dx_t_multiply + dn_dy_and_dn_dy_t_multiply) * k)

    matrix_h_sum = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            temp = 0
            for single_integration_point in range(4):
                temp += dn_dx_dn_dy_and_transponated_multiply_k[single_integration_point].item((i, j))
                matrix_h_sum.itemset((i, j), temp)
    print(matrix_h_sum)
    return matrix_h_sum


def matrix_c(element):
    matrix_n_matrix_n_t_multiply = []
    integration_point_set = create_integration_point_set()
    for i in range(4):
        integration_point = integration_point_set[i]

        n = shapes_functions_and_deratives.shape_functions_create(integration_point.ksi, integration_point.eta)

        jakobian = jakobian_matrix(element, integration_point_set[i])

        det_jakobian = linalg.det(jakobian)
        inv_jakobian = jakobian.I
        matrix_n = zeros(4)
        for j in range(4):
            matrix_n.itemset(j, n[j])
        matrix_n_t = matrix_n.reshape(4, 1)
        matrix_n_multiplied = matrix_n * matrix_n_t * ro * c * det_jakobian
        matrix_n_matrix_n_t_multiply.append(matrix_n_multiplied)
    matrix_c_sum = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            temp = 0
            for single_integration_point in range(4):
                temp += matrix_n_matrix_n_t_multiply[single_integration_point].item((i, j))
                matrix_c_sum.itemset((i, j), temp)
    print(matrix_c_sum)
    return matrix_c_sum


def create_integration_point_bc_set():
    integration_point_set = [Point(-1 / sqrt(3), -1), Point(1 / sqrt(3), -1), Point(1, -1 / sqrt(3)),
                             Point(1, 1 / sqrt(3)), Point(1 / sqrt(3), 1), Point(-1 / sqrt(3), 1),
                             Point(-1, 1 / sqrt(3)), Point(-1, -1 / sqrt(3))]
    return integration_point_set


# In element we've got area 1, 2, 3, 4 like this
#          3
# 4     element     2
#          1
#
# Every area has two integration points laying on it
# depending on from where the heat comes, we count integration points for this area/surface

def matrix_h_bc_2d(element):
    integration_point_set = create_integration_point_bc_set()

    det_j_area = [(element.nodes[1].x - element.nodes[0].x) * 0.5, (element.nodes[2].y - element.nodes[1].y) * 0.5,
                  (element.nodes[2].x - element.nodes[3].x) * 0.5, (element.nodes[3].y - element.nodes[0].y) * 0.5]

    areas = [{"area_number": 1, "det_j": det_j_area[0],
              "areas_integration_points": [integration_point_set[0], integration_point_set[1]],
              "heat": element.heated_areas[0]},

             {"area_number": 2, "det_j": det_j_area[1],
              "areas_integration_points": [integration_point_set[2], integration_point_set[3]],
              "heat": element.heated_areas[1]},

             {"area_number": 3, "det_j": det_j_area[2],
              "areas_integration_points": [integration_point_set[4], integration_point_set[5]],
              "heat": element.heated_areas[2]},

             {"area_number": 4, "det_j": det_j_area[3],
              "areas_integration_points": [integration_point_set[6], integration_point_set[7]],
              "heat": element.heated_areas[3]}]

    matrix_multiply = []
    for single_area in areas:
        matrix_n_sum = zeros([4, 4])
        for integration_point in single_area["areas_integration_points"]:
            vector_n = zeros(4)
            for i in range(4):
                vector_n.itemset(i, shapes_functions_and_deratives.shape_functions_create(
                    integration_point.ksi, integration_point.eta)[i])
            vector_n_t = vector_n.reshape(4, 1)
            matrix_n_sum += vector_n * vector_n_t * alfa

        matrix_n_sum *= single_area["det_j"] * single_area["heat"]
        matrix_multiply.append(matrix_n_sum)
    matrix_h_bc = zeros([4, 4])
    for matrix_index in range(len(matrix_multiply)):
        matrix_h_bc += matrix_multiply[matrix_index]
    return matrix_h_bc


def matrix_h_local(element):
    return matrix_h(element) + matrix_h_bc_2d(element)


def vector_p(element):
    integration_point_set = create_integration_point_bc_set()

    det_j_area = [(element.nodes[1].x - element.nodes[0].x) * 0.5, (element.nodes[2].y - element.nodes[1].y) * 0.5,
                  (element.nodes[2].x - element.nodes[3].x) * 0.5, (element.nodes[3].y - element.nodes[0].y) * 0.5]

    areas = [{"area_number": 1, "det_j": det_j_area[0],
              "areas_integration_points": [integration_point_set[0], integration_point_set[1]],
              "heat": element.heated_areas[0]},

             {"area_number": 2, "det_j": det_j_area[1],
              "areas_integration_points": [integration_point_set[2], integration_point_set[3]],
              "heat": element.heated_areas[1]},

             {"area_number": 3, "det_j": det_j_area[2],
              "areas_integration_points": [integration_point_set[4], integration_point_set[5]],
              "heat": element.heated_areas[2]},

             {"area_number": 4, "det_j": det_j_area[3],
              "areas_integration_points": [integration_point_set[6], integration_point_set[7]],
              "heat": element.heated_areas[3]}]

    vector_p_multiply = []
    for single_area in areas:
        vector_p_sum = zeros(4)
        for integration_point in single_area["areas_integration_points"]:
            vector_p_ = zeros(4)
            for i in range(4):
                vector_p_.itemset(i, shapes_functions_and_deratives.shape_functions_create(
                    integration_point.ksi, integration_point.eta)[i])
            vector_p_sum += vector_p_

        vector_p_sum *= single_area["det_j"] * single_area["heat"] * (-alfa) * ambient_temperature
        vector_p_multiply.append(vector_p_sum)
    vector_p_final = zeros(4)
    for i in range(len(vector_p_multiply)):
        vector_p_final += vector_p_multiply[i]
    print(vector_p_final)
    return vector_p_final
