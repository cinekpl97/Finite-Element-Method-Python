from structures.structures_constructors import Structures


def main():
    structures = Structures()
    structures.mesh_create()
    structures.mesh.heat_equation_solving()


if __name__ == '__main__':
    main()
