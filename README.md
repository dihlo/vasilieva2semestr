### Необходимо добавить свои файлы в

mesh = Mesh("mesh.xml")

subdomains = MeshFunction("size_t", mesh, "mesh_physical_region.xml")

boundaries = MeshFunction("size_t", mesh, "mesh_facet_region.xml")

для удобства можете называть также