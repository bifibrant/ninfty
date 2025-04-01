L0 = Poset({0:[1,2,3,4,],1:[4,],2:[3,4,],3:[4,],4:[],})
lattice_name = "Pentagon"


L = L0.relabel()
R = L.relations()

N = L.cardinality()

Q = []
for edge in R:
    if edge[0] != edge[1]:
        Q.append(edge)

print("//Starting Computation. There are ", len(Q), " edges")

conjugate_string = "std::vector<std::vector<unsigned>> conjugates{\n"

for i in range(0,len(Q)):
    conjugate_string += "{" + str(i) + "},\n"
conjugate_string = conjugate_string[:-2]
conjugate_string += "};\n \n"

dictionary_string = "std::vector<std::string> subgroup_dictionary{\n"
for i in range(0,N):
    dictionary_string += "\"" + str(i) + "\",\n"
dictionary_string = dictionary_string[:-2]
dictionary_string += "\n};\n"

lattice_string = "std::vector<std::pair<unsigned, unsigned>> lattice{\n"
for el in Q:
    lattice_string += "{" + str(el[0]) + "," + str(el[1]) + "},\n"
#        print(sub[el[0]])
lattice_string = lattice_string[:-2]
lattice_string += "\n}; \n \n"

# This takes in i = (K,H) from the lattice and spits out the unique (L \cap K, L)) = j where L<=H
intersection_string = "std::vector<std::vector<unsigned>> intersections{\n"
for i in Q:
    print("//Intersection: Checking edge ", Q.index(i), "out of", len(Q)-1)
    temp_store = []
    for l in L.list():
        if L.is_less_than(l,i[1]):
            p1 = L.meet(l, i[0])
            p2 = l
            if p1 != p2:
                temp_store.append([p1,p2])
    intersection_string += "{"

    for el in temp_store:
        intersection_string += str(Q.index(el)) + ","
    if len(temp_store) > 0:
        intersection_string = intersection_string[:-1]
    intersection_string += "},\n"
intersection_string = intersection_string[:-2]
intersection_string += "\n}; \n \n"

# This takes in i = (K,H) from the lattice and spits out the unique (L \cap K, L)) = j where L<=H in the dual lattice
# This is essential for the computation of saturated transfer systems (in the non Dedekind case) (Need to add a dedekind check)
cointersection_string = "std::vector<std::vector<unsigned>> cointersections{\n"
for i in Q:
    print("//CoIntersection: Checking edge ", Q.index(i), "out of", len(Q)-1)
    temp_store = []
    for l in L.list():
        if L.is_greater_than(l,i[0]):
            p1 = L.join(l, i[1])
            p2 = l
            if p1 != p2:
                temp_store.append([p2,p1])
    cointersection_string += "{"

    for el in temp_store:
        cointersection_string += str(Q.index(el)) + ","
    if len(temp_store) > 0:
        cointersection_string = cointersection_string[:-1]
    cointersection_string += "},\n"
cointersection_string = cointersection_string[:-2]
cointersection_string += "\n}; \n \n"

with open(lattice_name + ".h" , 'w') as f:
    f.write("#include <vector>\n\n")
    f.write(dictionary_string)
    f.write("\n")
    f.write(lattice_string)
    f.write("\n")
    f.write(conjugate_string)
    f.write("\n")
    f.write(intersection_string)
    f.write("\n")
    f.write(cointersection_string)
    f.write("\n")
    f.write("std::vector<std::string> pretty_subgroup_dictionary{};\nstd::vector<std::string> vertex_layout{};\nstd::vector<std::vector<std::string>> edge_options{};")
