G0 = gap.SmallGroup(60,5)
group_name = "A5"


G = PermutationGroup(gap_group = gap.Image(gap.IsomorphismPermGroup(G0)))
sub = G.subgroups()
f = lambda h,k: h.is_subgroup(k)
P = Poset((sub,f))
L = LatticePoset(P)
R = L.relations()

Q = []
for edge in R:
    if edge[0] != edge[1]:
        Q.append(edge)

print("//Starting Computation. There are ", len(Q), " edges")

conjugate_string = "std::vector<std::vector<unsigned>> conjugates{\n"

conj_data = {}
gap_dict_subgroups = {}
gap_dict_elements = {}

for H in sub:
    gap_dict_subgroups[H] = gap(H)
    
for g in G:
    gap_dict_elements[g] = gap(g)

if not gap.IsAbelian(gap(G)):
    for i, H in enumerate(sub):
        print("//Conjugation: Building subgroup ", i, " out of", len(sub))
        for g in G:
            Hg = gap_dict_subgroups[H]^gap_dict_elements[g]
            for K in sub:
                if gap_dict_subgroups[K] == Hg:
                    sageHg = K
                    break
            conj_data[(H,g)] = K
else:
    for i, H in enumerate(sub):
        for g in G:
            conj_data[(H,g)] = H
    
conj_complex = {}

# Compute the conjugates

if not gap.IsAbelian(gap(G)):
    for i, edge in enumerate(Q):
        print("//Conjugation: Checking edge ", i, "out of", len(Q))
        conjugate_string += "{"
        if i not in conj_complex:
            conj_store = [Q.index(edge)]
            for g in G:
                sageHg = conj_data[(edge[0],g)]
                sageKg = conj_data[(edge[1],g)]
                HgKgedge = [sageHg, sageKg]
                if edge != HgKgedge:
                    conj_store.append(Q.index(HgKgedge))
            conj_store = list(dict.fromkeys(conj_store))
            conj_complex[i] = conj_store
            for el in conj_complex[i]:
                if el not in conj_complex:
                    conj_complex[el] = conj_complex[i].copy()
        for p in conj_complex[i]:
            conjugate_string += str(p) + ","
        if conjugate_string[-1] != '{':
            conjugate_string = conjugate_string[:-1]
        conjugate_string += "},\n"
    conjugate_string = conjugate_string[:-2]
    conjugate_string += "};\n \n"

else:
    for i, edge in enumerate(Q):
        conjugate_string += "{" + str(i) + "},\n"
    conjugate_string = conjugate_string[:-2]
    conjugate_string += "};\n \n"

dictionary_string = "std::vector<std::string> subgroup_dictionary{\n"
for i in range(0,len(sub)):
    
    dictionary_string += "\"" + str(sub[i].structure_description()) + "\",\n"
dictionary_string = dictionary_string[:-2]
dictionary_string += "\n};\n"

lattice_string = "std::vector<std::pair<unsigned, unsigned>> lattice{\n"
for el in Q:
    lattice_string += "{" + str(sub.index(el[0])) + "," + str(sub.index(el[1])) + "},\n"
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

with open(group_name + ".h" , 'w') as f:
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
