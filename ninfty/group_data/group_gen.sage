G = PermutationGroup(gap_group = gap.SmallGroup(12,3).AsPermGroup())

# Set the file name
group_name = "A4"

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
    
for i, H in enumerate(sub):
    print("//Conjugation: Building subgroup ", i, " out of", len(sub))
    for g in G:
        Hg = gap_dict_subgroups[H]^gap_dict_elements[g]
        for K in sub:
            if gap_dict_subgroups[K] == Hg:
                sageHg = K
                break
        conj_data[(H,g)] = K

conj_complex = {}

# Compute the conjugates
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
conjugate_string += "\n};\n \n"

dictionary_string = "std::vector<std::string> subgroup_dictionary{\n"
for i in range(0,len(sub)):
    
    dictionary_string += "\"" + str(sub[i].structure_description()) + "\",\n"
dictionary_string = dictionary_string[:-2]
dictionary_string += "\n};\n\n"

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

# Want to find the transitive closures for i=(A,B) and j=(C,D) the entry [i][j] gives the k which is the transitive closure, or min(i,j) else.
composite_string = "std::vector<std::vector<unsigned>> transitive_closure{\n"
for num, i in enumerate(Q):
    print("//Composite: Checking edge ", num, "out of", len(Q)-1)
    composite_string += "{"
    for j in Q:
        ind = min(Q.index(i),Q.index(j))
        if i[1] == j[0] and i != j:
            test_el = [i[0],j[1]]
            if test_el in Q:
                ind = Q.index(test_el)
        composite_string += str(ind) + ","
    composite_string = composite_string[:-1]
    composite_string += "},\n"
composite_string = composite_string[:-2]
composite_string += "\n};"

with open(group_name + ".h" , 'w') as f:
    f.write("#include <vector>\n")
    f.write(dictionary_string)
    f.write("\n")
    f.write(lattice_string)
    f.write("\n")
    f.write(conjugate_string)
    f.write("\n")
    f.write(intersection_string)
    f.write("\n")
    f.write(composite_string)
