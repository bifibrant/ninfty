//
//  ninfty.h
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#ifndef ninfty_h
#define ninfty_h

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <chrono>
#include <unordered_set>
#include <future>
#include <random>
#include <regex>
#include <iostream>

enum GenerationType{ ALL, SATURATED, COSATURATED, UNDERLYING, CONJUGACY };

// A function which determines if one std::vector<T> is contained in another
template <typename T>
bool isSubsetOrEqual(std::vector<T> const& a, std::vector<T> const& b) {
    for(auto const& av:a){
        if(std::find(b.begin(),b.end(),av)==b.end())
            return false;
    }
    return true;
}

// Hash function for vectors of unsigned int.
// Taken from https://stackoverflow.com/a/12996028
class unsigned_vector_hasher {
public:
    std::size_t operator()(std::vector<unsigned> const& vec) const {
        std::size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Variable which stores the conjugacy classes of subgroups
std::vector<std::vector<unsigned>> CONJUGACY_CLASSES;

// Variable which stores the array of meets
std::vector<std::vector<unsigned>> meetArray;

// Variables which stores all requested transfer systems
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> RESULT;
std::vector<std::vector<unsigned>> ALL_STORE;
std::vector<std::vector<unsigned>> SATURATED_STORE;
std::vector<std::vector<unsigned>> OPPOSITE_SATURATED_STORE;
std::vector<std::vector<unsigned>> COSATURATED_STORE;
std::vector<std::vector<unsigned>> UNDERLYING_STORE;
std::vector<std::vector<unsigned>> CONJUGACY_STORE;

// Varaible which stores the maximally generated transfer systems
std::vector<std::vector<unsigned>> MAXIMALLY_GENERATED;
std::vector<std::vector<unsigned>> FLAT_STORE;
std::vector<std::vector<unsigned>> LSP_STORE;

// Varible which stores the inclusion lattice for transfer systems
std::vector<std::pair<unsigned,unsigned>> TRANSFER_LATTICE;

// Variables which store various pairwise checks
// Stores indices of TRANSFER_LATTICE
std::vector<unsigned> COMPATIBLE_PAIRS;
std::vector<unsigned> CCLOSED_PAIRS;
std::vector<unsigned> QUILLEN_PAIRS;

// Variable which stores the generation statistics for transfer systems
std::vector<unsigned> GENERATION_STATISTICS;

// Variable which stores the complexity of the requested computation
unsigned ALL_COMPLEXITY = 0;
unsigned SATURATED_COMPLEXITY = 0;
unsigned COSATURATED_COMPLEXITY = 0;
unsigned UNDERLYING_COMPLEXITY = 0;
unsigned CONJUGACY_COMPLEXITY = 0;

// Setting up global variables and constants
unsigned NUM_THREADS = unsigned(lattice.size());
std::vector<std::vector<unsigned>> NEW_EDGES;
std::vector<std::shared_future<std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher>>> THREAD_STORE;
std::vector<unsigned> GOOD_EDGES;
std::vector<std::vector<unsigned>> transitive_closure;
std::vector<unsigned> edgesFromE;
std::vector<unsigned> edgesToG;

// A function to check if the given group is Dedekind
bool isDedekind(){
    for(unsigned i = 0; i < conjugates.size(); ++i){
        if(conjugates[i].size() > 1){
            return false;
        }
    }
    return true;
}

// An algorithm which finds (and updates) the edges which are used for cosaturation generation
void findCosaturationEdges(){
    edgesToG.clear();
    for(unsigned i=0; i<lattice.size(); ++i){
        if(lattice[i].second == subgroup_dictionary.size()-1){
            edgesToG.push_back(i);
        }
    }
}

// An algorithm which finds (and updates) the edges which are used for saturation generation
void findSaturationEdges(){
    edgesFromE.clear();
    for(unsigned i=0; i<lattice.size(); ++i){
        if(lattice[i].first == 0){
            edgesFromE.push_back(i);
        }
    }
}

// An algorithm which finds (and updates) the transitive closure of the lattice elements
void findTransitiveClosure(){
    for(unsigned i=0; i<lattice.size(); ++i){
        std::vector<unsigned> i_row;
        for(unsigned j=0; j<lattice.size(); ++j){
            unsigned pos = std::min(i,j);
            if(lattice[i].second == lattice[j].first & i != j){
                std::pair<unsigned,unsigned> test_el{lattice[i].first, lattice[j].second};
                auto it = std::find(lattice.begin(), lattice.end(), test_el);
                if(it != lattice.end()){
                    pos = unsigned(std::distance(lattice.begin(), it));
                }
            }
            i_row.push_back(pos);
        }
        transitive_closure.push_back(i_row);
    }
}

// An algorithm to comptue the meet of two subgroups A and B
unsigned computeMeet(const unsigned& A, const unsigned& B){
    if(A == B){
        return A;
    }
    std::pair<unsigned,unsigned> test_element{A,B};
    if(std::find(lattice.begin(),lattice.end(),test_element) != lattice.end()){
        return A;
    }
    test_element = {B,A};
    if(std::find(lattice.begin(),lattice.end(),test_element) != lattice.end()){
        return B;
    }
    
    std::vector<unsigned> common_lower_elements;
    
    for(unsigned i=0; i<subgroup_dictionary.size(); ++i){
        std::pair<unsigned,unsigned> temp1{i,A};
        std::pair<unsigned,unsigned> temp2{i,B};
        if(std::find(lattice.begin(),lattice.end(),temp1) != lattice.end() & std::find(lattice.begin(),lattice.end(),temp2) != lattice.end()){
            common_lower_elements.push_back(i);
        }
    }
    
    //The unique maximal element will be the one that doesn't appear as the target of any of the edges
    for(unsigned i=0; i<common_lower_elements.size(); ++i){
        bool is_maximal = true;
        for(unsigned j=0; j<common_lower_elements.size(); ++j){
            std::pair<unsigned,unsigned> temp3{common_lower_elements[i], common_lower_elements[j]};
            if(std::find(lattice.begin(),lattice.end(),temp3) != lattice.end()){
                is_maximal = false;
                break;
            }
        }
        if(is_maximal){
            return common_lower_elements[i];
        }
    }
    return 0;
}

// A function which constructs all meets
void computeMeetArray(){
    meetArray.clear();
    for(unsigned i = 0; i < subgroup_dictionary.size(); ++i){
        std::vector<unsigned> temp;
        for(unsigned j = 0; j < subgroup_dictionary.size(); ++j){
            temp.push_back(computeMeet(i, j));
        }
        meetArray.push_back(temp);
    }
}

// An implementation of Rubin's algorithm to find the closure of a collection of norm maps
// Adapted from https://arxiv.org/pdf/1903.08723 Construction B.1
// The enum is used to see what type of generation is used
std::vector<unsigned> transferClosure(const std::vector<unsigned>& r0, const GenerationType& gen_type = ALL){
    std::set<unsigned> r1_set, r2_set, r4_set;
    
    // We check to see if we have computed the transitive closure or not
    if(transitive_closure.size() == 0){
        findTransitiveClosure();
    }
    
    // Close under conjugation
    for(auto it = r0.begin(); it != r0.end(); ++it){
        if(gen_type != UNDERLYING & gen_type != CONJUGACY){
            r1_set.insert(conjugates[(*it)].begin(), conjugates[(*it)].end());
        }
        else{
            r1_set.insert(*it);
        }
    }
    
    
    // Close under intersection
    r2_set = r1_set;
    if(gen_type == SATURATED){
        for(auto it = r1_set.begin(); it != r1_set.end(); ++it){
            r2_set.insert(cointersections[(*it)].begin(), cointersections[(*it)].end());
        }
    }
    else{
        for(auto it = r1_set.begin(); it != r1_set.end(); ++it){
            r2_set.insert(intersections[(*it)].begin(), intersections[(*it)].end());
        }
    }
    
    //Close under composition
    unsigned difference = 1;
    unsigned old_size = 0;
    
    while(difference != 0){
        r4_set.clear();
        old_size = unsigned(r2_set.size());
        
        for(auto it1 = r2_set.begin(); it1 != r2_set.end(); ++it1){
            for(auto it2 = it1; it2 != r2_set.end(); ++it2){
                r4_set.insert(transitive_closure[(*it1)][(*it2)]);
            }
        }
        r2_set.insert(r4_set.begin(), r4_set.end());
        
        difference = unsigned(r2_set.size()) - old_size;
    }
    
    std::vector<unsigned> result;
    result.assign(r2_set.begin(), r2_set.end());
    return result;
}

//The thread function which adds new norms to existing transfer systems
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> threadProcess(const unsigned& start_index, const unsigned& end_index, const GenerationType& gen_type = ALL){
    
    std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> temp_new_edges;
    
    for(auto it = NEW_EDGES.begin(); it != NEW_EDGES.end(); ++it){
        for(unsigned i=start_index; i<end_index; ++i){
            auto test = (*it);
            
            if(std::find(test.begin(), test.end(), GOOD_EDGES[i]) == test.end()){
                test.push_back(GOOD_EDGES[i]);
                
                auto closed = transferClosure(test, gen_type);
                
                if(!RESULT.contains(closed)){
                    temp_new_edges.insert(closed);
                }
            }
        }
    }
    return temp_new_edges;
}

// A function which finds all transfer systems for the given group as defined in https://arxiv.org/abs/1905.03797
// Algorithm adapted from A B B S W Z-C
// The verbose toggle can be used to supress generation statistics
void transferFind(const bool verbose = true, const GenerationType& gen_type = ALL){
    
    findSaturationEdges();
    findCosaturationEdges();
    computeMeetArray();
    
    unsigned long old_size = 0;
    unsigned long diff = 1;
    unsigned gen_step = 1;
    
    // Initalise the algorithm with the empty transfer system
    std::vector<unsigned> empty_transfer{};
    NEW_EDGES.push_back(empty_transfer);
    RESULT.insert(empty_transfer);
    GOOD_EDGES.clear();
    GENERATION_STATISTICS.clear();
    
    // Manually do the first iteration to find the basis edges, we call these edges "good"
    std::vector<std::vector<unsigned>> temp_new_edges;
    
    auto it = NEW_EDGES.begin();
    for(unsigned i=0; i<lattice.size(); ++i){
        if(gen_type == ALL || gen_type == UNDERLYING || gen_type == CONJUGACY){
            std::vector<unsigned> test{i};
            test.insert(std::end(test), std::begin((*it)), std::end((*it)));
            auto closed = transferClosure(test, gen_type);
            auto storage_size_old = RESULT.size();
            RESULT.insert(closed);
            auto storage_size_new = RESULT.size();
            
            if(storage_size_new > storage_size_old){
                temp_new_edges.push_back(closed);
                GOOD_EDGES.push_back(i);
            }
        }
        else if((gen_type == SATURATED) & std::find(edgesFromE.begin(), edgesFromE.end(), i) != edgesFromE.end()){
            std::vector<unsigned> test{i};
            test.insert(std::end(test), std::begin((*it)), std::end((*it)));
            auto closed = transferClosure(test, gen_type);
            auto storage_size_old = RESULT.size();
            RESULT.insert(closed);
            auto storage_size_new = RESULT.size();
            
            if(storage_size_new > storage_size_old){
                temp_new_edges.push_back(closed);
                GOOD_EDGES.push_back(i);
            }
        }
        else if(((gen_type == COSATURATED) & std::find(edgesToG.begin(), edgesToG.end(), i) != edgesToG.end())){
            std::vector<unsigned> test{i};
            test.insert(std::end(test), std::begin((*it)), std::end((*it)));
            auto closed = transferClosure(test, gen_type);
            auto storage_size_old = RESULT.size();
            RESULT.insert(closed);
            auto storage_size_new = RESULT.size();
            
            if(storage_size_new > storage_size_old){
                temp_new_edges.push_back(closed);
                GOOD_EDGES.push_back(i);
            }
        }
    }
    
    NEW_EDGES = temp_new_edges;
    
    //We intialise a thread for each good edge
    NUM_THREADS = unsigned(GOOD_EDGES.size());
    unsigned step = ceil(double(GOOD_EDGES.size()) / double(NUM_THREADS));
    
    if(verbose){
        std::cout << "Algorithm Step 0: 1" << std::endl;
    }
    
    GENERATION_STATISTICS.push_back(unsigned(GOOD_EDGES.size()));
    
    //The recursive algorithm for finding all transfer systems as closed sets
    while(diff != 0){
        
        if(verbose){
            std::cout << "Algorithm Step " << gen_step << ": " << RESULT.size() << std::endl;
        }
        
        THREAD_STORE.clear();
        for(unsigned i=0; i<NUM_THREADS; ++i){
            std::shared_future<std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher>> new_thread = std::async(std::launch::async, threadProcess, i*step, std::min((i+1)*step, unsigned(lattice.size())), gen_type);
            THREAD_STORE.push_back(new_thread);
        }
        
        std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> new_transfer_systems;
        
        for(unsigned i=0; i<NUM_THREADS; ++i){
            new_transfer_systems.insert(THREAD_STORE[i].get().begin(), THREAD_STORE[i].get().end());
        }
        
        NEW_EDGES.clear();
        std::copy(new_transfer_systems.begin(), new_transfer_systems.end(), std::back_inserter(NEW_EDGES));
        old_size = RESULT.size();
        
        for(auto it = NEW_EDGES.begin(); it != NEW_EDGES.end(); ++it){
            RESULT.insert(*it);
        }
        
        diff = RESULT.size() - old_size;
        GENERATION_STATISTICS.push_back(unsigned(diff));
        gen_step++;
        
        // This will store the maximally generated transfer systems
        if(NEW_EDGES.size() < 100 & NEW_EDGES.size() != 0){
            MAXIMALLY_GENERATED = NEW_EDGES;
        }
        
    }
    // Store the generation step as the complexity, and sort the arrays for ease of browsing later
    if(gen_type == ALL){
        ALL_COMPLEXITY = gen_step-1;
        ALL_STORE.clear();
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(ALL_STORE));
        std::sort(ALL_STORE.begin(), ALL_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == SATURATED){
        SATURATED_COMPLEXITY = gen_step-1;
        OPPOSITE_SATURATED_STORE.clear();
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(OPPOSITE_SATURATED_STORE));
        std::sort(OPPOSITE_SATURATED_STORE.begin(), OPPOSITE_SATURATED_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == COSATURATED){
        COSATURATED_COMPLEXITY = gen_step-1;
        COSATURATED_STORE.clear();
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(COSATURATED_STORE));
        std::sort(COSATURATED_STORE.begin(), COSATURATED_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == UNDERLYING){
        UNDERLYING_COMPLEXITY = gen_step-1;
        UNDERLYING_STORE.clear();
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(UNDERLYING_STORE));
        std::sort(UNDERLYING_STORE.begin(), UNDERLYING_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == CONJUGACY){
        CONJUGACY_COMPLEXITY = gen_step-1;
        CONJUGACY_STORE.clear();
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(CONJUGACY_STORE));
        std::sort(CONJUGACY_STORE.begin(), CONJUGACY_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    RESULT.clear();
    NEW_EDGES.clear();
}

std::vector<std::vector<unsigned>> allTransfers(){
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    return ALL_STORE;
}

// Need to check for the two-out-of-three property for saturation
bool isSaturated(const std::vector<unsigned>& rhs){
    for(unsigned i=0; i<rhs.size(); ++i){
        for(unsigned j=0; j<rhs.size(); ++j){
            if(lattice[rhs[i]].first == lattice[rhs[j]].first and lattice[rhs[i]].second != lattice[rhs[j]].second){
                std::pair<unsigned,unsigned> to_find{lattice[rhs[i]].second, lattice[rhs[j]].second};
                if(std::find(lattice.begin(), lattice.end(), to_find) != lattice.end()){
                    auto it = find(lattice.begin(), lattice.end(), to_find);
                    unsigned index = unsigned(it - lattice.begin());
                    if(std::find(rhs.begin(), rhs.end(), index) == rhs.end()){
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

std::vector<std::vector<unsigned>> saturatedTransfers(){
    SATURATED_STORE.clear();
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    for(unsigned i = 0; i <ALL_STORE.size(); ++i){
        if(isSaturated(ALL_STORE[i])){
            SATURATED_STORE.push_back(ALL_STORE[i]);
        }
    }
    return SATURATED_STORE;
}

std::vector<std::vector<unsigned>> cosaturatedTransfers(){
    if(COSATURATED_STORE.size() == 0){
        transferFind(false, COSATURATED);
    }
    return COSATURATED_STORE;
}

// An implementation of Rubin's algorithm in reverse providing a set of minimal generators for a given transfer system
std::vector<unsigned> findBasis(const std::vector<unsigned>& rhs){
    auto r2 = rhs;
    
    // We check to see if we have computed the transitive closure or not
    if(transitive_closure.size() == 0){
        findTransitiveClosure();
    }
    
    // Remove any edges which can be expressed as a composite
    for(unsigned i = 0; i<rhs.size(); ++i){
        bool is_composite = false;
        for(unsigned k=0; k<rhs.size(); ++k){
            for(unsigned l=0; l<rhs.size(); ++l){
                //check to see if rhs[i] = rhs[k]rhs[l]
                if(rhs[i] == transitive_closure[rhs[k]][rhs[l]] and rhs[i] != rhs[k] and rhs[i] != rhs[l]){
                    is_composite = true;
                    break;
                }
            }
        }
        if(is_composite){
            auto it = std::find(r2.begin(), r2.end(), rhs[i]);
            r2.erase(it);
        }
    }
    auto r3 = r2;
    
    // Remove any edges which are induced via the restriction axiom
    for(unsigned i=0; i<r2.size(); ++i){
        auto ints = intersections[r2[i]];
        for(unsigned k=0; k<ints.size(); ++k){
            auto it = std::find(r3.begin(), r3.end(), ints[k]);
            if(it != r3.end()){r3.erase(it);}
        }
    }
    
    auto r4 = r3;
    
    // Make a choice of conjugacy classes for edges
    for(unsigned i=0; i<r3.size(); ++i){
        if(std::find(r4.begin(), r4.end(), r3[i]) != r4.end()){
            auto conjs = conjugates[r3[i]];
            for(unsigned k=1; k<conjs.size(); ++k){
                if(std::find(r4.begin(), r4.end(), conjs[k]) != r4.end()){
                    auto it = std::find(r4.begin(), r4.end(), conjs[k]);
                    r4.erase(it);
                }
            }
        }
    }
    return r4;
}

// A function which finds the width of a group
// This is the size of a minimal basis for the complete transfer system
unsigned width(){
    std::vector<unsigned> complete;
    for(unsigned i = 0; i < lattice.size(); i++)
        complete.push_back(i);
    return unsigned(findBasis(complete).size());
}

// A function which finds the complexity of a group
// This is the maxima, size of any minimal basis
unsigned complexity(const GenerationType& gen_type = ALL){
    if(gen_type == ALL){
        if(ALL_COMPLEXITY != 0){
            return ALL_COMPLEXITY;
        }
        else{
            transferFind(false, ALL);
            return ALL_COMPLEXITY;
        }
    }
    else if(gen_type == SATURATED){
        if(SATURATED_COMPLEXITY != 0){
            return SATURATED_COMPLEXITY;
        }
        else{
            transferFind(false, SATURATED);
            return SATURATED_COMPLEXITY;
        }
    }
    else if(gen_type == COSATURATED){
        if(COSATURATED_COMPLEXITY != 0){
            return COSATURATED_COMPLEXITY;
        }
        else{
            transferFind(false, COSATURATED);
            return COSATURATED_COMPLEXITY;
        }
    }
    return 0;
}


// A function which checks if a given transfer system is cosaturated
bool isCosaturated(const std::vector<unsigned>& rhs){
    std::vector<unsigned> test_basis;
    for(unsigned i=0; i<rhs.size(); ++i){
        if(lattice[rhs[i]].second == subgroup_dictionary.size()-1){
            test_basis.push_back(rhs[i]);
        }
    }
    return (transferClosure(test_basis, COSATURATED) == rhs);
}

// A thread function that checks inclusions of trasnfer systems
std::vector<std::pair<unsigned,unsigned>> batchTransferLattice(const unsigned& start_index, const unsigned& end_index){
    std::vector<std::pair<unsigned,unsigned>> resultant_lattice;
    
    for(unsigned i=start_index; i<end_index; ++i){
        for(unsigned j=i; j<ALL_STORE.size(); ++j){
            if(isSubsetOrEqual(ALL_STORE[i], ALL_STORE[j])){
                resultant_lattice.push_back({i,j});
            }
        }
    }
    return resultant_lattice;
}

// A function which returns pairs (i,j) such that ALL_TRANSFERS[i] <= ALL_TRANSFERS[j]. Required for model structures and compatible pairs
std::vector<std::pair<unsigned,unsigned>> transferLattice(){
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    
    std::vector<std::shared_future<std::vector<std::pair<unsigned,unsigned>>>> LATTICE_THREAD_STORE;
    unsigned step = ceil(double(ALL_STORE.size()) / double(NUM_THREADS));
    for(unsigned i=0; i<NUM_THREADS; ++i){
        std::shared_future<std::vector<std::pair<unsigned,unsigned>>> new_thread = std::async(std::launch::async, batchTransferLattice, i*step, std::min((i+1)*step, unsigned(ALL_STORE.size())));
        LATTICE_THREAD_STORE.push_back(new_thread);
    }
    std::vector<std::pair<unsigned,unsigned>> result;
    
    for(unsigned i = 0; i < NUM_THREADS; ++i){
        result.insert(result.end(), LATTICE_THREAD_STORE[i].get().begin(),LATTICE_THREAD_STORE[i].get().end());
    }
    TRANSFER_LATTICE = result;
    return result;
}

// A function which determines if a pair of transfer systems is compatible in the sense of Chan
bool isCompatible(const std::vector<unsigned>& transfer_m, const std::vector<unsigned>& transfer_a){
    
    if(transfer_m.size() == 0){return true;}
    if(transfer_a.size() == lattice.size()){return true;}
    
    for(unsigned i=0; i<transfer_m.size(); ++i){
        unsigned B = lattice[transfer_m[i]].first;
        unsigned A = lattice[transfer_m[i]].second;
        // As the subgroups are ordered by size we need only check those subgroups that appear before A
        for(unsigned C=0; C<A; ++C){
            // Construct the potential edge C->A
            std::pair<unsigned,unsigned> temp{C,A};
            if(std::find(lattice.begin(), lattice.end(), temp) != lattice.end()){
                
                // We now want to check if B \cap C is in transfer_a
                std::pair<unsigned, unsigned> test_edge{meetArray[C][B],B};
                unsigned test_edge_index = unsigned(std::find(lattice.begin(), lattice.end(), test_edge) - lattice.begin());
                
                if((std::find(transfer_a.begin(), transfer_a.end(), test_edge_index) != transfer_a.end()) || ((meetArray[B][C] == B))){
                    // Now check that C -> A is in transfer_a
                    std::pair<unsigned, unsigned> test_edge2{C,A};
                    unsigned test_edge_index2 = unsigned(std::find(lattice.begin(), lattice.end(), test_edge2) - lattice.begin());
                    // If we do not find this transfer then we are not compatible
                    if(std::find(transfer_a.begin(), transfer_a.end(), test_edge_index2) == transfer_a.end()){
                        return false;
                    }
                }
                
            }
        }
    }
    return true;
}

// A thread function that checks compatibility
std::vector<unsigned> batchCompatible(const unsigned& start_index, const unsigned& end_index){
    std::vector<unsigned> result;
    
    for(unsigned i=start_index; i<end_index; ++i){
        if(isCompatible(ALL_STORE[TRANSFER_LATTICE[i].first],ALL_STORE[TRANSFER_LATTICE[i].second])){
            result.push_back(i);
        }
    }
    return result;
}

// A function which computes all of the compatible pairs (using parallel)
// The return is the index of the pair in TRANSFER_LATTICE
std::vector<unsigned> compatiblePairs(){
    if(COMPATIBLE_PAIRS.size() != 0){
        return COMPATIBLE_PAIRS;
    }
    
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    std::vector<std::shared_future<std::vector<unsigned>>> COMPATIBLE_THREAD_STORE;
    unsigned step = ceil(double(TRANSFER_LATTICE.size()) / double(NUM_THREADS));
    for(unsigned i=0; i<NUM_THREADS; ++i){
        std::shared_future<std::vector<unsigned>> new_thread = std::async(std::launch::async, batchCompatible, i*step, std::min((i+1)*step, unsigned(TRANSFER_LATTICE.size())));
        COMPATIBLE_THREAD_STORE.push_back(new_thread);
    }
    std::vector<unsigned> result;
    
    for(unsigned i = 0; i < NUM_THREADS; ++i){
        result.insert(result.end(), COMPATIBLE_THREAD_STORE[i].get().begin(),COMPATIBLE_THREAD_STORE[i].get().end());
    }
    
    COMPATIBLE_PAIRS = result;
    return result;
}

// An algorithm which finds the minimal fibrant subgroup of a transfer system in the sense of the minimal normal subgroup which has a transfer to G itself
unsigned minimalFibrantSubgroup(const std::vector<unsigned>& rhs){
    unsigned min_fibrant = unsigned(subgroup_dictionary.size()-1);
    for(unsigned i = 0; i <rhs.size(); ++i){
        // Note that as we are looking for the minimal subgroup with transfer to G, this will coincide with the minimal index in the subgroup_dictionary given that this is ordered by size of subgroup
        if((lattice[rhs[i]].second == subgroup_dictionary.size() - 1) & (lattice[rhs[i]].first < min_fibrant)){
            min_fibrant = lattice[rhs[i]].first;
        }
    }
    return min_fibrant;
}

// An algorithm which determines if a given transfer system is flat or not
bool isFlat(const std::vector<unsigned>& rhs){
    unsigned F = minimalFibrantSubgroup(rhs);
    
    if((F == subgroup_dictionary.size() - 1) & (rhs.size() != 0)){
        return false;
    }
    if(F == 0){
        return true;
    }
    
    // Find all emements below F and store the edges
    std::vector<unsigned> edges_below_F;
    for(unsigned i=0; i < F; ++i){
        std::pair<unsigned,unsigned> test_edge{i,F};
        if(std::find(lattice.begin(), lattice.end(), test_edge) != lattice.end()){
            unsigned test_edge_index = unsigned(std::find(lattice.begin(), lattice.end(), test_edge) - lattice.begin());
            edges_below_F.push_back(test_edge_index);
        }
    }
    
    // By closing up under transferClosure we find the sublattice under F
    auto lattice_below_F = transferClosure(edges_below_F);
    
    //If the intersection of rhs and lattice_below_F is trivial then we are flat
    for(unsigned i=0; i < lattice_below_F.size(); ++i){
        if(std::find(rhs.begin(), rhs.end(), lattice_below_F[i]) != rhs.end()){
            return false;
        }
    }
    return true;
}

// A function which returns the index of the flat transfers in ALL_TRANSFERS
std::vector<unsigned> flatTransfers(){
    FLAT_STORE.clear();
    std::vector<unsigned> result;
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        if(isFlat(ALL_STORE[i])){
            result.push_back(i);
            FLAT_STORE.push_back(ALL_STORE[i]);
        }
        
    }
    
    return result;
}

// A function which computes the left set of a transfer system using the algorithm of FOOQW
std::vector<unsigned> leftSet(const std::vector<unsigned> & T){
    std::set<unsigned> temp;
    temp.insert(T.begin(), T.end());
    //Start by taking the downward extension
    for(unsigned p=0; p<T.size(); ++p){
        for(unsigned i=0; i<lattice.size(); ++i){
            if(lattice[T[p]].second == lattice[i].second){
                //Now check for if lattice[i].first < lattice[T[p]].first, if so we add the index
                std::pair<unsigned, unsigned> to_find{lattice[i].first, lattice[T[p]].first};
                auto it = std::find(lattice.begin(), lattice.end(), to_find);
                if(it != lattice.end()){
                    std::pair<unsigned, unsigned> to_find_2{lattice[i].first, lattice[T[p]].second};
                    auto it2 = std::find(lattice.begin(), lattice.end(), to_find_2);
                    temp.insert(unsigned(it2 - lattice.begin()));
                }
            }
        }
    }
    
    //Now take the complement
    std::vector<unsigned> res;
    for(unsigned i=0; i<lattice.size(); ++i){
        if(!temp.contains(i)){
            res.push_back(i);
        }
    }
    return res;
}

// A function that produces the weak equivalences of an AF AC pair
std::vector<unsigned> weakEquivalences(const std::vector<unsigned>& AF, const std::vector<unsigned>& F){
    std::set<unsigned> store;
    auto AC = leftSet(F);
    store.insert(AC.begin(), AC.end());
    store.insert(AF.begin(), AF.end());
    
    //AF \circ AC
    for(unsigned i=0; i<AC.size(); ++i){
        for(unsigned j=0; j<AF.size(); ++j){
            store.insert(transitive_closure[AC[i]][AF[j]]);
        }
    }
    
    std::vector<unsigned> res;
    res.assign(store.begin(), store.end());
    return res;
}

// A function which checks if an AF AC pair is a CC model structure or model strucutr
// The return is 0 is it is neither, 1 if it is CC closed, and 2 if it is Quillen
unsigned modelCheck(const std::vector<unsigned>& AF, const std::vector<unsigned>& F){
    auto W = weakEquivalences(AF, F);
    
    auto AC = leftSet(F);
    
    // Check it is closed under composition (Using Lemma 3.1 of CClosed)
    std::set<unsigned> store;
    store.insert(AC.begin(), AC.end());
    store.insert(AF.begin(), AF.end());
    
    //AF \circ AC
    for(unsigned i=0; i<AC.size(); ++i){
        for(unsigned j=0; j<AF.size(); ++j){
            store.insert(transitive_closure[AF[j]][AC[i]]);
        }
    }
    std::vector<unsigned> composition_closed_check;
    composition_closed_check.assign(store.begin(), store.end());
    
    if(!isSubsetOrEqual(composition_closed_check, W)){
        return 0;
    }
    
    //We now check for two-out-of-three
    for(unsigned f=0; f<lattice.size(); ++f){
        for(unsigned g=0; g<f; ++g){
            if(lattice[f].first == lattice[g].second){
                unsigned gf = transitive_closure[g][f];
                unsigned count = 0;
                if(std::find(W.begin(), W.end(), f) != W.end()){
                    count++;
                }
                if(std::find(W.begin(), W.end(), g) != W.end()){
                    count++;
                }
                if(std::find(W.begin(), W.end(), gf) != W.end()){
                    count++;
                }
                if(count == 2){
                    return 1;
                }
            }
        }
    }
    return 2;
}

// A thread function that checks model structure properties
// The first vector will contain the indicies in TRANSFER_LATTICE of the CClosed model structures and the second will return the Model structures
std::pair<std::vector<unsigned>,std::vector<unsigned>> batchModel(const unsigned& start_index, const unsigned& end_index){
    std::pair<std::vector<unsigned>,std::vector<unsigned>> result;
    std::vector<unsigned> c_closed;
    std::vector<unsigned> quillen;
    
    for(unsigned i=start_index; i<end_index; ++i){
        unsigned model_check = modelCheck(ALL_STORE[TRANSFER_LATTICE[i].first], ALL_STORE[TRANSFER_LATTICE[i].second]);
        
        if(model_check == 2){
            c_closed.push_back(i);
            quillen.push_back(i);
        }
        else if(model_check == 1){
            c_closed.push_back(i);
        }
    }
    result.first = c_closed;
    result.second = quillen;
    return result;
}

void modelPairs(){
    if(QUILLEN_PAIRS.size() != 0){
        return;
    }
    
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    std::vector<std::shared_future<std::pair<std::vector<unsigned>,std::vector<unsigned>>>> MODEL_THREAD_STORE;
    unsigned step = ceil(double(TRANSFER_LATTICE.size()) / double(NUM_THREADS));
    for(unsigned i=0; i<NUM_THREADS; ++i){
        std::shared_future<std::pair<std::vector<unsigned>,std::vector<unsigned>>> new_thread = std::async(std::launch::async, batchModel, i*step, std::min((i+1)*step, unsigned(TRANSFER_LATTICE.size())));
        MODEL_THREAD_STORE.push_back(new_thread);
    }
    
    
    std::vector<unsigned> result1;
    std::vector<unsigned> result2;
    
    for(unsigned i = 0; i < NUM_THREADS; ++i){
        result1.insert(result1.end(), MODEL_THREAD_STORE[i].get().first.begin(),MODEL_THREAD_STORE[i].get().first.end());
        result2.insert(result2.end(), MODEL_THREAD_STORE[i].get().second.begin(),MODEL_THREAD_STORE[i].get().second.end());
    }
    
    CCLOSED_PAIRS = result1;
    QUILLEN_PAIRS = result2;
}

// A function which computes all of the (unique!) different weak equivalence structures
std::vector<std::vector<unsigned>> weakEquivalenceTypes(){
    std::vector<std::vector<unsigned>> result;
    std::set<std::vector<unsigned>> result_set;
    
    if(QUILLEN_PAIRS.size() == 0){
        modelPairs();
    }
    
    for(unsigned i=0; i<QUILLEN_PAIRS.size(); ++i){
        result_set.insert(weakEquivalences(ALL_STORE[TRANSFER_LATTICE[QUILLEN_PAIRS[i]].first], leftSet(ALL_STORE[TRANSFER_LATTICE[QUILLEN_PAIRS[i]].second])));
    }
    
    result.assign(result_set.begin(), result_set.end());
    return result;
}

// A function which returns the largest cosaturated transfer system inside a given one
std::vector<unsigned> cosaturatedCore(const std::vector<unsigned>& rhs){
    if(isCosaturated(rhs)){
        return rhs;
    }
    
    std::vector<unsigned> cosaturated_edges;
    for(unsigned i = 0; i < rhs.size(); ++i){
        if(lattice[rhs[i]].second == subgroup_dictionary.size() - 1){
            cosaturated_edges.push_back(i);
        }
    }
    return transferClosure(cosaturated_edges);
}

// A function which returns the smallest saturated transfer system containing a given one
std::vector<unsigned> saturatedHull(const std::vector<unsigned>& rhs){
    if(isSaturated(rhs)){
        return rhs;
    }
    
    std::set<unsigned> saturated_edges(rhs.begin(), rhs.end());
    
    // We add in all edges to saturated the result
    for(unsigned i=0; i<rhs.size(); ++i){
        for(unsigned j=0; j<rhs.size(); ++j){
            if(lattice[rhs[i]].first == lattice[rhs[j]].first and lattice[rhs[i]].second != lattice[rhs[j]].second){
                std::pair<unsigned,unsigned> to_find{lattice[rhs[i]].second, lattice[rhs[j]].second};
                if(std::find(lattice.begin(), lattice.end(), to_find) != lattice.end()){
                    auto it = find(lattice.begin(), lattice.end(), to_find);
                    unsigned index = unsigned(it - lattice.begin());
                    saturated_edges.insert(index);
                }
            }
        }
    }
    
    std::vector<unsigned> result;
    result.assign(saturated_edges.begin(), saturated_edges.end());
    return result;
}

// A function which returns the indices of the LSPs in the sense of https://arxiv.org/pdf/2401.13523
std::vector<unsigned> lspFind(){
    std::vector<unsigned> result;
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    LSP_STORE.clear();
    
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        unsigned count = 0;
        for(unsigned j=i; j<ALL_STORE.size()-1; ++j){
            if((isCompatible(ALL_STORE[i], ALL_STORE[j])) & (j != ALL_STORE.size()) & (saturatedHull(ALL_STORE[i]) != ALL_STORE[j])){
                count++;
            }
        }
        if(count == 0){
            result.push_back(i);
            LSP_STORE.push_back(ALL_STORE[i]);
        }
    }
    return result;
}


// A function which returns if a given transfer system is connected or not
bool isConnected(const std::vector<unsigned>& rhs){
    return(saturatedHull(rhs).size() == lattice.size());
}

// A function which returns the dual of a transfer system
// This will only work in the case that G=Cyclic group
// We use Theorem 4.21 of FOOQW to compute this as the dual of the opposite of the left set
std::vector<unsigned> dualTransferSystem(const std::vector<unsigned>& rhs){
    // If we are not cyclic then we return what we started with
    if(!(std::regex_search(subgroup_dictionary[subgroup_dictionary.size()-1], std::regex("^C\\d{1,}$")))){
        return rhs;
    }
    // Else we can compute the dual
    else{
        //We need to compute the involution by knowing the order of each subgroup
        std::vector<unsigned> subgroup_orders{1};
        for(unsigned i=1; i<subgroup_dictionary.size(); ++i){
            std::string subgroup = subgroup_dictionary[i];
            subgroup.erase(0, 1);
            subgroup_orders.push_back(std::stoi(subgroup));
        }
        unsigned group_order = subgroup_orders.back();
        auto left_set = leftSet(rhs);
        std::vector<unsigned> result;
        
        for(unsigned i=0; i<left_set.size(); ++i){
            unsigned subgroup1 = unsigned(group_order/subgroup_orders[lattice[left_set[i]].second]);
            unsigned subgroup2 = unsigned(group_order/subgroup_orders[lattice[left_set[i]].first]);
            unsigned subgroup1_index = unsigned(std::find(subgroup_orders.begin(), subgroup_orders.end(), subgroup1) - subgroup_orders.begin());
            unsigned subgroup2_index = unsigned(std::find(subgroup_orders.begin(), subgroup_orders.end(), subgroup2) - subgroup_orders.begin());
            std::pair<unsigned, unsigned> new_edge{subgroup1_index,subgroup2_index};
            result.push_back(unsigned(std::find(lattice.begin(), lattice.end(), new_edge) - lattice.begin()));
        }
        std::sort(result.begin(), result.end());
        return result;
    }
}

// The beginning of a function which will produce a numerical data sheet for a given group
void dataSheet(){
    std::string output;
    output += "G=" + subgroup_dictionary[subgroup_dictionary.size()-1] + "\n";
    
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    output += "#Transfer Systems=" + std::to_string(ALL_STORE.size()) + "\n";
    output += "Complexity=" + std::to_string(ALL_COMPLEXITY) + "\n";
    
    output += "Generation Statistics={1";
    
    for(unsigned i=0; i<GENERATION_STATISTICS.size()-1; ++i){
        output += "," + std::to_string(GENERATION_STATISTICS[i]);
    }
    
    output += "}\n";
    
    if(OPPOSITE_SATURATED_STORE.size() == 0){
        transferFind(false, SATURATED);
    }
    output += "#Saturated Transfer Systems=" + std::to_string(OPPOSITE_SATURATED_STORE.size()) + "\n";
    output += "Saturated Complexity=" + std::to_string(SATURATED_COMPLEXITY) + "\n";
    
    if(COSATURATED_STORE.size() == 0){
        transferFind(false, COSATURATED);
    }
    output += "#Cosaturated Transfer Systems=" + std::to_string(COSATURATED_STORE.size()) + "\n";
    output += "Cosaturated Complexity=" + std::to_string(COSATURATED_COMPLEXITY) + "\n";
    
    output += "Width=" + std::to_string(width()) + "\n";
    
    output += "#Flat transfers=" + std::to_string(flatTransfers().size()) + "\n";
    
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    
    output += "#Premodel structures=" + std::to_string(TRANSFER_LATTICE.size()) + "\n";
    
    if(QUILLEN_PAIRS.size() == 0){
        modelPairs();
    }
    
    output += "#Composition closed structures=" + std::to_string(CCLOSED_PAIRS.size()) + "\n";
    output += "#Quillen structures=" + std::to_string(QUILLEN_PAIRS.size()) + "\n";
    output += "#Weak equivalence types=" + std::to_string(weakEquivalenceTypes().size()) + "\n";
    
    output += "#Compatible pairs=" + std::to_string(compatiblePairs().size()) + "\n";
    
    std::cout << output << std::endl;
}

// A copy of the above function, but this time outputting to a LaTeX table (only suitable for small groups)
void dataSheetLatex(){
    std::string output;
    
    
    output += "\\begin{table}[]\n\\begin{tabular}{|cc|}\n\\hline\n\\multicolumn{2}{|c|}{$G = ";
    output += subgroup_dictionary[subgroup_dictionary.size()-1];
    output += "$} \\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Transfer systems} & ";
    
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    output += std::to_string(ALL_STORE.size());
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Complexity} & " + std::to_string(ALL_COMPLEXITY);
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Width} & " + std::to_string(width());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Generation values} & \\{1";
    for(unsigned i=0; i<GENERATION_STATISTICS.size()-1; ++i){
        output += "," + std::to_string(GENERATION_STATISTICS[i]);
    }
    output += "\\}";
    
    if(OPPOSITE_SATURATED_STORE.size() == 0){
        transferFind(false, SATURATED);
    }
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Saturated} & " + std::to_string(OPPOSITE_SATURATED_STORE.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Saturated complexity} & " + std::to_string(SATURATED_COMPLEXITY);
    
    if(COSATURATED_STORE.size() == 0){
        transferFind(false, COSATURATED);
    }
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Cosaturated} & " + std::to_string(COSATURATED_STORE.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Cosaturated complexity} & " + std::to_string(COSATURATED_COMPLEXITY);
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Flat} & " + std::to_string(flatTransfers().size());
    
    
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Premodel structures} & " + std::to_string(TRANSFER_LATTICE.size());
    
    if(QUILLEN_PAIRS.size() == 0){
        modelPairs();
    }
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#C.closed structures} & " + std::to_string(CCLOSED_PAIRS.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Quillen structures} & " + std::to_string(QUILLEN_PAIRS.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Weak equivalence types} & " + std::to_string(weakEquivalenceTypes().size());
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Compatible pairs} & " + std::to_string(compatiblePairs().size());
    
    output += " \\\\ \\hline\n\\end{tabular}\n\\end{table}";
    std::cout << output << std::endl;
}

// Another data sheet which only computes basic facts (no interval satatistics)
void dataSheetLatexRedux(){
    std::string output;
    
    
    output += "\\begin{table}[]\n\\begin{tabular}{|cc|}\n\\hline\n\\multicolumn{2}{|c|}{$G = ";
    output += subgroup_dictionary[subgroup_dictionary.size()-1];
    output += "$} \\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Transfer systems} & ";
    
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    output += std::to_string(ALL_STORE.size());
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Complexity} & " + std::to_string(ALL_COMPLEXITY);
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Width} & " + std::to_string(width());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Generation values} & \\{1";
    for(unsigned i=0; i<GENERATION_STATISTICS.size()-1; ++i){
        output += "," + std::to_string(GENERATION_STATISTICS[i]);
    }
    output += "\\}";
    
    if(OPPOSITE_SATURATED_STORE.size() == 0){
        transferFind(false, SATURATED);
    }
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Saturated} & " + std::to_string(OPPOSITE_SATURATED_STORE.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Saturated complexity} & " + std::to_string(SATURATED_COMPLEXITY);
    
    if(COSATURATED_STORE.size() == 0){
        transferFind(false, COSATURATED);
    }
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Cosaturated} & " + std::to_string(COSATURATED_STORE.size());
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{Cosaturated complexity} & " + std::to_string(COSATURATED_COMPLEXITY);
    
    output += "\\\\ \\hline\n\\multicolumn{1}{|c|}{\\#Flat} & " + std::to_string(flatTransfers().size());
    output += " \\\\ \\hline\n\\end{tabular}\n\\end{table}";
    std::cout << output << std::endl;
}

// A function which returns a Sage command for the poset of transfer systems
void printSageTransferPoset(){
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    std::string sage_string = "P = Poset({";
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        sage_string += std::to_string(i) + ":[";
        for(unsigned l=0; l<TRANSFER_LATTICE.size(); ++l){
            if(TRANSFER_LATTICE[l].first == i & TRANSFER_LATTICE[l].second != i){
                sage_string += std::to_string(TRANSFER_LATTICE[l].second) + ",";
            }
        }
        sage_string += "],";
    }
    sage_string += "})";
    std::cout << sage_string << std::endl;
}

// A function which returns a Sage command for the poset of transfer systems whose intervals are CClosed
void printSageCClosedPoset(){
    if(CCLOSED_PAIRS.size() == 0){
        modelPairs();
    }
    std::string sage_string = "P = Poset({";
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        sage_string += std::to_string(i) + ":[";
        for(unsigned l=0; l<CCLOSED_PAIRS.size(); ++l){
            if(TRANSFER_LATTICE[CCLOSED_PAIRS[l]].first == i & TRANSFER_LATTICE[CCLOSED_PAIRS[l]].second != i){
                sage_string += std::to_string(TRANSFER_LATTICE[CCLOSED_PAIRS[l]].second) + ",";
            }
        }
        sage_string += "],";
    }
    sage_string += "})";
    std::cout << sage_string << std::endl;
}

// A function which returns a Sage command for the poset of transfer systems whose intervals are Quillen
void printSageQuillenPoset(){
    if(QUILLEN_PAIRS.size() == 0){
        modelPairs();
    }
    std::string sage_string = "P = Poset({";
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        sage_string += std::to_string(i) + ":[";
        for(unsigned l=0; l<QUILLEN_PAIRS.size(); ++l){
            if(TRANSFER_LATTICE[QUILLEN_PAIRS[l]].first == i & TRANSFER_LATTICE[QUILLEN_PAIRS[l]].second != i){
                sage_string += std::to_string(TRANSFER_LATTICE[QUILLEN_PAIRS[l]].second) + ",";
            }
        }
        sage_string += "],";
    }
    sage_string += "})";
    std::cout << sage_string << std::endl;
}

// A function which returns the subgroup dictionary
std::string subgroupDictionary(){
    std::string result;
    unsigned group_counter = 0;
    unsigned group_counter2 = 0;
    for(unsigned i = 0; i < subgroup_dictionary.size(); ++i){
        if( i > 0 && subgroup_dictionary[i] == subgroup_dictionary[i+1]){
            group_counter += 1;
        }
        else{
            group_counter2 = group_counter;
            group_counter = 0;
        }
        result += "{" + std::to_string(i) + ":" + subgroup_dictionary[i];
        if(group_counter > 0){
            result += "(" + std::to_string(group_counter) + ")";
        }
        if(group_counter2 > 0){
            result += "(" + std::to_string(group_counter2+1) + ")";
            group_counter2 = 0;
        }
        result += "}\n";
    }
    
    std::string conj_string;
    
    
    CONJUGACY_CLASSES.clear();
    conj_string += "\nConjugacy Classes:\n[0]\n";
    std::vector<unsigned> seen_subgroups{0};
    std::vector<unsigned> conj_class{0};
    CONJUGACY_CLASSES.push_back(conj_class);
    for(unsigned i = 1; i<subgroup_dictionary.size(); ++i){
        if(std::find(seen_subgroups.begin(), seen_subgroups.end(), i) == seen_subgroups.end()){
            conj_class.clear();
            std::string conjs = "[";
            std::pair<unsigned, unsigned> lattice_edge{0, i};
            unsigned index_of_edge = unsigned(std::find(lattice.begin(), lattice.end(), lattice_edge) - lattice.begin());
            std::set<unsigned> conj_store;
            for(unsigned p=0; p<conjugates[index_of_edge].size(); ++p){
                seen_subgroups.push_back(lattice[conjugates[index_of_edge][p]].second);
                conj_store.insert(seen_subgroups.back());
            }
            for(auto it = conj_store.begin(); it != conj_store.end(); ++it){
                conjs += std::to_string(*it) + ",";
                conj_class.push_back(*it);
            }
            CONJUGACY_CLASSES.push_back(conj_class);
            conjs.pop_back();
            conjs += "]\n";
            conj_string += conjs;
        }
        
    }
    
    //If we are not Dedkind then we need to give information about conjugates
    if(!isDedekind()){
        result += conj_string;
    }
    
    return result;
}

// A function which will return the data of a given transfer system
std::string stringTransferSystem(const std::vector<unsigned>& rhs) {
    std::string result;
    if(rhs.size() == 0){
        return result += "(Ø)";
    }
    result += "{[";
    for(unsigned i = 0; i<rhs.size(); ++i){
        result += std::to_string(lattice[rhs[i]].first) + "," + std::to_string(lattice[rhs[i]].second) + "],[";
    }
    result.pop_back();
    result.pop_back();
    result += "}";
    return result;
}

// A collection of user-friendly helper functions to simply access numerical results
void printNumberOfTransfers(){
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    std::cout << ALL_STORE.size() << std::endl;
}

void printNumberOfSaturatedTransfers(){
    if(OPPOSITE_SATURATED_STORE.size() == 0){
        transferFind(false, SATURATED);
    }
    std::cout << OPPOSITE_SATURATED_STORE.size() << std::endl;
}

void printNumberOfCosaturatedTransfers(){
    if(COSATURATED_STORE.size() == 0){
        transferFind(false, COSATURATED);
    }
    std::cout << COSATURATED_STORE.size() << std::endl;
}

void printNumberOfUnderlyingTransfers(){
    if(UNDERLYING_STORE.size() == 0){
        transferFind(false, UNDERLYING);
    }
    std::cout << UNDERLYING_STORE.size() << std::endl;
}

void printGenerationStatistics(){
    if(GENERATION_STATISTICS.size() == 0){
        transferFind(false, ALL);
    }
    std::string output = "1";
    for(unsigned i=0; i<GENERATION_STATISTICS.size()-1; ++i){
        output += "," + std::to_string(GENERATION_STATISTICS[i]);
    }
    std::cout << output << std::endl;
}

void printNumberOfFlatTransfers(){
    std::cout << flatTransfers().size() << std::endl;
}

void printWidth(){
    std::cout << width() << std::endl;
}

void printComplexity(){
    if(ALL_COMPLEXITY == 0){
        transferFind(false, ALL);
    }
    std::cout << ALL_COMPLEXITY << std::endl;
}

void printNumberOfMaximallyGenerated(){
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    std::cout << MAXIMALLY_GENERATED.size() << std::endl;
}

void printNumberOfTransferPairs(){
    if(TRANSFER_LATTICE.size() == 0){
        transferLattice();
    }
    std::cout << TRANSFER_LATTICE.size() << std::endl;
}

void printNumberOfLSPs(){
    std::cout << lspFind().size() << std::endl;
}

void printNumberOfCompatiblePairs(){
    if(COMPATIBLE_PAIRS.size() == 0){
        compatiblePairs();
    }
    std::cout << COMPATIBLE_PAIRS.size() << std::endl;
}

void printNumberOfCClosedPairs(){
    if(CCLOSED_PAIRS.size() == 0){
        modelPairs();
    }
    std::cout << CCLOSED_PAIRS.size() << std::endl;
}

void printNumberOfQuillenPairs(){
    if(QUILLEN_PAIRS.size() == 0){
        modelPairs();
    }
    std::cout << QUILLEN_PAIRS.size() << std::endl;
}
void printNumberOfWeakEquivalenceTypes(){
    std::cout << weakEquivalenceTypes().size() << std::endl;
}

void printSubgroupDictionary(){
    std::cout << subgroupDictionary() << std::endl;
}

void printTransferSystem(const std::vector<unsigned>& rhs){
    std::cout << stringTransferSystem(rhs) << std::endl;
}

void printAllTransfers(){
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    for(unsigned i = 0; i < ALL_STORE.size(); ++i){
        std::cout << stringTransferSystem(ALL_STORE[i]) << std::endl;
    }
}

// A function which takes in a collection of edges and produces a TikZ string to plot it
// This will work with the poset up to conjugacy and will auto assign coords in a circle
// Unless explicit coordinates and labels are provided in the header file
void edgesToTikz(const std::vector<unsigned>& rhs){
    // Populate the conjugacy classes
    if(CONJUGACY_CLASSES.size() == 0){
        subgroupDictionary();
    }
    std::string output = "\\begin{tikzpicture}\n";
    unsigned num_nodes = unsigned(CONJUGACY_CLASSES.size());
    
    if(vertex_layout.size() != num_nodes){
        if(pretty_subgroup_dictionary.size() == num_nodes){
            for(unsigned i=0; i<num_nodes; ++i){
                double theta = 2.0*double(i)*3.14159/double(num_nodes);
                output += "\\node[inner sep=0cm] (" + std::to_string(i) + ") at (" + std::to_string(3.0*sin(theta)).substr(0,6) + "," + std::to_string(3.0*cos(theta)).substr(0,6) + ") {$" + pretty_subgroup_dictionary[i] + "$};\n";
            }
        }
        else{
            for(unsigned i=0; i<num_nodes; ++i){
                double theta = 2.0*double(i)*3.14159/double(num_nodes);
                output += "\\node[inner sep=0cm] (" + std::to_string(i) + ") at (" + std::to_string(3.0*sin(theta)).substr(0,6) + "," + std::to_string(3.0*cos(theta)).substr(0,6) + ") {$" + subgroup_dictionary[CONJUGACY_CLASSES[i][0]] + "$};\n";
            }
        }
    }
    else{
        if(pretty_subgroup_dictionary.size() == num_nodes){
            for(unsigned i=0; i<num_nodes; ++i){
                output += "\\node[inner sep=0cm] (" + std::to_string(i) + ") at " + vertex_layout[i] + "{$" + pretty_subgroup_dictionary[i] + "$};\n";
            }
        }
        else{
            for(unsigned i=0; i<num_nodes; ++i){
                output += "\\node[inner sep=0cm] (" + std::to_string(i) + ") at " + vertex_layout[i] + "{$" + subgroup_dictionary[CONJUGACY_CLASSES[i][0]] + "$};\n";
            }
        }
    }
    
    std::vector<std::string> transfer_edges;
    for(unsigned i=0; i<num_nodes; ++i){
        std::pair<unsigned,unsigned> test_inclusion{CONJUGACY_CLASSES[i][0],CONJUGACY_CLASSES[i][0]};
        for(unsigned j=0; j<num_nodes; ++j){
            for(unsigned p=0; p<CONJUGACY_CLASSES[j].size(); ++p){
                test_inclusion.second = CONJUGACY_CLASSES[j][p];
                if(std::find(lattice.begin(), lattice.end(), test_inclusion) != lattice.end()){
                    unsigned index = unsigned(std::find(lattice.begin(), lattice.end(), test_inclusion) - lattice.begin());
                    if(std::find(rhs.begin(), rhs.end(), index) != rhs.end()){
                        if(edge_options.size() == num_nodes){
                            transfer_edges.push_back("\\draw[red,->] (" + std::to_string(i) + ") edge" + edge_options[i][j] + " (" + std::to_string(j) + ");\n");
                        }
                        else{
                            transfer_edges.push_back("\\draw[red,->] (" + std::to_string(i) + ") edge (" + std::to_string(j) + ");\n");
                        }
                    }
                    else{
                        if(edge_options.size() == num_nodes){
                            output += "\\draw[black!10,->] (" + std::to_string(i) + ") edge" + edge_options[i][j] + " (" + std::to_string(j) + ");\n";
                        }
                        else{
                            output += "\\draw[black!10,->] (" + std::to_string(i) + ") edge (" + std::to_string(j) + ");\n";
                        }
                    }
                    break;
                }
            }
        }
    }
    
    for(unsigned i=0; i<transfer_edges.size(); ++i){
        output += transfer_edges[i];
    }
    
    output += "\\end{tikzpicture}";
    std::cout << output << std::endl;
    
}

// A function which returns the poset Sub(G)/G
std::vector<std::pair<unsigned, unsigned>> latticeUpToConjugacy(){
    std::vector<std::pair<unsigned, unsigned>> result;
    if(CONJUGACY_CLASSES.size() == 0){
        subgroupDictionary();
    }
    
    std::pair<unsigned, unsigned> new_edge;
    for(unsigned i=0; i<CONJUGACY_CLASSES.size(); ++i){
        for(unsigned j=0; j<CONJUGACY_CLASSES.size(); ++j){
            new_edge.first = i;
            new_edge.second = j;
            bool found = false;
            for(unsigned a = 0; a < CONJUGACY_CLASSES[i].size(); ++a){
                for(unsigned b = 0; b < CONJUGACY_CLASSES[j].size(); ++b){
                    std::pair<unsigned,unsigned> to_search{CONJUGACY_CLASSES[i][a], CONJUGACY_CLASSES[j][b]};
                    if(std::find(lattice.begin(), lattice.end(), to_search) != lattice.end() & !found){
                        result.push_back(new_edge);
                        found = true;
                        break;
                    }
                    
                }
            }
        }
    }
    
    return result;
}

// A function to find the maximal w in the intersection of the downset of x with the downset of z
std::vector<unsigned> maxInIntersection(const unsigned& x, const unsigned& z){
    std::vector<unsigned> result;
    auto lattice_up_to_conjugacy = latticeUpToConjugacy();
    
    if(x == z){
        result.push_back(x);
        return result;
    }
    
    // Start by finding the downset
    std::vector<unsigned> below_x{x};
    std::vector<unsigned> below_z{z};
    
    for(unsigned i=0; i<lattice_up_to_conjugacy.size(); ++i){
        if(lattice_up_to_conjugacy[i].second == x){
            below_x.push_back(lattice_up_to_conjugacy[i].first);
        }
        if(lattice_up_to_conjugacy[i].second == z){
            below_z.push_back(lattice_up_to_conjugacy[i].first);
        }
    }

    
    std::sort(below_x.begin(), below_x.end());
    std::sort(below_z.begin(), below_z.end());
    

    std::vector<unsigned> in_intersection;
    std::set_intersection(below_x.begin(),below_x.end(),below_z.begin(),below_z.end(),back_inserter(in_intersection));
    
    // Now find the maximal elements of in_intersection;
    for(unsigned m = 0; m < in_intersection.size(); ++m){
        bool m_is_maximal = true;
        for(unsigned a = 0; a < in_intersection.size(); ++ a){
            std::pair<unsigned,unsigned> test_edge{in_intersection[m],in_intersection[a]};
            if(std::find(lattice_up_to_conjugacy.begin(), lattice_up_to_conjugacy.end(), test_edge) != lattice_up_to_conjugacy.end()){
                m_is_maximal = false;
            }
        }
        
        if(m_is_maximal){
            result.push_back(in_intersection[m]);
        }
    }
    
    return result;
}

// A function which computes all the maximal w as above in Sub(G)/G
std::vector<std::vector<std::vector<unsigned>>> conjugationMaximalElements(){
    if(CONJUGACY_CLASSES.size() == 0){
        subgroupDictionary();
    }
    std::vector<std::vector<std::vector<unsigned>>> result;
    for(unsigned i=0; i<CONJUGACY_CLASSES.size(); ++i){
        std::vector<std::vector<unsigned>> row;
        for(unsigned j=0; j<CONJUGACY_CLASSES.size(); ++j){
            row.push_back(maxInIntersection(i, j));
        }
        result.push_back(row);
    }
    return result;
}

// A function which computes the intersections in Sub(G)/G
std::vector<std::vector<unsigned>> intersectionsUpToConjugacy(){
    std::vector<std::vector<unsigned>>  result;
    
    auto lattice_up_to_conjugacy = latticeUpToConjugacy();
    auto conjugation_maximal_elements = conjugationMaximalElements();
    
    for(unsigned i=0; i<lattice_up_to_conjugacy.size(); ++i){
        std::vector<unsigned> row{i};
        // So we have some (x,y) and we want to get all (w,z) for all z < y and w maximal in x \cap z.
        unsigned x = lattice_up_to_conjugacy[i].first;
        unsigned y = lattice_up_to_conjugacy[i].second;
        
        // Browse though all z < y
        std::vector<unsigned> z_array;
        for(unsigned z = 0; z < CONJUGACY_CLASSES.size(); ++z){
            std::pair<unsigned,unsigned> test_edge{z,y};
            if(std::find(lattice_up_to_conjugacy.begin(), lattice_up_to_conjugacy.end(), test_edge) != lattice_up_to_conjugacy.end()){
                z_array.push_back(z);
            }
        }
        
        for(unsigned k=0; k<z_array.size(); ++k){
            auto maxEls = conjugation_maximal_elements[z_array[k]][x];
            for(unsigned w=0; w<maxEls.size(); ++w){
                //We now find the (w,z)
                if(maxEls[w] != z_array[k]){
                    std::pair<unsigned,unsigned> new_edge{maxEls[w], z_array[k]};
                    //Now we find the index of new_edge in the lattice and we append it
                    unsigned index = unsigned(std::find(lattice_up_to_conjugacy.begin(), lattice_up_to_conjugacy.end(), new_edge) - lattice_up_to_conjugacy.begin());
                    row.push_back(index);
                }
            }
        }
        result.push_back(row);
    }
    return result;
}

// A function which computes the transitivity matrix for Sub(G)/G
std::vector<std::vector<unsigned>> findTransitiveClosureConjugation(){
    std::vector<std::vector<unsigned>> transitive_closure_conjugation;
    auto lattice_up_to_conjugacy = latticeUpToConjugacy();
    for(unsigned i=0; i<lattice_up_to_conjugacy.size(); ++i){
        std::vector<unsigned> i_row;
        for(unsigned j=0; j<lattice_up_to_conjugacy.size(); ++j){
            unsigned pos = std::min(i,j);
            if(lattice_up_to_conjugacy[i].second == lattice_up_to_conjugacy[j].first & i != j){
                std::pair<unsigned,unsigned> test_el{lattice_up_to_conjugacy[i].first, lattice_up_to_conjugacy[j].second};
                auto it = std::find(lattice_up_to_conjugacy.begin(), lattice_up_to_conjugacy.end(), test_el);
                if(it != lattice_up_to_conjugacy.end()){
                    pos = unsigned(std::distance(lattice_up_to_conjugacy.begin(), it));
                }
            }
            i_row.push_back(pos);
        }
        transitive_closure_conjugation.push_back(i_row);
    }
    return transitive_closure_conjugation;
}

// This function which returns all the transfer systems on Sub(G)/G
std::vector<std::vector<unsigned>> conjugacyTransfers(){
    if(CONJUGACY_STORE.size() == 0){
        
        auto res = intersectionsUpToConjugacy();
        auto new_lattice = latticeUpToConjugacy();
        auto new_intersections = intersectionsUpToConjugacy();
        auto new_transitive_closure = findTransitiveClosureConjugation();
        
        auto old_lattice = lattice;
        auto old_intersections = intersections;
        auto old_transitive_closure = transitive_closure;
        
        lattice = new_lattice;
        intersections = new_intersections;
        transitive_closure = new_transitive_closure;
        
        transferFind(false, CONJUGACY);
        
        //Restore the old data
        lattice = old_lattice;
        intersections = old_intersections;
        transitive_closure = old_transitive_closure;
        
    }
    return CONJUGACY_STORE;
}

void printNumberOfConjugacyTransfers(){
    if(CONJUGACY_STORE.size() == 0){
        conjugacyTransfers();
    }
    std::cout << CONJUGACY_STORE.size() << std::endl;
}
#endif /* ninfty_h */
