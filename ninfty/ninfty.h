//
//  ninfty.h
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#ifndef ninfty_h
#define ninfty_h

#include <iostream>
#include <vector>
#include <set>
#include <math.h>
#include <chrono>
#include <unordered_set>
#include <future>
#include <random>

enum GenerationType{ ALL, SATURATED, COSATURATED };

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

// Variable which stores the array of meets
std::vector<std::vector<unsigned>> meetArray;

// Variables which stores all requested transfer systems
// Note that OPPOSITE_SATURATED_STORE stores the "cosatuated transfer systems" on the opposite lattice and does not actually give the saturated transfer systems
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> RESULT;
std::vector<std::vector<unsigned>> ALL_STORE;
std::vector<std::vector<unsigned>> SATURATED_STORE;
std::vector<std::vector<unsigned>> OPPOSITE_SATURATED_STORE;
std::vector<std::vector<unsigned>> COSATURATED_STORE;

// A varaible which stores the maximally generated transfer systems
std::vector<std::vector<unsigned>> MAXIMALLY_GENERATED;

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

// Setting up global variables and constants
unsigned NUM_THREADS = unsigned(lattice.size());
std::vector<std::vector<unsigned>> NEW_EDGES;
std::vector<std::shared_future<std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher>>> THREAD_STORE;
std::vector<unsigned> GOOD_EDGES;
std::vector<std::vector<unsigned>> transitive_closure;
std::vector<unsigned> edgesFromE;
std::vector<unsigned> edgesToG;

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

// An algorithm which finds (and updates) the transitive closre of the lattice elements
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

// An algorithm to comptue the meet (i.e., intersection) of two subgroups A and B
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
// The enum is used to see if we want to use intersection (for normal generation) or cointersections (for saturated generation)
std::vector<unsigned> transferClosure(const std::vector<unsigned>& r0, const GenerationType& gen_type = ALL){
    std::set<unsigned> r1_set, r2_set, r4_set;

    // We check to see if we have computed the transitive closure or not
    if(transitive_closure.size() == 0){
        findTransitiveClosure();
    }
        
    // Close under conjugation
    for(auto it = r0.begin(); it != r0.end(); ++it){
        r1_set.insert(conjugates[(*it)].begin(), conjugates[(*it)].end());
//                This will allow for the underlying lattice calculation
//                r1_set.insert(*it);
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
// Algorithm adapted from CITE
// The verbose toggle can be used to supress generation statistics
// The gen_type adjusts if it is all, saturated, or cosaturated as required
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
    // Note that this only has an effect when the group G is non-abelian
    std::vector<std::vector<unsigned>> temp_new_edges;
    
    auto it = NEW_EDGES.begin();
    for(unsigned i=0; i<lattice.size(); ++i){
        if(gen_type == ALL){
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
        // The choice of 100 here is arbitrary, but should work given all seen examples
        if(NEW_EDGES.size() < 100 & NEW_EDGES.size() != 0){
            MAXIMALLY_GENERATED = NEW_EDGES;
        }
        
    }
    // Store the generation step as the complexity, and sort the arrays for ease of browsing later
    if(gen_type == ALL){
        ALL_COMPLEXITY = gen_step-1;
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(ALL_STORE));
        std::sort(ALL_STORE.begin(), ALL_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == SATURATED){
        SATURATED_COMPLEXITY = gen_step-1;
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(OPPOSITE_SATURATED_STORE));
        std::sort(OPPOSITE_SATURATED_STORE.begin(), OPPOSITE_SATURATED_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
    }
    else if(gen_type == COSATURATED){
        COSATURATED_COMPLEXITY = gen_step-1;
        std::copy(RESULT.begin(), RESULT.end(), std::back_inserter(COSATURATED_STORE));
        std::sort(COSATURATED_STORE.begin(), COSATURATED_STORE.end(), [](const std::vector<unsigned> & a, const std::vector<unsigned> & b){ return a.size() < b.size(); });
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

// An implementation of Rubin's algorithm in reverse
// This provides a set of minimal generators for a given transfer system
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

bool isCosaturated(const std::vector<unsigned>& rhs){
    std::vector<unsigned> test_basis;
    for(unsigned i=0; i<rhs.size(); ++i){
        if(lattice[rhs[i]].second == subgroup_dictionary.size()-1){
            test_basis.push_back(rhs[i]);
        }
    }
    return (transferClosure(test_basis, COSATURATED) == rhs);
}

// A thread function that checks inclusions
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

// A function which determines if a pair of transfer systems is compatible in the sense of CITE
bool isCompatible(const std::pair<unsigned,unsigned> rhs){
    
    // rhs.first = multiplicative
    // rhs.second = additive
    auto transfer_m = ALL_STORE[rhs.first];
    auto transfer_a = ALL_STORE[rhs.second];
    
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
        if(isCompatible(TRANSFER_LATTICE[i])){
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
        // Note that this element is unique so this is well defined
        if((lattice[rhs[i]].second == subgroup_dictionary.size() - 1) & (lattice[rhs[i]].first < min_fibrant)){
            min_fibrant = lattice[rhs[i]].first;
        }
    }
    return min_fibrant;
}

// An algorithm which determines if a given transfer system is flat or not in the sense of CITE
// This asks that if F is the minimal fibrant node then the transfer system resticted to F is trivial.
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
    // As F is necessarily normal this is well-founded
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
    std::vector<unsigned> result;
    if(ALL_STORE.size() == 0){
        transferFind(false, ALL);
    }
    
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        if(isFlat(ALL_STORE[i])){
            result.push_back(i);
        }
        
    }
    return result;
}

/* Start of Model Structure functions */
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
                    //std::cout << "adding edge " << (it2 - lattice.begin()) << std::endl;
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
std::vector<unsigned> weakEquivalences(const std::vector<unsigned>& AF, const std::vector<unsigned>& AC){
    std::set<unsigned> store;
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
unsigned modelCheck(const std::vector<unsigned>& AF, const std::vector<unsigned>& AC){
    auto W = weakEquivalences(AF, AC);
    
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
        // Need to make sure that g<f works here as opposed to g<lattice.size()
        for(unsigned g=0; g<f; ++g){
            if(lattice[f].first == lattice[g].second){
                unsigned gf = transitive_closure[g][f];
                //If exacly 2 of them are in W then we fail.
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
        unsigned model_check = modelCheck(ALL_STORE[TRANSFER_LATTICE[i].first], leftSet(ALL_STORE[TRANSFER_LATTICE[i].second]));
        
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
    // (note if we have computed all of the saturated things there may be a more efficient algorithm TODO)
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

// The beginning of a function which will produce a numerical data sheet for a given group (eventually to be made into a LaTeX table
// This may run very slowly for large/compelx groups!
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
    output += "#Saturated Transfer Systems=" + std::to_string(COSATURATED_STORE.size()) + "\n";
    output += "Saturated Complexity=" + std::to_string(COSATURATED_COMPLEXITY) + "\n";
    
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
    output += "#Weak equivalence types" + std::to_string(weakEquivalenceTypes().size()) + "\n";
    
    output += "#Compatible pairs=" + std::to_string(compatiblePairs().size()) + "\n";
    
    std::cout << output << std::endl;
}

// A copy of the above function, but this time outputting to a LaTeX table
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


// Implementation ToDo
// General functions for people to access results
// Non equivariant (ie underlying) transfer systems and also maybe quotient poset
// return maximally generated things
// Involution of transfer systems (only in the cyclic case)
// Saving data from the code
// TikZ diagrams from the code (maybe just the generating sets? would require an input of the lattice of subgroups)

#endif /* ninfty_h */
