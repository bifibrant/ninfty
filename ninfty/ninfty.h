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

// Variable which stores all requested transfer systems
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> RESULT;
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> ALL_STORE;
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> SATURATED_STORE;
std::unordered_set<std::vector<unsigned>, unsigned_vector_hasher> COSATURATED_STORE;

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
    for(unsigned i=0; i<lattice.size(); ++i){
        if(lattice[i].second == subgroup_dictionary.size()-1){
            edgesToG.push_back(i);
        }
    }
}

// An algorithm which finds (and updates) the edges which are used for saturation generation
void findSaturationEdges(){
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
    
    unsigned long old_size = 0;
    unsigned long diff = 1;
    unsigned gen_step = 1;
    
    // Initalise the algorithm with the empty transfer system
    std::vector<unsigned> empty_transfer{};
    NEW_EDGES.push_back(empty_transfer);
    RESULT.insert(empty_transfer);
    GOOD_EDGES.clear();
    
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
        gen_step++;
    }
    // Store the generation step as the complexity
    if(gen_type == ALL){
        ALL_COMPLEXITY = gen_step-1;
        ALL_STORE = RESULT;
    }
    else if(gen_type == SATURATED){
        SATURATED_COMPLEXITY = gen_step-1;
        SATURATED_STORE = RESULT;
    }
    else if(gen_type == COSATURATED){
        COSATURATED_COMPLEXITY = gen_step-1;
        COSATURATED_STORE = RESULT;
    }
    RESULT.clear();
    NEW_EDGES.clear();
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

// Implementation ToDo
// Compatible pairs:
    // Compatibility check
    // Intervals of xfer systems
    // All compatible intervals
// Model structures:
    // Intervals of xfer systems
    // Left set
    // Extension + complements
    // Weak equivalences
    // Model check
    // Composition closed check
// bool isSaturated
// bool isCosaturated
// return maximally generated things
// Involution of transfer systems (only in the cyclic case)
// (co)saturated hull
// Minimal "fibrant" node
// Flat transfer systems


#endif /* ninfty_h */
