//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/Cp2q.h"
#include "ninfty.h"


int main() {
//    transferFind(false, ALL);
    dataSheet();
    
    unsigned counter = 0;
    std::vector<unsigned> cc_store;
    for(unsigned i=0; i<TRANSFER_LATTICE.size(); ++i){
        if(modelCheck(ALL_STORE[TRANSFER_LATTICE[i].first], leftSet(ALL_STORE[TRANSFER_LATTICE[i].second])) >= 1){
            counter++;
            cc_store.push_back(i);
        }
    }

    
    std::string sage_string = "P = Poset({";
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        sage_string += std::to_string(i) + ":[";
        for(unsigned l=0; l<cc_store.size(); ++l){
            if(TRANSFER_LATTICE[cc_store[l]].first == i & TRANSFER_LATTICE[cc_store[l]].second != i){
                sage_string += std::to_string(TRANSFER_LATTICE[cc_store[l]].second) + ",";
            }
        }
        sage_string += "],";
    }
    sage_string += "})";
    std::cout << sage_string << std::endl;
    
    //sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
    std::cout << counter << std::endl;
    
//    std::cout << ALL_STORE[37].size() << std::endl;
    for(unsigned i=0; i<ALL_STORE[65].size(); ++i){
        std::cout << subgroup_dictionary[lattice[ALL_STORE[65][i]].first] << " --> " << subgroup_dictionary[lattice[ALL_STORE[65][i]].second] << std::endl;
    }
    
    std::cout << ALL_STORE[65].size() << std::endl;

}
