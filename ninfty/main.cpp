//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/Cp3q.h"
#include "ninfty.h"


int main() {
//    transferFind(false, ALL);

//    std::cout << weakEquivalenceTypes().size() << std::endl;
    
    dataSheetLatex();
    
    return 0;
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


}
