//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/S4.h"
#include "ninfty.h"


int main() {
    transferFind(false, ALL);
    
    auto lattice_storage = transferLattice();
    unsigned counter = 0;
    

    std::cout << lattice_storage.size() << std::endl;
    std::cout << compatiblePairs().size() << std::endl;
    //
//    for(unsigned i=0; i<lattice_storage.size(); ++i){
//        if(isCompatible(lattice_storage[i])){
//            counter++;
//        }
//
//    }
//    std::cout << counter << std::endl;
}
