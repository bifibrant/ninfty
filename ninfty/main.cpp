//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/Cpq.h"
#include "ninfty.h"


int main() {
    transferFind(false, ALL);
    auto lattice_storage = transferLattice();
    unsigned counter = 0;
    

    
    for(unsigned i=0; i<lattice_storage.size(); ++i){
        //std::cout << lattice_storage[i].first << " " << lattice_storage[i].second <<  std::endl;
        if(isCompatible(lattice_storage[i])){
            counter++;
            //std::cout << lattice_storage[i].first << " " << lattice_storage[i].second << std::endl;
        }

    }
    std::cout << counter << std::endl;
}
