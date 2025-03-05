//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/Cqprs.h"
#include "ninfty.h"


int main() {
    transferFind(false, ALL);
    
    unsigned counter = 0;
    
    for(unsigned i=0; i<ALL_STORE.size(); ++i){
        if(isFlat(ALL_STORE[i])){
            counter++;
        }
    }
    
    std::cout << counter << " Flat transfer systems" << std::endl;
}
