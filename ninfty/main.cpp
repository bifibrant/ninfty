//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/Q8.h"
#include "ninfty.h"


int main() {
    //std::cout << sageCClosedPoset() << std::endl;
    transferFind(false, ALL);
    subgroupDictionary();
    std::cout << CONJUGACY_CLASSES.size() << std::endl;
    edgesToTikz(ALL_STORE[57]);
    edgesToTikz(ALL_STORE[62]);
    //printSubgroupDictionary();
    //    transferFind(true, ALL);
    return 0;
}
