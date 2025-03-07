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
    std::cout << subgroupDictionary() << std::endl;
    return 0;
}
