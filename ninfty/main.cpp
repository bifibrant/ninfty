//
//  main.cpp
//  ninfty
//
//  Created by Scott Balchin on 19/02/2025.
//

#include "group_data/S3.h"
#include "ninfty.h"


int main() {
//    dataSheetLatexRedux();
    transferFind(false, ALL);
    transferFind(false, UNDERLYING);
    
    std::cout << ALL_STORE.size() << std::endl;
    std::cout << UNDERLYING_STORE.size() << std::endl;
    
    return 0;
}
