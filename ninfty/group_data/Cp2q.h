#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C3",
"C4",
"C6",
"C12"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,2},
{0,1},
{0,4},
{0,3},
{0,5},
{2,4},
{2,5},
{1,4},
{1,3},
{1,5},
{4,5},
{3,5}
}; 
 

std::vector<std::vector<unsigned>> conjugates{
{0},
{1},
{2},
{3},
{4},
{5},
{6},
{7},
{8},
{9},
{10},
{11}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{0,1},
{1},
{0,1,2,3},
{1},
{1,5,3},
{0},
{},
{0,7,8},
{8},
{0,7}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{7,11},
{5},
{5,7,11},
{6,8,10},
{6,9,10,11},
{},
{10},
{11},
{10},
{10,11},
{},
{}
}; 
 
