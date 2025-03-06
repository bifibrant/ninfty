#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C2",
"C2",
"C2 x C2"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,3},
{0,4},
{1,4},
{2,4},
{3,4}
}; 
 

std::vector<std::vector<unsigned>> conjugates{
{0},
{1},
{2},
{3},
{4},
{5},
{6}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{},
{0,1,2},
{1,2},
{0,2},
{0,1}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{5,6},
{4,6},
{4,5},
{4,5,6},
{},
{},
{}
}; 
 
