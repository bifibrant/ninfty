#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C4"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{1,2}
}; 
 

std::vector<std::vector<unsigned>> conjugates{
{0},
{1},
{2}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{0},
{}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{},
{2},
{}
}; 
 
