#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C3",
"C9",
"C27",
"C81"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,3},
{0,4},
{1,2},
{1,3},
{1,4},
{2,3},
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
{6},
{7},
{8},
{9}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{0},
{0,1},
{0,1,2},
{},
{4},
{4,5},
{},
{7},
{}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{},
{4},
{5,7},
{6,8,9},
{},
{7},
{8,9},
{},
{9},
{}
}; 
 

std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};