#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C3",
"C9",
"C27",
"C81",
"C243"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,3},
{0,4},
{0,5},
{1,2},
{1,3},
{1,4},
{1,5},
{2,3},
{2,4},
{2,5},
{3,4},
{3,5},
{4,5}
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
{11},
{12},
{13},
{14}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{0},
{0,1},
{0,1,2},
{0,1,2,3},
{},
{5},
{5,6},
{5,6,7},
{},
{9},
{9,10},
{},
{12},
{}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{},
{5},
{6,9},
{7,10,12},
{8,11,13,14},
{},
{9},
{10,12},
{11,13,14},
{},
{12},
{13,14},
{},
{14},
{}
}; 
 

std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};