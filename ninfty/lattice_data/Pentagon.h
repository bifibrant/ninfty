#include <vector>

std::vector<std::string> subgroup_dictionary{
"0",
"1",
"2",
"3",
"4"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,4},
{0,3},
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
{7}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{0,1,3},
{1},
{1,3},
{},
{0,5},
{0}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{6,7},
{4},
{4,6,7},
{4,5},
{},
{},
{7},
{}
}; 
 

std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};