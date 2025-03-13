#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C4",
"C8"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,3},
{1,2},
{1,3},
{2,3}
}; 
 

std::vector<std::vector<unsigned>> conjugates{
{0},
{1},
{2},
{3},
{4},
{5}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{0},
{0,1},
{},
{3},
{}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{},
{3},
{4,5},
{},
{5},
{}
}; 
 
std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};
