#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C2",
"C2",
"C2 x C2",
"C4",
"C4",
"C4 x C2"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,3},
{0,5},
{0,6},
{0,4},
{0,7},
{1,5},
{1,6},
{1,4},
{1,7},
{5,7},
{6,7},
{2,4},
{2,7},
{3,4},
{3,7},
{4,7}
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
{14},
{15},
{16},
{17}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{},
{0},
{0},
{0,1,2},
{0,3,4,1,2,5},
{},
{},
{1,2},
{7,8,1,2,9},
{8,1,2,9},
{7,1,2,9},
{0,2},
{0,3,4,2,13},
{0,1},
{0,3,4,1,15},
{7,8}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{13,15},
{9,11,12,15},
{9,11,12,13},
{7,12,14,16,17},
{8,11,14,16,17},
{9,11,12,13,15},
{10,11,12,14,16,17},
{12,17},
{11,17},
{11,12},
{11,12,17},
{},
{},
{},
{17},
{},
{17},
{}
}; 
 

std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};