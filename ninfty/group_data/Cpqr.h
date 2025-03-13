#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C3",
"C5",
"C6",
"C10",
"C15",
"C30"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,3},
{0,2},
{0,1},
{0,6},
{0,5},
{0,4},
{0,7},
{3,6},
{3,5},
{3,7},
{2,6},
{2,4},
{2,7},
{6,7},
{1,5},
{1,4},
{1,7},
{5,7},
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
{17},
{18}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{},
{0,1},
{0,2},
{1,2},
{0,1,3,2,4,5},
{1},
{2},
{1,7,2,8,5},
{0},
{2},
{0,10,2,4,11},
{2,8,11},
{0},
{1},
{0,1,3,14,15},
{1,7,15},
{0,10,14}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{10,14,18},
{7,15,17},
{8,11,13},
{7,10,16,17,18},
{8,12,13,14,18},
{9,11,13,15,17},
{9,12,13,16,17,18},
{17},
{13},
{13,17},
{18},
{13},
{13,18},
{},
{18},
{17},
{17,18},
{},
{}
}; 
 
std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};
