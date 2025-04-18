#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C3",
"C4",
"C6",
"C8",
"C12",
"C24"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,2},
{0,1},
{0,4},
{0,3},
{0,6},
{0,5},
{0,7},
{2,4},
{2,6},
{2,7},
{1,4},
{1,3},
{1,6},
{1,5},
{1,7},
{4,6},
{4,7},
{3,6},
{3,5},
{3,7},
{6,7},
{5,7}
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
{18},
{19},
{20},
{21}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{0,1},
{1},
{0,1,2,3},
{1,3},
{0,1,2,3,4,5},
{1},
{1,7,3},
{1,7,3,8,5},
{0},
{},
{0,10,11},
{11},
{0,10,11,12,13},
{11},
{11,15,13},
{0,10},
{},
{0,10,17,18},
{18},
{0,10,17}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{10,17,21},
{7},
{7,10,17,21},
{8,11,15},
{8,12,15,17,21},
{9,13,16,18,20},
{9,14,16,19,20,21},
{},
{15},
{16,20},
{17,21},
{15},
{15,17,21},
{16,18,20},
{16,19,20,21},
{},
{20},
{21},
{20},
{20,21},
{},
{}
}; 
 
std::vector<std::string> pretty_subgroup_dictionary{};
std::vector<std::string> vertex_layout{};
std::vector<std::vector<std::string>> edge_options{};
