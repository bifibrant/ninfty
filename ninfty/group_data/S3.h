#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C2",
"C2",
"C3",
"S3"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,2},
{0,4},
{0,3},
{0,5},
{1,5},
{2,5},
{4,5},
{3,5}
}; 
 

std::vector<std::vector<unsigned>> conjugates{
{0,1,3},
{0,1,3},
{2},
{0,1,3},
{4},
{5,6,8},
{5,6,8},
{7},
{5,6,8}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{},
{},
{},
{0,1,2,3},
{1,2,3},
{0,2,3},
{0,1,3},
{0,1,2}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{6,7,8},
{5,7,8},
{5,6,8},
{5,6,7},
{5,6,7,8},
{},
{},
{},
{}
}; 
 
std::vector<std::string> pretty_subgroup_dictionary{
    "1",
    "C_2",
    "C_3",
    "S_3"
};

std::vector<std::string> vertex_layout{
    "(1,0)",
    "(1.88, 0.803)",
    "(0.125, 0.803)",
    "(1, 1.61)"
};
std::vector<std::vector<std::string>> edge_options{};
