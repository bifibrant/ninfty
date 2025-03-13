#include <vector>

std::vector<std::string> subgroup_dictionary{
"1",
"C2",
"C4",
"C4",
"C4",
"Q8"
};

std::vector<std::pair<unsigned, unsigned>> lattice{
{0,1},
{0,3},
{0,2},
{0,4},
{0,5},
{1,3},
{1,2},
{1,4},
{1,5},
{3,5},
{2,5},
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
{11}};
 

std::vector<std::vector<unsigned>> intersections{
{},
{0},
{0},
{0},
{0,1,2,3},
{},
{},
{},
{5,6,7},
{6,7},
{5,7},
{5,6}
}; 
 

std::vector<std::vector<unsigned>> cointersections{
{},
{5,10,11},
{6,9,11},
{7,9,10},
{8,9,10,11},
{10,11},
{9,11},
{9,10},
{9,10,11},
{},
{},
{}
}; 
 
std::vector<std::string> pretty_subgroup_dictionary{
    "1",
    "C_2",
    "C_4",
    "C_4",
    "C_4",
    "Q_8"
};

std::vector<std::string> vertex_layout{
    "(2,0)",
    "(2,0.803)",
    "(2,1.76)",
    "(3.88,1.76)",
    "(0.125,1.76)",
    "(2,2.71)"
};

std::vector<std::vector<std::string>> edge_options{
    {"","","[bend right]","","","[bend left]"},
    {"","","","","","[bend left]"},
    {"","","","","",""},
    {"","","","","",""},
    {"","","","","",""},
    {"","","","","",""},
};
