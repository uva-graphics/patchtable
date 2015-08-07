
#include <string>
#include "filter.h"
using std::string;

int name_to_index(const string &a, int n) {
    if (a == "Input") { return INPUT_SOURCE_IMAGE; }
    else if (a == "OUT") { return n-1; }
    string sub = "Filter";
    if (a.substr(0, sub.size()) != sub) {
        ASSERT(0, "name_to_index expected Input|OUT|Filter");
    }
    return std::stoi(a.substr(sub.size()));
}

