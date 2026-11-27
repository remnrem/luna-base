#pragma once
#include <string>
#include <utility>
#include <vector>

// Standard factor declarations, in writer/catalog declaration order.
// 'i': integer, 'r': real, 't': text. T is an interval in typed stores.
std::vector<std::pair<std::string,char>> luna_standard_factors();
