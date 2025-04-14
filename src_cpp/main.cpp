
#include <iostream>

#include "renyi_entropy.h"
#include "safe_btree_set.h"
#include <absl/container/internal/btree_container.h>

int main(int argc, char *argv[])
{
    std::vector<std::vector<float>> dataset{{1, 1, 1}, {2, 1, 1}, {1, 1, 2},
        {2, 1, 2}, {3, 1, 1}, {1, 1, 3},
        {3, 3, 1}, {3, 3, 3}};
    renyi_entropy::renyi_entropy<float> calculator;
    auto tree = calculator.create_KDtree(dataset);

    auto indices = tree.nearest_points({0, 0, 0}, 3);

    // absl::container_internal::btree_container<absl::container_internal::set_params<int,
    // > btree_inst;
    btree::safe_btree_set<int> btree_inst;
    for (auto i = 0; i < 10; ++i) {
        btree_inst.insert(i);
    }

    btree::safe_btree_set<std::array<int, 3>> btree_inst2;
    btree_inst2.insert({0, 1, 2});
    btree_inst2.insert({0, 1, 3});
    auto iter = btree_inst2.find({0, 1, 2});
    if (iter != btree_inst2.end()) {
        std::cout << "Found" << std::endl;
    }

    return 0;
}
