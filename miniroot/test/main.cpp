#include "typeutils.hpp"
#include "miniroot.hpp"
#include <iostream>
#include <memory>

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Usage\n\tMiniroot_Test <path/to/root>\n";
        return 0;
    }
    auto        miniroot = std::make_unique<Miniroot>(argv[1]);
    const auto& branches = miniroot->GetBranches();

    std::cout << " --- Branches --- \n";
    for (const auto& b : branches) { std::cout << b << '\n'; }

    auto vec = miniroot->Get<float>(branches.at(0));

    print_mem((uint8_t*)vec.data(), 10);
    std::cout << "length = " << vec.size() << std::endl;

    return 0;
}
