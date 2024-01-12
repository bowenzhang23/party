#include "miniroot.hpp"
#include <iostream>
#include <memory>

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Usage\n\tMiniroot_Test <path/to/root>\n";
        return 0;
    }
    auto miniroot = std::make_unique<Miniroot>(argv[1]);

    miniroot->Read();
    const auto& data = miniroot->GetData();

    auto print_mem = [](uint8_t* mem, unsigned short rows) {
        unsigned short i = 0;
        do {
            printf("%02x%02x%02x%02x\n", mem[4 * i + 0], mem[4 * i + 1],
                   mem[4 * i + 2], mem[4 * i + 3]);
        } while (i++ < rows);
        std::cout << std::dec;
    };

    for (const auto& p : data) {
        std::cout << p.first << ", " << p.second.size() << '\n';
        auto decomp = std::make_unique<uint8_t>();
        decomp.reset(miniroot->GetUncompressedBytes(p.second.at(0)));
        print_mem((uint8_t*) decomp.get(), 5);
    }

    return 0;
}