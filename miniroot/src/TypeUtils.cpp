#include "TypeUtils.hpp"
#include <cstdio>

void print_mem(uint8_t* mem, unsigned rows)
{
    unsigned short i = 0;
    do {
        printf("%02x%02x %02x%02x %02x%02x %02x%02x\n", 
            mem[8 * i + 0], mem[8 * i + 1], 
            mem[8 * i + 2], mem[8 * i + 3], 
            mem[8 * i + 4], mem[8 * i + 5], 
            mem[8 * i + 6], mem[8 * i + 7]);
    } while (i++ < rows);
};
