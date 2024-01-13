#include <cinttypes>

union float_union
{
    float f_val;
    int8_t bytes[4];
};

union double_union
{
    double f_val;
    int8_t bytes[8];
};

union char_union
{
    char f_val;
    int8_t bytes;
};

union int_union
{
    int f_val;
    int8_t bytes[4];
};

union long_union
{
    long f_val;
    int8_t bytes[8];
};

void print_mem(uint8_t* mem, unsigned short rows);
