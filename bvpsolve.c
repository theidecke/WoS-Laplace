#include <inttypes.h>
#include <stdio.h>

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

union uint64_double_union {
    uint64_t i;
    double d;
};

double uint64_to_double(uint64_t x) {
    union uint64_double_union res;
    res.i = 0x3ff0000000000000ULL + (x >> 12); // take first 52 bits of x as mantisse of our double value
    return res.d - 1.0;
}

double uniform(pcg32_random_t* rng) {
    pcg32_random_r(rng);
    return uint64_to_double(rng->state);
}

int main(int argc, char const *argv[])
{
    pcg32_random_t rngstate = {42ULL, 0ULL};
    /* code */
    printf("State: %" PRIu64 "\n", rngstate.state);
    pcg32_random_r(&rngstate);
    printf("State: %" PRIu64 "\n", rngstate.state);
    pcg32_random_r(&rngstate);
    printf("State: %" PRIu64 "\n", rngstate.state);
    pcg32_random_r(&rngstate);
    printf("State: %" PRIu64 "\n", rngstate.state);

    printf("One = %f\n", uint64_to_double(0ULL));
    printf("One = %f\n", uint64_to_double(18446744073709551615ULL));

    printf("uniform random double: %f\n", uniform(&rngstate));
    printf("uniform random double: %f\n", uniform(&rngstate));
    printf("uniform random double: %f\n", uniform(&rngstate));
    printf("uniform random double: %f\n", uniform(&rngstate));
    printf("uniform random double: %f\n", uniform(&rngstate));

    return 0;
}