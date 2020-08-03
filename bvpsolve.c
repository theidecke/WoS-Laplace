#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define NR_OF_CIRCLES 4 // nr of circles that make up our boundary conditions
#define MAX_WALK_LENGTH 128 // maximum nr of steps we try without hitting a boundary before stopping
#define EPSILON 0.01f // distance we consider close enough to count as boundary hit
#define RES 64 // simulation and image export resolution
#define SPP 256 // number of WoS runs (samples) per pixel

// RANDOM NUMBER GENERATION /////////////////////////////////////////// 
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

union uint32_float_union {
    uint32_t i;
    float f;
};

float float_from_uint32(uint32_t x) {
    union uint32_float_union res;
    res.i = 0x3f800000 + (x >> 9); // take first 23 bits of x as mantisse of our float value
    return res.f - 1.0f;
}

float uniform(pcg32_random_t* rng) {
    return float_from_uint32(pcg32_random_r(rng));
}

// IMAGE EXPORT ///////////////////////////////////////////////////////

void writeImageFile(char* fn, float* imagedata, uint32_t w, uint32_t h) {
    FILE* file = fopen(fn, "w+");
    double maxvalue = 255.0;
    fprintf(file, "P2\n%u %u\n%u\n", w, h, (uint32_t)maxvalue);
    for (int row = 0; row < h; ++row) {
        for (int column = 0; column < w; ++column) {
            double pixel_value = (maxvalue + 1.0) * imagedata[w * row + column]; // we scale by 256 but clamp to 255
            if(pixel_value > maxvalue) pixel_value = maxvalue;
            fprintf(file, "%u ", (uint32_t)(pixel_value));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

// DEFINE ENVIRONMENT / BOUNDARY CONDITIONS ///////////////////////////

struct point {
    float x, y;
};

struct circle {
    float center_x, center_y, radius, boundary_value_inside, boundary_value_outside;
};

void closestCirclePoint(struct point* retval, float* retbval, struct circle* c, float px, float py) {
    float dx = px - c->center_x;
    float dy = py - c->center_y;
    float r_sq = dx * dx + dy * dy;
    if (r_sq==0.0f) // in case we are right in the center of the circle pick a direction
    {
        dx = 1.0f;
        dy = 0.0f;
        r_sq = 1.0f;
    }
    float scale_to_radius = c->radius / sqrtf(r_sq);
    retval->x = c->center_x + scale_to_radius * dx;
    retval->y = c->center_y + scale_to_radius * dy;
    if (r_sq < c->radius * c->radius)
    { // inside circle
        *retbval = c->boundary_value_inside;
    } else { // outside circle
        *retbval = c->boundary_value_outside;
    }

}

void closestBoundary(struct circle* circles, float* distance, float* boundary_value, float px, float py) {
    struct point closest_point;
    float closest_distance = DBL_MAX;
    float ret_boundary_value = 0.0f;

    for (int i = 0; i < NR_OF_CIRCLES; ++i)
    {
        float current_distance;
        float current_boundary_value;
        closestCirclePoint(&closest_point, &current_boundary_value, &circles[i], px, py);
        current_distance = sqrtf( (closest_point.x-px)*(closest_point.x-px) + (closest_point.y-py)*(closest_point.y-py) );
        if (current_distance < closest_distance)
        {
            closest_distance = current_distance;
            ret_boundary_value = current_boundary_value;
        }
    }

    *distance = closest_distance;
    *boundary_value = ret_boundary_value;
}

// MAIN SIMULATION ////////////////////////////////////////////////////

float walkOnSpheres(pcg32_random_t* rng, struct circle* env, float sx, float sy) {
    float boundary_distance, boundary_value;
    float x = sx, y = sy;
    int n = 0;
    closestBoundary(env, &boundary_distance, &boundary_value, x, y);
    while(boundary_distance > EPSILON && n < MAX_WALK_LENGTH) {
        float phi = 2.0f * M_PI * uniform(rng); // choose a random angle
        float nx, ny;
        nx = x + boundary_distance * cosf(phi);
        ny = y + boundary_distance * sinf(phi);
        x = nx; y = ny;
        ++n;
        closestBoundary(env, &boundary_distance, &boundary_value, x, y);
    }
    return boundary_value;
}

float averageWalkOnSpheres(pcg32_random_t* rng, struct circle* env, int repetitions, float sx, float sy) {
    float wosSum = 0.0f;
    for(int i = 0; i<repetitions; ++i) {
        wosSum += walkOnSpheres(rng, env, sx, sy);
    }
    return wosSum / repetitions;
}

int main(int argc, char const *argv[])
{
    pcg32_random_t rngstate = {42ULL, 0ULL};

    struct circle environment[NR_OF_CIRCLES] = {
        {-1.0f, -1.0f, 1.0f, 0.0f, 0.0f}, // x, y, radius, inside-bv, outside-bv
        { 1.0f, -1.0f, 1.0f, 1.0f, 1.0f},
        { 1.0f,  1.0f, 1.0f, 0.0f, 0.0f},
        {-1.0f,  1.0f, 1.0f, 1.0f, 1.0f}
    };

    float imagedata[RES*RES];
    for (int j = 0; j < RES; ++j)
    {
        for (int i = 0; i < RES; ++i)
        {
            float sx = -4.0f + 8.0f*((float)i+0.5f)/(float)RES;
            float sy =  4.0f - 8.0f*((float)j+0.5f)/(float)RES;
            *(imagedata+RES*j+i) = averageWalkOnSpheres(&rngstate, &environment[0], SPP, sx, sy);
        }
    }

    writeImageFile("test.pgm", imagedata, RES, RES);

    return 0;
}