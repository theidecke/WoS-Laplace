#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define NR_OF_CIRCLES 4 // nr of circles that make up our boundary conditions
#define MAX_WALK_LENGTH 128 // maximum nr of steps we try without hitting a boundary before stopping
#define EPSILON 0.01 // distance we consider close enough to count as boundary hit
#define RES 64 // simulation and image export resolution

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

// IMAGE EXPORT ///////////////////////////////////////////////////////

void writeImageFile(char* fn, double* imagedata, uint32_t w, uint32_t h) {
    FILE* file = fopen(fn, "w+");
    double maxvalue = 255.0;
    fprintf(file, "P2\n%u %u\n%u\n", w, h, (uint32_t)maxvalue);
    for (int row = 0; row < h; ++row) {
        for (int column = 0; column < w; ++column) {
            fprintf(file, "%u ", (uint32_t)(maxvalue * imagedata[w * row + column]));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

// DEFINE ENVIRONMENT / BOUNDARY CONDITIONS ///////////////////////////

struct point {
    double x, y;
};

struct circle {
    double center_x, center_y, radius, boundary_value_inside, boundary_value_outside;
};

void closestCirclePoint(struct point* retval, double* retbval, struct circle* c, double px, double py) {
    double dx = px - c->center_x;
    double dy = py - c->center_y;
    double r_sq = dx * dx + dy * dy;
    if (r_sq==0.0) // in case we are right in the center of the circle pick a direction
    {
        dx = 1.0;
        dy = 0.0;
        r_sq = 1.0;
    }
    double scale_to_radius = c->radius / sqrt(r_sq);
    retval->x = c->center_x + scale_to_radius * dx;
    retval->y = c->center_y + scale_to_radius * dy;
    if (r_sq < c->radius * c->radius)
    { // inside circle
        *retbval = c->boundary_value_inside;
    } else { // outside circle
        *retbval = c->boundary_value_outside;
    }

}

void closestBoundary(struct circle* circles, double* distance, double* boundary_value, double px, double py) {
    struct point closest_point;
    double closest_distance = DBL_MAX;
    double ret_boundary_value = 0.0;

    for (int i = 0; i < NR_OF_CIRCLES; ++i)
    {
        double current_distance;
        double current_boundary_value;
        closestCirclePoint(&closest_point, &current_boundary_value, &circles[i], px, py);
        current_distance = sqrt( (closest_point.x-px)*(closest_point.x-px) + (closest_point.y-py)*(closest_point.y-py) );
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

double walkOnSpheres(pcg32_random_t* rng, struct circle* env, double sx, double sy) {
    double boundary_distance, boundary_value;
    double x = sx, y = sy;
    int n = 0;
    closestBoundary(env, &boundary_distance, &boundary_value, x, y);
    while(boundary_distance > EPSILON && n < MAX_WALK_LENGTH) {
        double phi = 2.0 * M_PI * uniform(rng); // choose a random angle
        double nx, ny;
        nx = x + boundary_distance * cos(phi);
        ny = y + boundary_distance * sin(phi);
        x = nx; y = ny;
        ++n;
        closestBoundary(env, &boundary_distance, &boundary_value, x, y);
    }
    return boundary_value;
}

double averageWalkOnSpheres(pcg32_random_t* rng, struct circle* env, int repetitions, double sx, double sy) {
    double wosSum = 0.0;
    for(int i = 0; i<repetitions; ++i) {
        wosSum += walkOnSpheres(rng, env, sx, sy);
    }
    return wosSum / repetitions;
}

int main(int argc, char const *argv[])
{
    pcg32_random_t rngstate = {42ULL, 0ULL};

    struct circle environment[NR_OF_CIRCLES] = {
        {-1.0, -1.0, 1.0, 0.0, 0.0}, // x, y, radius, inside-bv, outside-bv
        { 1.0, -1.0, 1.0, 1.0, 1.0},
        { 1.0,  1.0, 1.0, 0.0, 0.0},
        {-1.0,  1.0, 1.0, 1.0, 1.0}
    };

    double imagedata[RES*RES];
    for (int j = 0; j < RES; ++j)
    {
        for (int i = 0; i < RES; ++i)
        {
            double sx = -4.0 + 8.0*((double)i+0.5)/(double)RES;
            double sy =  4.0 - 8.0*((double)j+0.5)/(double)RES;
            *(imagedata+RES*j+i) = averageWalkOnSpheres(&rngstate, &environment[0], 256, sx, sy);
        }
    }

    writeImageFile("test.pgm", imagedata, RES, RES);

    return 0;
}