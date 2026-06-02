/*
 * lattice.c — 2D lattice of doubles
 */

#include "renorm_learn.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

lattice_t *lattice_create(int L)
{
    assert(L > 0);

    lattice_t *lat = (lattice_t *)calloc(1, sizeof(lattice_t));
    if (!lat) return NULL;

    lat->L = L;
    lat->data = (double *)calloc((size_t)(L * L), sizeof(double));
    if (!lat->data) {
        free(lat);
        return NULL;
    }
    return lat;
}

void lattice_random(lattice_t *lat, double lo, double hi, unsigned seed)
{
    if (!lat || !lat->data) return;
    srand(seed);
    int N = lat->L * lat->L;
    double range = hi - lo;
    for (int i = 0; i < N; i++) {
        lat->data[i] = lo + range * ((double)rand() / (double)RAND_MAX);
    }
}

lattice_t *lattice_copy(const lattice_t *src)
{
    if (!src) return NULL;
    lattice_t *dst = lattice_create(src->L);
    if (!dst) return NULL;
    memcpy(dst->data, src->data, (size_t)(src->L * src->L) * sizeof(double));
    return dst;
}

void lattice_free(lattice_t *lat)
{
    if (!lat) return;
    free(lat->data);
    free(lat);
}

double lattice_magnetization(const lattice_t *lat)
{
    if (!lat || !lat->data) return 0.0;
    int N = lat->L * lat->L;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += lat->data[i];
    }
    return sum / (double)N;
}

double lattice_energy(const lattice_t *lat)
{
    if (!lat || !lat->data || lat->L < 2) return 0.0;
    int L = lat->L;
    double e = 0.0;

    for (int y = 0; y < L; y++) {
        for (int x = 0; x < L; x++) {
            double s = lat->data[y * L + x];
            /* right neighbor (periodic) */
            int xr = (x + 1) % L;
            double sr = lat->data[y * L + xr];
            e -= s * sr;
            /* bottom neighbor (periodic) */
            int yb = (y + 1) % L;
            double sb = lat->data[yb * L + x];
            e -= s * sb;
        }
    }
    return e;
}

double lattice_correlation(const lattice_t *lat, int r)
{
    if (!lat || !lat->data || r < 0 || r >= lat->L) return 0.0;
    int L = lat->L;
    double mag = lattice_magnetization(lat);
    double sum = 0.0;
    int count = 0;

    for (int y = 0; y < L; y++) {
        for (int x = 0; x < L; x++) {
            double s_i = lat->data[y * L + x];
            int xr = (x + r) % L;
            double s_j = lat->data[y * L + xr];
            sum += s_i * s_j;
            count++;
            /* also do the vertical direction to improve statistics */
            int yr = (y + r) % L;
            s_j = lat->data[yr * L + x];
            sum += s_i * s_j;
            count++;
        }
    }
    return (sum / (double)count) - (mag * mag);
}
