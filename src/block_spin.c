/*
 * block_spin.c — Wilsonian RG: coarse-graining and rescaling
 */

#include "renorm_learn.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

lattice_t *block_spin_coarse_grain(const lattice_t *lat, int b)
{
    if (!lat || !lat->data) return NULL;
    int L = lat->L;
    assert(L % b == 0);

    int newL = L / b;
    lattice_t *result = lattice_create(newL);
    if (!result) return NULL;

    for (int y = 0; y < newL; y++) {
        for (int x = 0; x < newL; x++) {
            double sum = 0.0;
            for (int dy = 0; dy < b; dy++) {
                for (int dx = 0; dx < b; dx++) {
                    int sx = x * b + dx;
                    int sy = y * b + dy;
                    sum += lat->data[sy * L + sx];
                }
            }
            result->data[y * newL + x] = sum / (double)(b * b);
        }
    }

    return result;
}

void block_spin_rescale(lattice_t *lat, double zeta)
{
    if (!lat || !lat->data) return;
    int N = lat->L * lat->L;
    for (int i = 0; i < N; i++) {
        lat->data[i] *= zeta;
    }
}

lattice_t *rg_step(const lattice_t *lat, int b, double zeta)
{
    lattice_t *coarse = block_spin_coarse_grain(lat, b);
    if (!coarse) return NULL;
    block_spin_rescale(coarse, zeta);
    return coarse;
}

rl_rgflow_t *rgflow_create(const lattice_t *lat, int b, double zeta,
                            int nsteps)
{
    if (!lat || nsteps < 0) return NULL;

    rl_rgflow_t *flow = (rl_rgflow_t *)calloc(1, sizeof(rl_rgflow_t));
    if (!flow) return NULL;

    flow->b = b;
    flow->zeta = zeta;
    flow->nlevels = nsteps + 1;  /* include original */

    flow->levels = (lattice_t **)calloc((size_t)flow->nlevels,
                                         sizeof(lattice_t *));
    if (!flow->levels) {
        free(flow);
        return NULL;
    }

    /* level 0: copy of original */
    flow->levels[0] = lattice_copy(lat);
    if (!flow->levels[0]) {
        free(flow->levels);
        free(flow);
        return NULL;
    }

    for (int i = 1; i < flow->nlevels; i++) {
        flow->levels[i] = rg_step(flow->levels[i - 1], b, zeta);
        if (!flow->levels[i]) {
            /* free what we have so far */
            for (int j = 0; j < i; j++) {
                lattice_free(flow->levels[j]);
            }
            free(flow->levels);
            free(flow);
            return NULL;
        }
    }

    return flow;
}

void rgflow_free(rl_rgflow_t *flow)
{
    if (!flow) return;
    for (int i = 0; i < flow->nlevels; i++) {
        lattice_free(flow->levels[i]);
    }
    free(flow->levels);
    free(flow);
}

const lattice_t *rgflow_level(const rl_rgflow_t *flow, int i)
{
    if (!flow || i < 0 || i >= flow->nlevels) return NULL;
    return flow->levels[i];
}

int compute_beta_function(double g, const lattice_t *lat, int b,
                          double zeta, double eps,
                          double *out_beta)
{
    (void)eps;  /* reserved for future use */
    if (!lat || !out_beta) return -1;

    /* Scale the lattice to coupling g (uniform scaling of all sites) */
    lattice_t *scaled = lattice_copy(lat);
    if (!scaled) return -1;

    int N = lat->L * lat->L;
    double cur_mag = lattice_magnetization(lat);
    if (fabs(cur_mag) < 1e-15) cur_mag = 1.0;  /* avoid division by zero */
    for (int i = 0; i < N; i++) {
        scaled->data[i] *= g / cur_mag;
    }

    /* One RG step: coarse-grain + rescale */
    lattice_t *cg = rg_step(scaled, b, zeta);
    lattice_free(scaled);

    if (!cg) return -1;

    /* Coupling after RG = magnetization of coarse-grained lattice */
    double g_prime = lattice_magnetization(cg);
    lattice_free(cg);

    /* beta(g) = (g' - g) / ln(b)
     * This is the RG beta function: dg/dt where t = ln(b).
     * For a fixed point, g' = g, so beta(g) = 0.
     */
    *out_beta = (g_prime - g) / log((double)b);

    return 0;
}

int find_fixed_point(const lattice_t *lat0, int b, double zeta,
                     double tol, int max_steps,
                     int *out_steps, double *out_mag)
{
    if (!lat0) return -1;

    lattice_t *lat = lattice_copy(lat0);
    if (!lat) return -1;

    double prev_mag = lattice_magnetization(lat);
    int steps = 0;
    int result = 1; /* default: max_steps reached */

    for (steps = 0; steps < max_steps; steps++) {
        lattice_t *next = rg_step(lat, b, zeta);
        if (!next) {
            lattice_free(lat);
            return -1;
        }
        lattice_free(lat);
        lat = next;

        double mag = lattice_magnetization(lat);
        if (fabs(mag - prev_mag) < tol) {
            result = 0; /* converged */
            prev_mag = mag;
            steps++;
            break;
        }
        prev_mag = mag;
    }

    if (out_steps) *out_steps = steps;
    if (out_mag) *out_mag = prev_mag;

    lattice_free(lat);
    return result;
}
