/*
 * skill_model.c — Skill acquisition as renormalization group
 *
 * Analogy:
 *   - Skill level s_i in [0,1] = spin value on lattice
 *   - Coarse-graining (block average) = learning process
 *   - Fixed point (uniform s ≈ 1) = mastered skill
 *   - Critical exponent nu = learning rate
 *   - Universality class = skills that learn the same way
 */

#include "renorm_learn.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int skill_is_fixed_point(const skill_lattice_t *s, double threshold)
{
    if (!s || !s->data) return 0;
    double mag = lattice_magnetization(s);
    if (mag < threshold) return 0;

    /* Check uniformity: max deviation from mean */
    int N = s->L * s->L;
    for (int i = 0; i < N; i++) {
        if (fabs(s->data[i] - mag) > 1.0 - threshold) return 0;
    }
    return 1;
}

int skill_critical_exponent(const double *lengths, int n, int b,
                            double *out_nu)
{
    if (!lengths || n < 2 || b < 2 || !out_nu) return -1;

    /* Fit length[i] ~ C * b^{i * nu}
     * ln(length[i]) = ln(C) + i * nu * ln(b)
     * So slope of ln(length[i]) vs i gives nu * ln(b)
     */
    double sum_x = 0.0, sum_y = 0.0;
    double sum_xx = 0.0, sum_xy = 0.0;

    for (int i = 0; i < n; i++) {
        double x = (double)i;
        double y = log(lengths[i]);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    double denom = (double)n * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) return -1;

    double slope = ((double)n * sum_xy - sum_x * sum_y) / denom;
    *out_nu = slope / log((double)b);

    if (*out_nu < 0) *out_nu = 0; /* exponent must be non-negative */

    return 0;
}

int skill_same_universality_class(double nu1, double nu2, double tol)
{
    return fabs(nu1 - nu2) <= tol;
}

double skill_phase_transition_magnetization(const skill_lattice_t *s)
{
    return lattice_magnetization(s);
}
