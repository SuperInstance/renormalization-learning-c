/*
 * learning_curve.c — Predict learning curves from RG critical exponents
 *
 * Learning follows a power law near criticality:
 *   p(t) = p_max - A * t^{-alpha}
 *
 * where alpha is related to the correlation-length critical exponent nu.
 * In the RG framework, the learning rate exponent alpha = 1/nu.
 * This follows from: correlation length xi ~ (g - g_c)^{-nu}
 * and the number of coarse-graining steps ~ log(t), giving the relation.
 */

#include "renorm_learn.h"
#include <math.h>
#include <stdlib.h>

int fit_power_law(const double *t, const double *p, int n, double p_max,
                  double *out_a, double *out_alpha, double *out_rmse)
{
    if (!t || !p || n < 3 || !out_a || !out_alpha) return -1;
    if (p_max <= 0.0) return -1;

    /* Fit p_max - p(t) = A * t^{-alpha}
     * ln(p_max - p) = ln(A) - alpha * ln(t)
     *
     * Line fitting: Y = ln(p_max - p), X = ln(t)
     * Y = ln(A) - alpha * X
     */
    double sum_X = 0.0, sum_Y = 0.0, sum_XX = 0.0, sum_XY = 0.0;
    int used = 0;

    for (int i = 0; i < n; i++) {
        double diff = p_max - p[i];
        if (diff <= 1e-15) continue; /* skip saturated points */
        double X = log(t[i]);
        double Y = log(diff);
        sum_X += X;
        sum_Y += Y;
        sum_XX += X * X;
        sum_XY += X * Y;
        used++;
    }

    if (used < 2) return -1;

    double denom = (double)used * sum_XX - sum_X * sum_X;
    if (fabs(denom) < 1e-15) return -1;

    double slope = ((double)used * sum_XY - sum_X * sum_Y) / denom;
    double intercept = (sum_Y - slope * sum_X) / (double)used;

    *out_alpha = -slope;
    if (*out_alpha < 0) *out_alpha = 1e-6; /* clamp to positive */
    *out_a = exp(intercept);

    /* Compute RMSE */
    if (out_rmse) {
        double ss = 0.0;
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            double diff = p_max - p[i];
            if (diff <= 1e-15) continue;
            double predicted = p_max - *out_a * pow(t[i], -*out_alpha);
            ss += (p[i] - predicted) * (p[i] - predicted);
            cnt++;
        }
        *out_rmse = (cnt > 0) ? sqrt(ss / cnt) : 0.0;
    }

    return 0;
}

int predict_mastery_time(double p_max, double a, double alpha,
                         double target, double *out_t)
{
    if (!out_t) return -1;
    if (target >= p_max) return -1;
    if (a <= 0.0 || alpha <= 0.0) return -1;

    /* p(t) = p_max - A * t^{-alpha}
     * => t = (A / (p_max - target))^{1/alpha}
     */
    double diff = p_max - target;
    if (diff <= 0.0) return -1;

    *out_t = pow(a / diff, 1.0 / alpha);
    return 0;
}

double compare_learning_curves(double alpha1, double alpha2, double tol)
{
    double diff = fabs(alpha1 - alpha2);
    if (tol <= 0.0) return diff;
    return diff / tol;
}
