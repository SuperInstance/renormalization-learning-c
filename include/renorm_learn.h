/*
 * renorm_learn.h — Renormalization-Learning C Library
 *
 * Wilsonian renormalization group as a model for agent skill acquisition.
 *   - coarse-graining = learning (integrating out microscopic fluctuations)
 *   - fixed points = perfected skills
 *   - critical exponents = learning rate predictions
 */

#ifndef RENORM_LEARN_H
#define RENORM_LEARN_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Lattice system ───────────────────────────────────────────────── */

/** A 2D lattice of double-valued "spins" or "skill levels" */
typedef struct {
    int L;           /**< linear dimension (L x L lattice) */
    double *data;    /**< row-major array of size L*L */
} lattice_t;

/**
 * Create an L x L lattice initialized to zero.
 * Returns NULL on allocation failure.
 */
lattice_t *lattice_create(int L);

/**
 * Fill lattice with uniform random values in [lo, hi).
 * Seed is deterministic for reproducibility.
 */
void lattice_random(lattice_t *lat, double lo, double hi, unsigned seed);

/**
 * Deep-copy a lattice.
 * Returns NULL on allocation failure.
 */
lattice_t *lattice_copy(const lattice_t *src);

/** Free all memory associated with a lattice. */
void lattice_free(lattice_t *lat);

/**
 * Magnetization (order parameter): average value over all sites.
 * Sum(lat->data[i]) / (L*L)
 */
double lattice_magnetization(const lattice_t *lat);

/**
 * Energy: sum over all pairs of nearest-neighbor product.
 * H = -sum_{<ij>} s_i * s_j   (Ising-like interaction)
 */
double lattice_energy(const lattice_t *lat);

/**
 * Spatial correlation function C(r) for 0 <= r < L.
 * C(r) = <s_i * s_{i+r}> - <s_i>^2  averaged over all i.
 */
double lattice_correlation(const lattice_t *lat, int r);


/* ─── Block spin transformation (Wilsonian RG) ────────────────────── */

/**
 * Coarse-grain a lattice by a block factor b using block averaging:
 * group b x b blocks, average each block → one site on reduced lattice.
 * Returns new lattice of size (L/b) x (L/b).
 * L must be divisible by b.
 */
lattice_t *block_spin_coarse_grain(const lattice_t *lat, int b);

/**
 * Rescale coupling / amplitudes back to original lattice scale.
 * Under RG, spin fields rescale as s' = s * zeta where zeta = b^{(2-d)/2}
 * for d=2: zeta = 1 (trivial rescaling for the Gaussian fixed point).
 * More generally: apply a multiplicative factor to all sites.
 */
void block_spin_rescale(lattice_t *lat, double zeta);

/** Convenience: one full RG step = coarse-grain + rescale */
lattice_t *rg_step(const lattice_t *lat, int b, double zeta);

/**
 * A multi-step RG flow stored as an array of lattices at decreasing
 * resolution.
 */
typedef struct {
    lattice_t **levels;   /**< level[0] = original lattice */
    int nlevels;          /**< number of levels stored */
    int b;                /**< block factor used */
    double zeta;          /**< rescale factor used */
} rl_rgflow_t;

/** Run RG flow for nsteps steps starting from lat. */
rl_rgflow_t *rgflow_create(const lattice_t *lat, int b, double zeta,
                            int nsteps);

/** Free all levels and the flow. */
void rgflow_free(rl_rgflow_t *flow);

/** Return lattice at level i (0 = original). */
const lattice_t *rgflow_level(const rl_rgflow_t *flow, int i);


/* ─── Beta function / fixed point ────────────────────────────────── */

/**
 * Compute the beta function for a given coupling parameter.
 * beta(g) = dg/dt  where t = ln(b), approximated by finite difference
 * over one RG step.
 *
 * @param g          current coupling value
 * @param lat        reference lattice prepared at coupling g
 * @param b          block factor
 * @param zeta       scale factor
 * @param eps        (unused, reserved for future use)
 * @param out_beta   output: beta(g)
 * @return 0 on success, -1 on error
 */
int compute_beta_function(double g, const lattice_t *lat, int b,
                          double zeta, double eps,
                          double *out_beta);

/**
 * Find a fixed point by iterating RG until the coupling stabilizes.
 * Starting from lat0, run RG steps until |magnetization(lat) - magnetization(prev)| < tol,
 * or until max_steps is reached.
 *
 * @param lat0       initial lattice
 * @param b          block factor
 * @param zeta       scale factor
 * @param tol        convergence tolerance for magnetization
 * @param max_steps  maximum number of RG steps
 * @param out_steps  output: number of steps taken
 * @param out_mag    output: fixed-point magnetization
 * @return 0 if converged, 1 if max_steps reached, -1 on error
 */
int find_fixed_point(const lattice_t *lat0, int b, double zeta,
                     double tol, int max_steps,
                     int *out_steps, double *out_mag);


/* ─── Skill acquisition as RG ─────────────────────────────────────── */

/**
 * A "SkillLattice" is a lattice where each site holds a skill level
 * in [0,1] (0 = unlearned, 1 = fully mastered).
 */
typedef lattice_t skill_lattice_t;

/**
 * Convert a skill lattice to the standard spin lattice (alias).
 * All skill_lattice functions operate on lattice_t internally.
 */

/** Create a skill lattice (wrapper). */
static inline skill_lattice_t *skill_lattice_create(int L)
{
    return lattice_create(L);
}

/** Free a skill lattice (wrapper). */
static inline void skill_lattice_free(skill_lattice_t *s)
{
    lattice_free(s);
}

/** Copy a skill lattice. */
static inline skill_lattice_t *skill_lattice_copy(const skill_lattice_t *s)
{
    return lattice_copy(s);
}

/**
 * Initialize a skill lattice with skill levels.
 * Random initialization in [0, 1] with given seed.
 */
static inline void skill_lattice_random(skill_lattice_t *s, unsigned seed)
{
    lattice_random(s, 0.0, 1.0, seed);
}

/**
 * Coarse-graining a skill lattice = learning.
 * Integrate out microscopic fluctuations by block averaging.
 * Returns a new, coarser skill lattice.
 */
static inline skill_lattice_t *skill_coarse_grain(const skill_lattice_t *s, int b)
{
    return block_spin_coarse_grain(s, b);
}

/**
 * After learning (coarse-graining), a mastered skill is an RG fixed point:
 * coarse-graining does not change it.  Return 1 if the lattice is a
 * fixed point (uniform and >= 0.999 average), 0 otherwise.
 */
int skill_is_fixed_point(const skill_lattice_t *s, double threshold);

/**
 * Critical exponent extraction.
 * Given a sequence of correlation lengths at successive RG steps,
 * fit length[i] ~ b^{i * nu} to extract nu (the correlation-length exponent).
 *
 * @param lengths   array of correlation lengths at RG steps 0..n-1
 * @param n         number of data points
 * @param b         block factor
 * @param out_nu    output: fitted critical exponent nu
 * @return 0 on success, -1 on error
 */
int skill_critical_exponent(const double *lengths, int n, int b,
                            double *out_nu);

/**
 * Determine universality class.
 * Two skills are in the same universality class if their critical
 * exponents agree within `tol`.
 *
 * @return 1 if same class, 0 otherwise
 */
int skill_same_universality_class(double nu1, double nu2, double tol);

/**
 * Predict the phase-transition behavior near the novice→expert boundary.
 * At criticality, the learning curve follows a power law.
 * Return the magnetization (order parameter).
 */
double skill_phase_transition_magnetization(const skill_lattice_t *s);


/* ─── Learning curve analysis ─────────────────────────────────────── */

/**
 * Fit performance data to a power law: p(t) = p_max - A * t^{-alpha}
 * using log-log least squares.
 *
 * @param t       array of time / practice steps
 * @param p       array of performance values (same length)
 * @param n       number of data points
 * @param p_max   known asymptotic performance limit
 * @param out_a   output: amplitude A (fitted)
 * @param out_alpha  output: exponent alpha (fitted)
 * @param out_rmse   output: root-mean-square error of fit
 * @return 0 on success, -1 on error
 */
int fit_power_law(const double *t, const double *p, int n, double p_max,
                  double *out_a, double *out_alpha, double *out_rmse);

/**
 * Predict the time (number of practice steps) to reach a target performance
 * level, using the power-law fit p(t) = p_max - A * t^{-alpha}.
 *
 * @param p_max        asymptotic maximum performance
 * @param a            amplitude (from fit)
 * @param alpha        exponent (from fit)
 * @param target        target performance < p_max
 * @param out_t        output: predicted time to reach target
 * @return 0 on success, -1 if target >= p_max or other error
 */
int predict_mastery_time(double p_max, double a, double alpha,
                         double target, double *out_t);

/**
 * Compare two learning curves by their power-law parameters.
 * Returns a similarity metric: difference in alpha normalized by tolerance.
 * Value < 1.0 means they're in the same universality class (similar learning).
 */
double compare_learning_curves(double alpha1, double alpha2, double tol);

#ifdef __cplusplus
}
#endif

#endif /* RENORM_LEARN_H */
