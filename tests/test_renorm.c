/*
 * test_renorm.c — Comprehensive test suite
 *
 * Tests cover: lattice operations, block spin RG, fixed points,
 * skill model, learning curve analysis, and edge cases.
 * Minimum 25 tests.
 */

#include "renorm_learn.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    printf("  TEST %-55s ", #name); \
    if (test_##name()) { \
        printf("PASS\n"); \
        tests_passed++; \
    } else { \
        printf("FAIL\n"); \
        tests_failed++; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, eps) do { \
    if (fabs((a) - (b)) > (eps)) { \
        fprintf(stderr, "  ASSERT_NEAR failed at %s:%d: " \
                "%.10f vs %.10f (eps %.10f)\n", \
                __FILE__, __LINE__, (double)(a), (double)(b), (double)(eps)); \
        return 0; \
    } \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "  ASSERT_TRUE failed at %s:%d: %s\n", \
                __FILE__, __LINE__, #cond); \
        return 0; \
    } \
} while(0)

/* ================================================================
 * Test 1: Lattice creation
 * ================================================================ */
static int test_lattice_create(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    ASSERT_TRUE(lat->L == 8);
    ASSERT_TRUE(lat->data != NULL);
    /* Should be zero-initialized */
    for (int i = 0; i < 64; i++) {
        ASSERT_NEAR(lat->data[i], 0.0, 1e-15);
    }
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 2: Lattice magnetization of uniform lattice
 * ================================================================ */
static int test_magnetization_uniform(void)
{
    lattice_t *lat = lattice_create(4);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 16; i++) lat->data[i] = 1.0;
    ASSERT_NEAR(lattice_magnetization(lat), 1.0, 1e-15);
    for (int i = 0; i < 16; i++) lat->data[i] = -0.5;
    ASSERT_NEAR(lattice_magnetization(lat), -0.5, 1e-15);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 3: Lattice energy of ferromagnetic state
 * ================================================================ */
static int test_energy_ferromagnetic(void)
{
    int L = 8;
    lattice_t *lat = lattice_create(L);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < L * L; i++) lat->data[i] = 1.0;
    /* For 2D square lattice with periodic BC:
     * Each site has 4 neighbors, so total bonds = L*L*2 (because each
     * bond counted once = L*L*4/2). Energy = -sum(s_i * s_j) over bonds.
     * For all spins = 1: energy = - (L*L*2) * 1 = -128 for L=8
     */
    double e = lattice_energy(lat);
    ASSERT_NEAR(e, -128.0, 1e-10);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 4: Lattice energy of anti-ferromagnetic alternating
 * ================================================================ */
static int test_energy_antiferro(void)
{
    int L = 4;
    lattice_t *lat = lattice_create(L);
    ASSERT_TRUE(lat != NULL);
    /* Checkerboard pattern: +1, -1 alternating */
    for (int y = 0; y < L; y++) {
        for (int x = 0; x < L; x++) {
            lat->data[y * L + x] = ((x + y) % 2 == 0) ? 1.0 : -1.0;
        }
    }
    double e = lattice_energy(lat);
    /* In 2D square lattice with NN+1, -1 alternating:
     * Each +1 has four -1 neighbors => product = -1, so contribution +1 for each bond.
     * Number of bonds = L*L*2 = 32. Expected energy = 32.
     */
    ASSERT_NEAR(e, 32.0, 1e-10);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 5: Lattice correlation C(0) = variance
 * ================================================================ */
static int test_correlation_zero(void)
{
    lattice_t *lat = lattice_create(4);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 16; i++) lat->data[i] = (i % 2 == 0) ? 1.0 : -1.0;
    double c0 = lattice_correlation(lat, 0);
    /* C(0) = <s_i^2> - <s_i>^2 = 1 - 0 = 1 */
    ASSERT_NEAR(c0, 1.0, 1e-10);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 6: Lattice correlation on uniform lattice = 0
 * ================================================================ */
static int test_correlation_uniform(void)
{
    lattice_t *lat = lattice_create(4);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 16; i++) lat->data[i] = 1.0;
    double c = lattice_correlation(lat, 1);
    /* <s_i * s_j> = 1, <s_i>^2 = 1, so C(r) = 0 */
    ASSERT_NEAR(c, 0.0, 1e-10);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 7: Lattice copy
 * ================================================================ */
static int test_lattice_copy(void)
{
    lattice_t *lat = lattice_create(6);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 42);
    lattice_t *cpy = lattice_copy(lat);
    ASSERT_TRUE(cpy != NULL);
    ASSERT_TRUE(cpy->L == lat->L);
    int N = lat->L * lat->L;
    for (int i = 0; i < N; i++) {
        ASSERT_NEAR(cpy->data[i], lat->data[i], 1e-15);
    }
    lattice_free(lat);
    lattice_free(cpy);
    return 1;
}

/* ================================================================
 * Test 8: Lattice free with NULL (shouldn't crash)
 * ================================================================ */
static int test_lattice_free_null(void)
{
    lattice_free(NULL);
    return 1;
}

/* ================================================================
 * Test 9: Block spin reduces lattice by factor 2
 * ================================================================ */
static int test_block_spin_reduces(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 123);
    lattice_t *cg = block_spin_coarse_grain(lat, 2);
    ASSERT_TRUE(cg != NULL);
    ASSERT_TRUE(cg->L == 4);  /* 8/2 = 4 */
    lattice_free(lat);
    lattice_free(cg);
    return 1;
}

/* ================================================================
 * Test 10: Block spin reduces by factor 4
 * ================================================================ */
static int test_block_spin_reduces4(void)
{
    lattice_t *lat = lattice_create(16);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 456);
    lattice_t *cg = block_spin_coarse_grain(lat, 4);
    ASSERT_TRUE(cg != NULL);
    ASSERT_TRUE(cg->L == 4);  /* 16/4 = 4 */
    lattice_free(lat);
    lattice_free(cg);
    return 1;
}

/* ================================================================
 * Test 11: Magnetization preserved for uniform lattice under coarse-graining
 * ================================================================ */
static int test_mag_preserved_uniform(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 64; i++) lat->data[i] = 0.75;
    double mag0 = lattice_magnetization(lat);
    lattice_t *cg = block_spin_coarse_grain(lat, 2);
    ASSERT_TRUE(cg != NULL);
    double mag_cg = lattice_magnetization(cg);
    ASSERT_NEAR(mag_cg, mag0, 1e-15);
    lattice_free(lat);
    lattice_free(cg);
    return 1;
}

/* ================================================================
 * Test 12: Fixed point detection on uniform lattice -> 1-step convergence
 * ================================================================ */
static int test_fixed_point_uniform(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 64; i++) lat->data[i] = 1.0;

    int steps;
    double mag;
    int result = find_fixed_point(lat, 2, 1.0, 1e-10, 20, &steps, &mag);
    ASSERT_TRUE(result == 0);  /* converged */
    ASSERT_NEAR(mag, 1.0, 1e-10);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 13: Fixed point detection on nearly-uniform lattice
 * ================================================================ */
static int test_fixed_point_nearly_uniform(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 64; i++) lat->data[i] = 0.999;
    /* Small perturbation */
    lat->data[0] += 0.001;

    int steps;
    double mag;
    int result = find_fixed_point(lat, 2, 1.0, 1e-6, 20, &steps, &mag);
    ASSERT_TRUE(result == 0);  /* should converge */
    ASSERT_NEAR(mag, 1.0, 0.01);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 14: Beta function computation
 * ================================================================ */
static int test_beta_function(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 99);

    double beta;
    int ret = compute_beta_function(1.0, lat, 2, 1.0, 0.01, &beta);
    ASSERT_TRUE(ret == 0);
    /* Beta should return a finite number (no check on sign for random) */
    ASSERT_TRUE(isfinite(beta));
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 15: Beta function on uniform lattice should be zero
 * ================================================================ */
static int test_beta_function_uniform(void)
{
    lattice_t *lat = lattice_create(8);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 64; i++) lat->data[i] = 1.0;

    double beta;
    int ret = compute_beta_function(1.0, lat, 2, 1.0, 0.01, &beta);
    ASSERT_TRUE(ret == 0);
    /* For uniform lattice, beta should be near zero since uniform is a fixed point */
    ASSERT_NEAR(beta, 0.0, 0.01);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 16: RG flow creation
 * ================================================================ */
static int test_rgflow_create(void)
{
    lattice_t *lat = lattice_create(16);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 42);

    rl_rgflow_t *flow = rgflow_create(lat, 2, 1.0, 3);
    ASSERT_TRUE(flow != NULL);
    ASSERT_TRUE(flow->nlevels == 4);  /* original + 3 steps */
    ASSERT_TRUE(flow->levels[0]->L == 16);
    ASSERT_TRUE(flow->levels[1]->L == 8);
    ASSERT_TRUE(flow->levels[2]->L == 4);
    ASSERT_TRUE(flow->levels[3]->L == 2);

    rgflow_free(flow);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 17: Skill lattice is a fixed point when mastered
 * ================================================================ */
static int test_skill_fixed_point_mastered(void)
{
    skill_lattice_t *s = skill_lattice_create(8);
    ASSERT_TRUE(s != NULL);
    for (int i = 0; i < 64; i++) s->data[i] = 1.0;

    ASSERT_TRUE(skill_is_fixed_point(s, 0.99));

    skill_lattice_free(s);
    return 1;
}

/* ================================================================
 * Test 18: Skill lattice is NOT a fixed point when not uniform
 * ================================================================ */
static int test_skill_not_fixed(void)
{
    skill_lattice_t *s = skill_lattice_create(8);
    ASSERT_TRUE(s != NULL);
    for (int i = 0; i < 64; i++) s->data[i] = (i % 2 == 0) ? 0.8 : 0.6;

    ASSERT_TRUE(!skill_is_fixed_point(s, 0.99));

    skill_lattice_free(s);
    return 1;
}

/* ================================================================
 * Test 19: Skill coarse-graining preserves mastered skills
 * ================================================================ */
static int test_skill_coarse_grain_preserves(void)
{
    skill_lattice_t *s = skill_lattice_create(16);
    ASSERT_TRUE(s != NULL);
    for (int i = 0; i < 256; i++) s->data[i] = 1.0;

    skill_lattice_t *cg = skill_coarse_grain(s, 2);
    ASSERT_TRUE(cg != NULL);
    double mag_cg = lattice_magnetization(cg);
    ASSERT_NEAR(mag_cg, 1.0, 1e-15);

    /* Should still be a fixed point */
    ASSERT_TRUE(skill_is_fixed_point(cg, 0.99));

    skill_lattice_free(s);
    skill_lattice_free(cg);
    return 1;
}

/* ================================================================
 * Test 20: Skill critical exponent extraction from known data
 * ================================================================ */
static int test_skill_critical_exponent(void)
{
    /* Generate synthetic data: lengths[i] = b^{i * nu} with nu = 1.0 */
    int b = 2;
    double nu_true = 1.0;
    double lengths[4];
    for (int i = 0; i < 4; i++) {
        lengths[i] = pow((double)b, (double)i * nu_true);
    }

    double nu_fit;
    int ret = skill_critical_exponent(lengths, 4, b, &nu_fit);
    ASSERT_TRUE(ret == 0);
    ASSERT_NEAR(nu_fit, nu_true, 0.01);

    return 1;
}

/* ================================================================
 * Test 21: Universality class comparison
 * ================================================================ */
static int test_universality_class(void)
{
    /* Same class */
    ASSERT_TRUE(skill_same_universality_class(1.0, 1.01, 0.05));
    /* Different class */
    ASSERT_TRUE(!skill_same_universality_class(1.0, 2.0, 0.1));
    return 1;
}

/* ================================================================
 * Test 22: Power law fitting
 * ================================================================ */
static int test_power_law_fit(void)
{
    /* Generate synthetic data: p(t) = 0.95 - 0.5 * t^{-0.4}
     * All values positive and increasing toward p_max.
     * This is a well-behaved learning curve approaching mastery.
     */
    double p_max = 0.95;
    double a_true = 0.5;
    double alpha_true = 0.4;

    int n = 15;
    double t[15], p[15];
    for (int i = 0; i < n; i++) {
        t[i] = (double)(i * 10 + 5);  /* t = 5, 15, 25, ... 145 */
        double diff = a_true * pow(t[i], -alpha_true);
        p[i] = p_max - diff;
    }

    double a_fit, alpha_fit, rmse;
    int ret = fit_power_law(t, p, n, p_max, &a_fit, &alpha_fit, &rmse);
    ASSERT_TRUE(ret == 0);
    ASSERT_NEAR(a_fit, a_true, 0.01);
    ASSERT_NEAR(alpha_fit, alpha_true, 0.01);

    return 1;
}

/* ================================================================
 * Test 23: Predict mastery time
 * ================================================================ */
static int test_predict_mastery(void)
{
    double p_max = 1.0;
    double a = 2.0;
    double alpha = 0.5;
    double target = 0.9;

    /* t = (a / (p_max - target))^{1/alpha} = (2.0 / 0.1)^2 = 400 */
    double t_pred;
    int ret = predict_mastery_time(p_max, a, alpha, target, &t_pred);
    ASSERT_TRUE(ret == 0);
    ASSERT_NEAR(t_pred, 400.0, 1.0);

    return 1;
}

/* ================================================================
 * Test 24: Compare learning curves
 * ================================================================ */
static int test_compare_curves(void)
{
    double tol = 0.1;
    /* Similar curves: diff / tol < 1 */
    double sim1 = compare_learning_curves(0.5, 0.52, tol);
    ASSERT_TRUE(sim1 < 1.0);
    /* Dissimilar curves: diff / tol > 1 */
    double sim2 = compare_learning_curves(0.5, 0.8, tol);
    ASSERT_TRUE(sim2 > 1.0);
    /* Identical */
    double sim3 = compare_learning_curves(0.5, 0.5, 1.0);
    ASSERT_NEAR(sim3, 0.0, 1e-15);
    return 1;
}

/* ================================================================
 * Test 25: Phase transition magnetization
 * ================================================================ */
static int test_phase_transition_mag(void)
{
    lattice_t *lat = lattice_create(4);
    ASSERT_TRUE(lat != NULL);
    for (int i = 0; i < 16; i++) lat->data[i] = (i < 8) ? 0.8 : 0.2;
    double mag = skill_phase_transition_magnetization(lat);
    /* Average of (8*0.8 + 8*0.2) / 16 = 0.5 */
    ASSERT_NEAR(mag, 0.5, 1e-15);
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 26: Edge case — single-site lattice
 * ================================================================ */
static int test_single_site(void)
{
    lattice_t *lat = lattice_create(1);
    ASSERT_TRUE(lat != NULL);
    lat->data[0] = 1.0;
    ASSERT_NEAR(lattice_magnetization(lat), 1.0, 1e-15);
    ASSERT_NEAR(lattice_energy(lat), 0.0, 1e-15);           /* no neighbors */
    ASSERT_NEAR(lattice_correlation(lat, 0), 0.0, 1e-15);   /* C(0) = 1-1^2 = 0 */
    lattice_free(lat);
    return 1;
}

/* ================================================================
 * Test 27: Edge case — trivial RG step (no reduction)
 * ================================================================ */
static int test_rg_b1(void)
{
    lattice_t *lat = lattice_create(4);
    ASSERT_TRUE(lat != NULL);
    lattice_random(lat, -1.0, 1.0, 77);
    lattice_t *cg = block_spin_coarse_grain(lat, 1);
    ASSERT_TRUE(cg != NULL);
    ASSERT_TRUE(cg->L == 4);  /* block=1 means no actual reduction */
    double mag0 = lattice_magnetization(lat);
    double mag1 = lattice_magnetization(cg);
    ASSERT_NEAR(mag1, mag0, 1e-15);
    lattice_free(lat);
    lattice_free(cg);
    return 1;
}

/* ================================================================
 * Test 28: Skill lattice copy preserves values
 * ================================================================ */
static int test_skill_copy(void)
{
    skill_lattice_t *s = skill_lattice_create(5);
    ASSERT_TRUE(s != NULL);
    for (int i = 0; i < 25; i++) s->data[i] = 0.314 * (i + 1);
    skill_lattice_t *cpy = skill_lattice_copy(s);
    ASSERT_TRUE(cpy != NULL);
    for (int i = 0; i < 25; i++) {
        ASSERT_NEAR(cpy->data[i], s->data[i], 1e-15);
    }
    skill_lattice_free(s);
    skill_lattice_free(cpy);
    return 1;
}

/* ================================================================
 * Test 29: Two different lattices with same critical exponent
 * ================================================================ */
static int test_two_lattices_same_exponent(void)
{
    /* Create two synthetic datasets with same nu=1.5 */
    int b = 2;
    double lengths1[4], lengths2[4];
    for (int i = 0; i < 4; i++) {
        lengths1[i] = pow(b, 1.5 * i) * (1.0 + 0.01 * (rand() % 100 - 50) / 50.0);
        lengths2[i] = pow(b, 1.5 * i) * (1.0 + 0.01 * (rand() % 100 - 50) / 50.0);
    }

    double nu1, nu2;
    int r1 = skill_critical_exponent(lengths1, 4, b, &nu1);
    int r2 = skill_critical_exponent(lengths2, 4, b, &nu2);
    ASSERT_TRUE(r1 == 0 && r2 == 0);
    ASSERT_TRUE(skill_same_universality_class(nu1, nu2, 0.2));
    return 1;
}

/* ================================================================
 * Test 30: Mastery time error on invalid input
 * ================================================================ */
static int test_mastery_time_invalid(void)
{
    double t;
    int ret = predict_mastery_time(1.0, 2.0, 0.5, 1.5, &t);
    ASSERT_TRUE(ret == -1);  /* target > p_max should fail */
    return 1;
}

int main(void)
{
    printf("Renormalization-Learning-C Test Suite\n");
    printf("=====================================\n\n");

    TEST(lattice_create);
    TEST(magnetization_uniform);
    TEST(energy_ferromagnetic);
    TEST(energy_antiferro);
    TEST(correlation_zero);
    TEST(correlation_uniform);
    TEST(lattice_copy);
    TEST(lattice_free_null);
    TEST(block_spin_reduces);
    TEST(block_spin_reduces4);
    TEST(mag_preserved_uniform);
    TEST(fixed_point_uniform);
    TEST(fixed_point_nearly_uniform);
    TEST(beta_function);
    TEST(beta_function_uniform);
    TEST(rgflow_create);
    TEST(skill_fixed_point_mastered);
    TEST(skill_not_fixed);
    TEST(skill_coarse_grain_preserves);
    TEST(skill_critical_exponent);
    TEST(universality_class);
    TEST(power_law_fit);
    TEST(predict_mastery);
    TEST(compare_curves);
    TEST(phase_transition_mag);
    TEST(single_site);
    TEST(rg_b1);
    TEST(skill_copy);
    TEST(two_lattices_same_exponent);
    TEST(mastery_time_invalid);

    printf("\n=====================================\n");
    printf("Results: %d passed, %d failed out of %d\n",
           tests_passed, tests_failed, tests_passed + tests_failed);

    return tests_failed == 0 ? 0 : 1;
}
