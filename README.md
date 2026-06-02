# Renormalization-Learning-C

**Why do some skills take 10 hours and others take 10,000? The answer isn't complexity — it's the critical exponent.**

This is a C11 library that models skill acquisition using the Wilsonian renormalization group — the same mathematics physicists use to understand phase transitions in magnets and fluids. The payoff: you can *predict learning curves from first principles*.

## The Insight That Changes How You See Learning

**Learning IS coarse-graining.** When you practice, you're not adding detail — you're *removing* it. A mastered skill is a fixed point.

Think about learning to type. As a beginner, every keystroke is a separate, noisy event — your performance is a jagged, high-resolution field of random outcomes. As you improve, those individual fluctuations average out. Your fingers just *go*. The microscopic detail vanishes into smooth, automatic performance.

This is exactly what happens in a block-spin renormalization transformation on a lattice:

```
         Novice                              Expert
   ┌─┬─┬─┬─┬─┬─┬─┬─┐                 ┌─┬─┬─┬─┐
   │░│▓│░│░│▓│░│▓│░│                 │█│█│█│█│
   ├─┼─┼─┼─┼─┼─┼─┼─┤   block-spin    ├─┼─┼─┼─┤
   │▓│░│▓│▓│░│▓│░│▓│   ──────────►   │█│█│█│█│
   ├─┼─┼─┼─┼─┼─┼─┼─┤   coarse-grain  ├─┼─┼─┼─┤
   │░│▓│░│░│▓│░│▓│░│                 │█│█│█│█│
   ├─┼─┼─┼─┼─┼─┼─┼─┤                 ├─┼─┼─┼─┤
   │▓│░│▓│▓│░│▓│░│▓│                 │█│█│█│█│
   └─┴─┴─┴─┴─┴─┴─┴─┘                 └─┴─┴─┴─┘
    8×8, noisy, random              4×4, smooth, uniform

    After enough RG steps → fixed point: everywhere 1.0
    The skill is mastered. Coarse-graining changes nothing.
```

Each practice session is a block-spin transformation. You average over the noise of individual attempts and emerge with a smoother, more coherent performance. Repeat enough times and the lattice hits a **fixed point** — further practice doesn't change anything. That's mastery.

## The Critical Exponent Predicts Everything

Here's the part that should keep you up at night.

Every skill has a **correlation length exponent** ν that governs how fast it converges under coarse-graining. This single number determines the shape of the entire learning curve:

```
performance(t) = p_max − A · t^{−α}     where α = 1/ν
```

- **Small ν** (say 0.3) → **large α** → skill rockets toward mastery. Think: riding a bike.
- **Large ν** (say 3.0) → **small α** → skill crawls. Think: mastering chess.

The exponent is why some skills take 10 hours and others take 10,000. It's not that one is "harder" — it's that one has a steeper RG flow.

## Universality: The Deepest Cut

**All skills cluster into universality classes. Two completely different skills with the same critical exponent learn at the same rate.**

This is the Wilsonian insight that won Ken Wilson the 1982 Nobel Prize. In physics, systems that look nothing alike — magnets, fluids, polymers — can have *identical* critical behavior if they share the same symmetry and dimensionality. The microscopic details literally do not matter.

The same applies to learning:
- Typing and playing piano scales? Probably the same universality class (fine motor sequencing).
- Chess and Go? Likely the same class (strategic spatial reasoning).
- Chess and typing? Different classes, different exponents, completely different learning curves.

This means you can **transfer knowledge about learning rates between superficially unrelated skills**. If you know ν for one skill in a class, you know it for all of them.

## The Mathematical Framework

A skill is modeled as a 2D lattice where each site holds a skill level `s_i ∈ [0, 1]`. Under a block-spin RG transformation with block factor `b`:

```
s'(x,y) = (1/b²) Σ s(bx+dx, by+dy)
```

This is the real-space Wilsonian RG. After rescaling, a **fixed point** satisfies `s' = s` — the skill is invariant under further practice.

The **beta function** measures how the coupling (skill level) changes under RG:

```
β(g) = (g' − g) / ln(b)
```

At a fixed point, `β(g*) = 0`. The **correlation length** diverges near the novice→expert critical point:

```
ξ ~ |g − g_c|^{−ν}
```

This ν is the critical exponent that predicts everything.

## Applications

Once you see skills as RG flows, power laws appear *everywhere*:

### 🎮 Game Design — Skill Progression
A game designer fits power-law curves to player performance data. The critical exponent tells you whether a skill curve feels "snappy" (fast ν) or "grindy" (slow ν). Tune difficulty by tuning the effective ν — block size, practice variety, spacing.

### 🏢 Employee Onboarding
Track performance across tasks as a skill lattice. Compute ν for each role. Two positions with the same ν are in the same universality class — training materials transfer. Roles with wildly different ν need fundamentally different onboarding strategies.

### 🧠 Habit Formation
A habit is a skill lattice that has reached a fixed point. The question isn't "how do I build habits?" — it's "what's ν for this habit?" Small-ν habits form fast. Large-ν habits need environmental restructuring to effectively lower ν.

### 🤖 ML Learning Curves
Neural network training *is* a renormalization process — each layer coarse-grains the previous one. The test loss follows `L(t) ∝ t^{−α}` where α is set by the network's effective ν. Wide networks, narrow networks, transformers, CNNs — different architectures, potentially different universality classes.

## Quick Start

```bash
make
make test    # 30 tests, zero dependencies beyond libm
```

Requires: C11 compiler.

## API by Concept

### Build a Skill Lattice

```c
#include "renorm_learn.h"

// A 16×16 lattice of skill levels, each in [0, 1]
skill_lattice_t *skill = skill_lattice_create(16);
skill_lattice_random(skill, 42);        // random novice levels
```

### Coarse-Grain = Practice

```c
// One "practice session": average over 2×2 blocks
skill_lattice_t *after_practice = skill_coarse_grain(skill, 2);
```

### Check for Mastery (Fixed Point)

```c
if (skill_is_fixed_point(after_practice, 0.99)) {
    printf("Skill mastered. Further practice changes nothing.\n");
}
```

### Extract the Critical Exponent

```c
// Run multi-step RG flow and collect correlation lengths
rl_rgflow_t *flow = rgflow_create(skill, 2, 1.0, 5);
double lengths[5];
for (int i = 0; i < flow->nlevels; i++) {
    lengths[i] = /* correlation length at level i */;
}

double nu;
skill_critical_exponent(lengths, 5, 2, &nu);
printf("Critical exponent ν = %.3f\n", nu);
printf("Learning rate α = 1/ν = %.3f\n", 1.0 / nu);
```

### Predict Learning Curves

```c
// Fit performance data to p(t) = p_max − A · t^{−α}
double A, alpha, rmse;
fit_power_law(time_data, perf_data, n, p_max, &A, &alpha, &rmse);

// How long until 95% mastery?
double t_target;
predict_mastery_time(p_max, A, alpha, 0.95, &t_target);
printf("Reach 95%% mastery after %.0f practice units\n", t_target);
```

### Universality Class Comparison

```c
// Two different skills — do they learn at the same rate?
int same = skill_same_universality_class(nu_piano, nu_typing, 0.1);
if (same) {
    printf("Same universality class — transfer training materials!\n");
}
```

## Full API Reference

| Function | What It Does |
|---|---|
| `lattice_create(L)` | Create an L×L lattice |
| `lattice_random(lat, lo, hi, seed)` | Fill with random values |
| `lattice_magnetization(lat)` | Order parameter (mean skill) |
| `lattice_energy(lat)` | Nearest-neighbor interaction energy |
| `lattice_correlation(lat, r)` | Spatial correlation C(r) |
| `block_spin_coarse_grain(lat, b)` | RG transformation: b×b blocks → 1 site |
| `block_spin_rescale(lat, zeta)` | Rescale field amplitudes |
| `rg_step(lat, b, zeta)` | One complete RG step |
| `rgflow_create(lat, b, zeta, n)` | Multi-step RG flow |
| `compute_beta_function(g, lat, b, ...)` | β function β(g) |
| `find_fixed_point(lat0, b, ...)` | Iterate RG to convergence |
| `skill_coarse_grain(s, b)` | Learn by integrating practice |
| `skill_is_fixed_point(s, threshold)` | Check if mastered |
| `skill_critical_exponent(lengths, n, b, &nu)` | Extract ν from RG flow |
| `skill_same_universality_class(nu1, nu2, tol)` | Same class? |
| `fit_power_law(t, p, n, p_max, ...)` | Fit p = p_max − A·t^{−α} |
| `predict_mastery_time(p_max, a, alpha, target, &t)` | Time to reach target |
| `compare_learning_curves(α1, α2, tol)` | Similarity metric |

## Project Structure

```
renormalization-learning-c/
├── include/
│   └── renorm_learn.h      # Public API header
├── src/
│   ├── lattice.c           # 2D lattice system
│   ├── block_spin.c        # Wilsonian RG transformation
│   ├── skill_model.c       # Skill acquisition as RG
│   └── learning_curve.c    # Power-law learning curve prediction
├── tests/
│   └── test_renorm.c       # 30 tests
├── Makefile
└── README.md
```

## The Payoff

This library gives you three things:

1. **A language** for talking about learning that's precise enough to make predictions
2. **A single number** (ν) that determines the entire learning curve of a skill
3. **Universality** — the insight that superficially different skills can be *mathematically identical* in how they're learned

Once you see learning as renormalization, you can't unsee it. Every skill plateau is critical slowing down. Every breakthrough is a phase transition. Every expert is a fixed point.

## Author

Inspired by the [SuperInstance](https://github.com/SuperInstance) insight connecting Wilsonian RG to machine learning. Built in C11 with zero dependencies beyond libm.
