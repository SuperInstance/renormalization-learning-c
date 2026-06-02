# Renormalization-Learning-C

**Wilsonian renormalization group as a model for agent skill acquisition.**

> *coarse-graining = learning*
> *fixed points = perfected skills*
> *critical exponents predict learning curves*

A C11 library that maps the mathematical structure of the Wilsonian RG onto the process of skill acquisition. The core insight (credit: [Step-3.5-Flash](https://github.com/SuperInstance)): just as a physicist coarse-grains a lattice to discover universal behavior, a learning agent integrates over microscopic practice details to converge on universal skill patterns.

## The Analogy

| Renormalization Group | Skill Acquisition |
|---|---|
| Spin field ϕ(x) on lattice | Skill level s_i across tasks / contexts |
| Coarse-grain over 2×2 blocks | Learn by integrating over microscopic practice details |
| Block average + rescale | Generalize from specific examples |
| RG fixed point ϕ* | Mastered skill (invariant under practice) |
| Correlation length exponent ν | Predicts learning rate (how fast skill generalizes) |
| Universality class | Skills that learn the same way |
| Beta function β(g) | How skill level changes with practice |
| Critical slowing down near T_c | Plateaus near expert-level transitions |
| Phase transition | Novice → Expert transition (power-law learning curve) |

## Mathematical Framework

A skill lattice is a 2D array of values in [0,1]. Under coarse-graining with block factor b:

```
s'(x,y) = (1/b²) Σ_{dx=0}^{b-1} Σ_{dy=0}^{b-1} s(bx+dx, by+dy)
```

This is the **Wilsonian RG transformation** in real space. After rescaling, the fixed point satisfies s' = s.

The **beta function** measures how the "skill coupling" g changes under RG:

```
β(g) = (g' - g) / ln(b)  
```

At a fixed point, β(g*) = 0.

The **correlation length** grows as the system approaches criticality:

```
ξ ~ |g - g_c|^{-ν}
```

This ν exponent predicts the learning curve: performance follows p(t) ∝ t^{-α} with α = 1/ν in the asymptotic regime.

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
│   └── test_renorm.c       # 30+ tests
├── Makefile
└── README.md
```

## Build & Test

```bash
make
make test
```

Requires: C11 compiler, libm.

## API Overview

### Lattice
- `lattice_create(L)` — L × L lattice of doubles
- `lattice_random()` — Uniform random initialization
- `lattice_magnetization()` — Average value (order parameter)
- `lattice_energy()` — Pairwise interaction energy
- `lattice_correlation(r)` — Spatial correlation function

### RG Transformations
- `block_spin_coarse_grain(lat, b)` — Block average (b × b → 1)
- `block_spin_rescale(lat, zeta)` — Rescale field amplitudes
- `rg_step(lat, b, zeta)` — One complete RG step
- `rgflow_create(lat, b, zeta, nsteps)` — Multi-step RG flow
- `compute_beta_function(g, lat, b, ...)` — Beta function β(g)
- `find_fixed_point(lat, b, ...)` — Iterate to fixed point

### Skill Model
- `skill_coarse_grain(s, b)` — Learn by integrating practice
- `skill_is_fixed_point(s, threshold)` — Check if mastered
- `skill_critical_exponent(lengths, n, b)` — Extract ν from RG flow
- `skill_same_universality_class(nu1, nu2, tol)` — Same class?
- `skill_phase_transition_magnetization(s)` — Order parameter

### Learning Curves
- `fit_power_law(t, p, n, p_max, ...)` — Fit p = p_max - A·t^{-α}
- `predict_mastery_time(p_max, a, alpha, target)` — Time to reach performance target
- `compare_learning_curves(alpha1, alpha2, tol)` — Similarity metric

## Author

Inspired by the [SuperInstance](https://github.com/SuperInstance) insight connecting Wilsonian RG to machine learning. Built in C11 with zero dependencies beyond libm.
