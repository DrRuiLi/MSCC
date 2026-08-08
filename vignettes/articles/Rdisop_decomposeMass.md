# Rdisop mass decomposition (`decomposeMass`)

Quick reference for how Rdisop turns an exact mass (or isotope pattern) into candidate sum formulas, and how that fits a formula-calculator workflow.

Algorithms: Böcker et al. (DECOMP / Money Changing Problem). See Bioconductor vignette *Mass decomposition with the Rdisop package*.

---

## 1. High-level workflow (formula calculator)

```
m/z + polarity/adduct + ppm + element set
        ↓
1. Convert ion m/z → neutral (or charged) mass consistently
        ↓
2. Enumerate formulas (Rdisop::decomposeMass / decomposeIsotopes)
        ↓
3. Filter chemically (parity, H/C, heteroatom rules, DBE / Senior)
        ↓
4. Rank (ppm, RDBe, isotope fit if MS1 pattern available)
        ↓
5. Return tidy table of candidates
```

| Input | Role | Example |
|--------|------|---------|
| `mz` / `mass` | Measured ion m/z or neutral mass | `203.0526` |
| `ppm` / `mzabs` | Mass tolerance | `5` ppm |
| `z` | Charge (sign = polarity) | `+1` / `-1` |
| `adduct` | Optional ion type | `[M+H]+`, `[M-H]-` |
| `elements` | Allowed CHNOPS… | `C,H,N,O` |

Decide explicitly whether the input is **ion m/z** or **neutral mass**. Do not mix “neutral mass + `z = 1`” unless that is intentional.

In the MSdev / MSCC split:

| Layer | Responsibility |
|-------|----------------|
| **Rdisop** | Enumerate candidates from mass |
| **MSCC** | Adduct mass, formula format, isotope pattern |
| **MSdev** | Call calculator on features (`mzmed`), attach `candidate.formula` |

---

## 2. What `decomposeMass` does

`decomposeMass` is a thin wrapper:

```text
decomposeMass(mass, ppm, ...)
  → decomposeIsotopes(c(mass), intensity = 1, ...)
  → C++ Money-Changing / mass-decomposition solver
```

Core problem (Money Changing Problem / MCP):

> Find non-negative integers \((n_C, n_H, n_N, \ldots)\) such that  
> \(\sum n_e \cdot m_e\) lies within the mass tolerance window.

Typical R call:

```r
library(Rdisop)

elements <- initializeElements(c("C", "H", "N", "O"))
results <- decomposeMass(
  mass = 203.0526,
  ppm = 5,
  z = 1,
  elements = elements
)
getFormula(results)
```

Useful accessors: `getFormula()`, `getMass()`, `getScore()`, `getValid()`, `getIsotope()`.

With measured isotope peaks, prefer `decomposeIsotopes(masses, intensities, ...)` — same enumeration core, better ranking.

---

## 3. How MCP decomposes a mass into formulas

### 3.1 Integer problem (exact coins)

Classic Money Changing Problem: given coin denominations \(w_1,\ldots,w_k\) (here: **scaled** element masses) and target integer \(M\), find all non-negative integer vectors \(n = (n_1,\ldots,n_k)\) with

\[
\sum_{i=1}^{k} n_i\, w_i = M.
\]

In formula space that is “how many atoms of each element add up to this integer mass?”

### 3.2 Real masses → integers (what Rdisop / imslib actually do)

Monoisotopic masses are real numbers (e.g. \(m_\mathrm{C} = 12.000000\), \(m_\mathrm{H} = 1.007825\)). The pipeline is:

```text
element alphabet (mono masses m_i)
        │
        ▼
scale by precision (Rdisop default ~ 1e-5 Da)
  w_i = round(m_i / precision)   → integer Weights
        │
        ▼
optional divide-by-GCD (shrink integers without changing solutions)
        │
        ▼
RealMassDecomposer
  1. expand tolerance window for rounding error
  2. for each integer mass M in [M_lo, M_hi):
       IntegerMassDecomposer → all n with Σ n·w = M
  3. keep only those with |Σ n·m_i − mass| ≤ abs_error
        │
        ▼
formula strings + exact masses (real Σ n·m_i)
```

Tolerance on the real mass:

\[
\texttt{abs\_error} = \texttt{ppm} \cdot \texttt{mass} \cdot 10^{-6} + \texttt{mzabs}
\]

`RealMassDecomposer` widens the integer search range slightly using precomputed min/max weight rounding errors, then **re-checks** each candidate against the true real mass so false positives from scaling are dropped.

### 3.3 Integer engine (Böcker & Lipták residue table)

`IntegerMassDecomposer` (paper: *Efficient Mass Decomposition*, Böcker & Lipták, ACM SAC-BIO 2004) does **not** store a classical DP table of size \(O(k \cdot M)\).

Idea in one line: let \(a = \min_i w_i\) (smallest alphabet weight). Precompute, for each residue class modulo \(a\), the **smallest** mass of that residue that is still decomposable, plus witness info for reconstructing atom counts. Querying a mass \(M\) then becomes:

1. Look up residue \(M \bmod a\).
2. Decide existence / recover one or all decompositions by walking witnesses and branching over the remaining alphabet (backtracking guided by the residue structure).

Consequences:

| Property | Meaning for formula finding |
|----------|-----------------------------|
| Output-sensitive | Cost tracks **how many** formulas hit the window, not \(\prod (n_i^{\max}+1)\) |
| Memory | Residue table ≪ full knapsack DP for typical MS masses |
| Completeness | Returns **all** non-negative compositions for each integer \(M\) in range |
| Alphabet order | Weights are fixed for a cached decomposer; changing elements rebuilds the table |

### 3.4 Tiny worked sketch

Alphabet `{H, C}` (conceptually), target ≈ \(16.0313\) Da (methane), tight `abs_error`.

1. Scale \(m_\mathrm{H}, m_\mathrm{C}\) → integers \(w_H, w_C\); build residue table mod \(\min(w_H,w_C)\).
2. Convert target ± error → integer interval \([M_\mathrm{lo}, M_\mathrm{hi})\).
3. For each \(M\) in that interval, enumerate \(n_H, n_C \ge 0\) with \(n_H w_H + n_C w_C = M\).
4. Keep vectors whose **real** mass \(\lvert n_H m_H + n_C m_C - 16.0313\rvert \le\) error → e.g. `CH4`.

Same steps for CHNOPS; only the alphabet (and thus the residue table) grows.

### 3.5 Pipeline position (still after MCP)

```text
target mass + ppm (+ mzabs)
        │
        ▼
MCP enumeration (sections 3.1–3.4)
  → all integer vectors n with |Σ n·m − mass| ≤ tol
        │
        ▼
optional chemical validity (DBE, nitrogen rule, …)
        │
        ▼
molecule objects / tidy table (formula, exactmass, score, …)
```

MCP only answers “what compositions match the mass window?” Chemical filters and isotope ranking are separate stages.

### 3.6 How MSCC supports negative element counts

Classic MCP (and Rdisop) only allow **non-negative** atom counts. That is fine for molecules, but **not** for chemical replacements / mass differences, e.g.

| Change | Difference formula | Meaning |
|--------|--------------------|---------|
| COOH → COONa | `H-1Na` | lose 1 H, gain 1 Na |
| observed Δmass | \(M = m_\mathrm{Na} - m_\mathrm{H} \approx 21.98194\) Da | mass to decompose |

You want ranges such as \(n_H \in [-3,\, n_H^{\max}]\), not only \(n_H \ge 0\).

#### Why we do not change the MCP kernel

Vendored imslib `IntegerMassDecomposer` / `RealMassDecomposer` enumerate \(n'_i \ge 0\) only (unsigned counts, coin problem). Rewriting that for signed coins is unnecessary: any signed box \([\min_i, \max_i]\) reduces to a non-negative problem by a **mass offset**.

#### Variable substitution (offset)

For each element \(i\):

\[
n_i = n'_i + \min_i,\qquad 0 \le n'_i \le \max_i - \min_i
\]

Mass identity:

\[
\sum_i n_i\, m_i
= \sum_i n'_i\, m_i + \sum_i \min_i\, m_i
\]

so MCP is run on the **shifted** target

\[
M' = M - \sum_i \min_i\, m_i
\]

with relative bounds \(n'_i \in [0,\, \max_i - \min_i]\). After enumeration, recover \(n_i = n'_i + \min_i\) and emit signed formula strings (`H-1Na`).

```text
user: mass M, min_i (may be < 0), max_i
        │
        ▼
mass_offset = Σ min_i * m_i
M' = M - mass_offset
relative max = max_i - min_i
        │
        ▼
RealMassDecomposer.getDecompositions(M', abs_error)
  → only n'_i ≥ 0
        │
        ▼
filter n'_i ≤ max_i - min_i
n_i = n'_i + min_i          # can be negative
exactmass = Σ n_i * m_i
formula   = counts_to_formula(n)   # e.g. "H-1Na"
```

Implemented in `src/mcp_decompose.cpp` (`mcp_decompose_mass`); R API is `chemform_decompose_mass(..., min_elements = ...)`.

Default `min_elements = NULL` → all mins `0`, so ordinary molecule decompositions are unchanged.

#### Worked example: COOH → COONa (`H-1Na`)

```r
dM <- chemform_mz("Na") - chemform_mz("H")   # ≈ 21.98194
chemform_decompose_mass(
  mass = dM,
  elements = c("H", "Na"),
  min_elements = c(H = -3),   # or "H-3"
  max_elements = c(H = 3, Na = 3),
  check_rule = FALSE
)
# → formula H-1Na, exactmass ≈ 21.98194, ppm ≈ 0
```

Numbers for that hit:

| Quantity | Value |
|----------|--------|
| \(\min_H, \min_{Na}\) | \(-3\), \(0\) |
| \(\mathrm{mass\_offset}\) | \(-3 \cdot m_H \approx -3.023\) |
| \(M'\) | \(21.982 - (-3.023) \approx 25.005\) |
| relative solution \(n'\) | \(n'_H = 2\), \(n'_{Na} = 1\) |
| recovered \(n\) | \(n_H = 2 + (-3) = -1\), \(n_{Na} = 1\) |
| check | \((-1)m_H + m_{Na} = M\) |

#### API notes

- **Input:** named vector `c(H = -3)` or formula string `"H-3"` / `"C0H-3"`.
- **Output:** formulas keep the minus sign (`H-1`); count `1` is still omitted for positive atoms (`Na` not `Na1`).
- **`check_rule`:** Seven Golden Rules assume molecular formulas. If any `min_elements < 0`, MSCC turns `check_rule` off (with a message) even if you pass `TRUE`.
- **Target mass:** use the **observed Δmass** of the replacement (here \(m_{Na}-m_H\)), not a negative “molecular” mass.
- Rdisop `decomposeMass` itself does **not** expose this; negative mins are an MSCC extension around the same MCP core.

---

## 4. Brute force vs Rdisop / MSCC MCP

**Naive enumeration** (nested loops over atom counts) scales roughly as \(\prod_e (n_e^{\max}+1)\) — exponential in the number of element types when bounds grow with mass.

**Rdisop / DECOMP / MSCC’s vendored imslib** do **not** do that. They use the residue-table MCP above:

- Finds **all** compositions whose mass falls in the tolerance window
- Runtime is typically closer to **proportional to the number of solutions**, not to a full nested-loop lattice walk
- Far better memory profile than classical full DP tables for this problem

So: **yes**, it enumerates all mass-matching compositions for the allowed elements; **no**, it is not a dumb “try every combination” search.

---

## 5. Does adding elements explode the cost?

**Yes in practice** — mainly because the **solution space** grows, not because the search is naive.

| Factor | Effect |
|--------|--------|
| More element species | Larger alphabet → more ways to hit the same mass → more candidates |
| Higher mass | More atoms possible → more solutions |
| Wider `ppm` / `mzabs` | Wider window → more solutions |
| Algorithm | MCP-style; pay per solution + DP overhead, not \(k^{n}\) nested loops |

Rdisop’s vignette notes that the result list grows with mass, allowed ppm, and the allowed elements list.

---

## 6. Practical controls

1. Keep the alphabet small (CHNO or CHNOPS first).
2. Tighten `ppm` / `mzabs`.
3. Use `minElements` / `maxElements` to cap atom counts.
4. Prefer `decomposeIsotopes` when M+1 / M+2 intensities exist.
5. Post-filter (H/C ranges, RDBE, nitrogen rule) and keep top-\(N\) by ppm / score.
6. Expect a sharp jump in candidates when adding Cl, Br, metals, or many heteroatoms.

---
 
## 7. Path A in MSCC (formula-only MCP calculator)

MSCC vendors the same imslib MCP stack (see §3) behind two wrappers (no isotope ranking):

| Function | Input | Role |
|----------|--------|------|
| `chemform_decompose_mass()` | Neutral exact mass | Thin MCP wrapper only |
| `chemform_decompose_mz()` | Ion m/z + `charge` | Convert m/z → neutral, then call `chemform_decompose_mass()` |

Charge conversion (same electron-mass convention as `chemform_mz()`):

```text
M_neutral = mz * abs(charge) + e * charge   # when charge != 0
```

Suggested parity check (formula sets, ignoring isotope ranking):

```r
library(Rdisop)
library(MSCC)

els <- initializeElements(c("C","H","N","O","P"))

r <- decomposeMass(
  mass = 203.0526,
  ppm = 5,
  elements = els,
  z = 0,
  maxisotopes = 1
)$formula

m <- chemform_decompose_mass(
  mass = 203.0526,
  ppm = 5,
  elements = c("C","H","N","O","P")
)$formula

stopifnot(setequal(r, m))

# ion m/z path
chemform_decompose_mz(mz = 204.0599, charge = 1, ppm = 5)

# replacement / difference formula (signed mins; see §3.6)
dM <- chemform_mz("Na") - chemform_mz("H")
chemform_decompose_mass(
  mass = dM,
  elements = c("H", "Na"),
  min_elements = c(H = -3),
  max_elements = c(H = 3, Na = 3),
  check_rule = FALSE
)
```

---

## 8. Suggested return table

After decomposition, tidy to something like:

| Column | Meaning |
|--------|---------|
| `formula` | Sum formula string |
| `exactmass` | Calculated exact mass |
| `ppm` | \((exactmass - target) / target \times 10^6\) |
| `score` | Rdisop / isotope score |
| `DBE` / `valid` | Chemical plausibility hints |

---

## 9. Bottom line

- `decomposeMass` / `chemform_decompose_mass` return essentially **all** compositions for the allowed elements that fit the mass window.
- The engine is an **MCP algorithm** (scale → integer residue-table decompose → real-mass filter), not brute-force nested loops.
- Adding element types still hurts a lot: the **number of valid formulas** grows combinatorially, so CPU and memory grow with output size.
- A usable formula calculator = MCP enumeration + adduct/mass handling + chemical filters + ranking (isotope / DB).
- Signed `min_elements` (MSCC) cover replacement deltas via mass offset; leave `check_rule` off for those.