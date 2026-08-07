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

## 3. Internal steps

```text
target mass + ppm (+ mzabs)
        │
        ▼
build element alphabet (monoisotopic masses), sorted by mass
        │
        ▼
MCP / mass decomposition
  → all integer vectors n with |Σ n·m − mass| ≤ tol
        │
        ▼
optional chemical validity (DBE, nitrogen rule, …)
        │
        ▼
molecule objects (formula, exactmass, score, isotopes, …)
```

---

## 4. Brute force vs Rdisop

**Naive enumeration** (nested loops over atom counts) scales roughly as \(\prod_e (n_e^{\max}+1)\) — exponential in the number of element types when bounds grow with mass.

**Rdisop / DECOMP** does **not** do that. It uses efficient MCP algorithms (dynamic programming + smart backtracking; Böcker & Lipták):

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

MSCC provides two wrappers around the MCP enumeration (no isotope ranking):

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

- `decomposeMass` returns essentially **all** compositions for the allowed elements that fit the mass window.
- The engine is an **MCP algorithm**, not brute-force nested loops.
- Adding element types still hurts a lot: the **number of valid formulas** grows combinatorially, so CPU and memory grow with output size.
- A usable formula calculator = Rdisop enumeration + adduct/mass handling + chemical filters + ranking (isotope / DB).
