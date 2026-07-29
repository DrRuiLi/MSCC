# Max valences matching MassTools::getOption("MassTools.maxValences")
# (Kind & Fiehn allow N=5 in places; MassTools uses N=4 — we match MassTools.)
.golden_max_valences <- c(
  C = 4, H = 1, N = 4, O = 2, P = 5, S = 6,
  F = 1, Cl = 1, Br = 1, I = 1, Si = 4,
  B = 3, As = 5, Se = 6,
  Li = 1, Na = 1, K = 1, Mg = 2, Ca = 2,
  Fe = 3, Co = 4, Mn = 4, Cu = 2, Zn = 2
)

#' Kind & Fiehn / MassTools max element counts by mass bin (Rule #1)
#' @noRd
.golden_max_counts_limit <- function(mass) {
  Da <- abs(as.numeric(mass))
  if (is.na(Da) || Da < 500) {
    c(C = 29L, H = 72L, N = 10L, O = 18L, P = 4L, S = 7L, F = 15L, Cl = 8L, Br = 5L)
  } else if (Da < 1000) {
    c(C = 66L, H = 126L, N = 25L, O = 27L, P = 6L, S = 8L, F = 16L, Cl = 11L, Br = 8L)
  } else if (Da < 2000) {
    c(C = 115L, H = 236L, N = 32L, O = 63L, P = 6L, S = 8L, F = 16L, Cl = 11L, Br = 8L)
  } else if (Da < 3000) {
    c(C = 162L, H = 208L, N = 48L, O = 78L, P = 6L, S = 9L, F = 16L, Cl = 11L, Br = 8L)
  } else {
    NULL
  }
}

#' Senior's third theorem score (Rule #2)
#' @noRd
.golden_senior3 <- function(counts, valences = .golden_max_valences) {
  els <- names(counts)
  n <- as.numeric(counts)
  n[is.na(n)] <- 0
  v <- valences[els]
  v[is.na(v)] <- 0
  n_atoms <- sum(n)
  if (n_atoms <= 0) return(NA_real_)
  sum(n * v) - 2 * (n_atoms - 1)
}

#' HNOPS multi-element probability check (Rule #6); TRUE = pass
#' @noRd
.golden_element_heuristic <- function(counts) {
  get_n <- function(el) {
    if (el %in% names(counts)) as.numeric(counts[[el]]) else 0
  }
  r1 <- c(N = get_n("N"), O = get_n("O"), P = get_n("P"), S = get_n("S"))
  r2 <- r1[r1 > 1]
  if (length(r2) <= 2L) return(TRUE)
  # Match MassTools: r2 keeps order of N,O,P,S subset
  if (length(r2) == 4L &&
      (r2["N"] >= 10 || r2["O"] >= 20 || r2["P"] >= 4 || r2["S"] >= 3)) {
    return(FALSE)
  }
  key <- paste(names(r2), collapse = "")
  switch(
    key,
    NOP = {
      if (min(r2) > 3 && (r2["N"] >= 11 || r2["O"] >= 22 || r2["P"] >= 7)) return(FALSE)
    },
    OPS = {
      if (min(r2) > 1 && (r2["S"] >= 3 || r2["O"] >= 14 || r2["P"] >= 3)) return(FALSE)
    },
    NPS = {
      if (min(r2) > 1 && (r2["N"] >= 4 || r2["S"] >= 3 || r2["P"] >= 3)) return(FALSE)
    },
    NOS = {
      if (min(r2) > 6 && (r2["S"] >= 8 || r2["O"] >= 14 || r2["N"] >= 19)) return(FALSE)
    }
  )
  TRUE
}

#' Count vector for one formula row (uniso element symbols)
#' @noRd
.golden_counts_from_row <- function(row, colnames_mat) {
  els <- vapply(colnames_mat, get_ele_uniso, character(1))
  counts <- as.numeric(row)
  counts[is.na(counts)] <- 0
  # Aggregate isotope columns onto base elements (e.g. [13]C + C)
  out <- tapply(counts, els, sum)
  structure(as.numeric(out), names = names(out))
}

#' Monoisotopic mass from count named vector
#' @noRd
.golden_exact_mass <- function(counts) {
  et <- get_elem_table()
  et_uni <- et[!as.logical(et$is.isotope), , drop = FALSE]
  els <- names(counts)
  mono <- et_uni$mass[match(els, as.character(et_uni$symbol))]
  mono[is.na(mono)] <- 0
  sum(as.numeric(counts) * mono)
}

#' DBE from elem_table$unsaturation
#' @noRd
.golden_dbe <- function(counts) {
  et <- get_elem_table()
  et_uni <- et[!as.logical(et$is.isotope), , drop = FALSE]
  els <- names(counts)
  uns <- et_uni$unsaturation[match(els, as.character(et_uni$symbol))]
  uns[is.na(uns)] <- 0
  1 + sum(as.numeric(counts) * uns)
}

#' Check formulas against Kind & Fiehn seven golden rules (#1, #2, #4–#6)
#'
#' Reimplements the MassTools::calcMF post-filters (no MassTools dependency).
#' Rules #3 (isotope pattern) and #7 (TMS) are not applied.
#'
#' @param chemform Character vector of chemical formulas (neutral preferred).
#' @param mass Exact mass(es) for Rule #1 mass bins. Scalar recycled, or length
#'   equal to `chemform`. If `NULL`, monoisotopic mass is computed from the formula.
#' @param maxCounts Logical; apply element count caps by mass range (Rule #1).
#' @param SeniorRule Logical; apply Senior's third theorem (Rule #2).
#' @param HCratio Logical; apply H/C common range 0.2–3.1 (Rule #4).
#' @param moreRatios Logical; apply N/O/P/S vs C ratio caps (Rule #5).
#' @param elementHeuristic Logical; apply multi-HNOPS probability check (Rule #6).
#' @param DBErange Optional numeric length-2 `c(min, max)` for DBE filter; `NULL` disables.
#' @param parity Optional `"e"` or `"o"` (integer vs half-integer DBE); `NULL` disables.
#' @param SENIOR3_min Minimum Senior3 score to pass (default `0`).
#' @param return One of `"chemform"` (passing formulas), `"valid"` (logical vector),
#'   or `"all"` (diagnostic data.frame).
#'
#' @return Depends on `return`:
#'   - `"chemform"`: character vector of formulas that pass enabled rules
#'   - `"valid"`: logical vector, one per input
#'   - `"all"`: data.frame with input `chemform`, ratios/scores, per-rule flags, and `valid`
#'
#' @references
#' Kind T, Fiehn O (2007) Seven Golden Rules for heuristic filtering of molecular
#' formulas obtained by accurate mass spectrometry. BMC Bioinformatics 8:105.
#'
#' @seealso [chemform_decompose_mass()], [chemform_parse()]
#' @export
chemform_check_seven_golden_rules <- function(chemform,
                                              mass = NULL,
                                              maxCounts = TRUE,
                                              SeniorRule = TRUE,
                                              HCratio = TRUE,
                                              moreRatios = TRUE,
                                              elementHeuristic = TRUE,
                                              DBErange = NULL,
                                              parity = NULL,
                                              SENIOR3_min = 0,
                                              return = c("chemform", "valid", "all")) {
  return <- match.arg(return)
  chemform <- as.character(chemform)
  n <- length(chemform)

  empty_all <- function() {
    data.frame(
      chemform = character(),
      mass = numeric(),
      maxCounts = logical(),
      SENIOR3 = numeric(),
      SeniorRule = logical(),
      HtoC = numeric(),
      HCratio = logical(),
      NtoC = numeric(),
      OtoC = numeric(),
      PtoC = numeric(),
      StoC = numeric(),
      moreRatios = logical(),
      elementHeuristic = logical(),
      DBE = numeric(),
      DBE_ok = logical(),
      parity = character(),
      parity_ok = logical(),
      valid = logical(),
      stringsAsFactors = FALSE
    )
  }

  if (!n) {
    if (return == "chemform") return(character())
    if (return == "valid") return(logical())
    return(empty_all())
  }

  mat <- chemform_parse(chemform, return = "matrix")
  # Ensure row order matches input (chemform_parse uses formatted names as rownames)
  if (nrow(mat) != n) {
    stop("Internal error: parsed formula matrix row count does not match input.")
  }

  count_list <- lapply(seq_len(n), function(i) {
    .golden_counts_from_row(mat[i, , drop = TRUE], colnames(mat))
  })

  if (is.null(mass)) {
    mass_vec <- vapply(count_list, .golden_exact_mass, numeric(1))
  } else {
    mass_vec <- as.numeric(mass)
    if (length(mass_vec) == 1L) {
      mass_vec <- rep(mass_vec, n)
    } else if (length(mass_vec) != n) {
      stop("`mass` must be NULL, length 1, or length(chemform).")
    }
  }

  get_c <- function(counts, el) {
    if (el %in% names(counts)) as.numeric(counts[[el]]) else 0
  }

  maxCounts_ok <- rep(TRUE, n)
  SENIOR3 <- rep(NA_real_, n)
  SeniorRule_ok <- rep(TRUE, n)
  HtoC <- rep(NA_real_, n)
  HCratio_ok <- rep(TRUE, n)
  NtoC <- OtoC <- PtoC <- StoC <- rep(NA_real_, n)
  moreRatios_ok <- rep(TRUE, n)
  elementHeuristic_ok <- rep(TRUE, n)
  DBE <- rep(NA_real_, n)
  DBE_ok <- rep(TRUE, n)
  parity_val <- rep(NA_character_, n)
  parity_ok <- rep(TRUE, n)

  for (i in seq_len(n)) {
    counts <- count_list[[i]]
    C <- get_c(counts, "C")
    H <- get_c(counts, "H")
    N <- get_c(counts, "N")
    O <- get_c(counts, "O")
    P <- get_c(counts, "P")
    S <- get_c(counts, "S")

    # Rule #1
    lim <- .golden_max_counts_limit(mass_vec[i])
    if (is.null(lim)) {
      maxCounts_ok[i] <- TRUE
    } else {
      ok <- TRUE
      for (el in names(lim)) {
        n_el <- get_c(counts, el)
        if (n_el > lim[[el]]) {
          ok <- FALSE
          break
        }
      }
      maxCounts_ok[i] <- ok
    }

    # Rule #2
    SENIOR3[i] <- .golden_senior3(counts)
    SeniorRule_ok[i] <- !is.na(SENIOR3[i]) && SENIOR3[i] >= SENIOR3_min

    # Rule #4
    if (C > 0 && H > 0) {
      HtoC[i] <- H / C
      HCratio_ok[i] <- HtoC[i] > 0.2 && HtoC[i] < 3.1
    } else {
      HtoC[i] <- NA_real_
      HCratio_ok[i] <- TRUE
    }

    # Rule #5
    if (C > 0) {
      NtoC[i] <- N / C
      OtoC[i] <- O / C
      PtoC[i] <- P / C
      StoC[i] <- S / C
      moreRatios_ok[i] <-
        NtoC[i] < 1.3 && OtoC[i] < 1.2 && PtoC[i] < 0.3 && StoC[i] < 0.8
    } else {
      moreRatios_ok[i] <- TRUE
    }

    # Rule #6
    elementHeuristic_ok[i] <- .golden_element_heuristic(counts)

    # Optional DBE / parity
    DBE[i] <- .golden_dbe(counts)
    if (!is.null(DBErange) && length(DBErange) == 2L) {
      DBE_ok[i] <- DBE[i] >= DBErange[1] && DBE[i] <= DBErange[2]
    }
    # Integer DBE -> even-electron "e"; half-integer -> "o"
    parity_val[i] <- if (abs(DBE[i] - round(DBE[i])) < 1e-8) "e" else "o"
    if (!is.null(parity) && parity %in% c("e", "o")) {
      parity_ok[i] <- identical(parity_val[i], parity)
    }
  }

  valid <- rep(TRUE, n)
  if (isTRUE(maxCounts)) valid <- valid & maxCounts_ok
  if (isTRUE(SeniorRule)) valid <- valid & SeniorRule_ok
  if (isTRUE(HCratio)) valid <- valid & HCratio_ok
  if (isTRUE(moreRatios)) valid <- valid & moreRatios_ok
  if (isTRUE(elementHeuristic)) valid <- valid & elementHeuristic_ok
  if (!is.null(DBErange)) valid <- valid & DBE_ok
  if (!is.null(parity) && parity %in% c("e", "o")) valid <- valid & parity_ok

  out <- data.frame(
    chemform = chemform,
    mass = mass_vec,
    maxCounts = maxCounts_ok,
    SENIOR3 = SENIOR3,
    SeniorRule = SeniorRule_ok,
    HtoC = HtoC,
    HCratio = HCratio_ok,
    NtoC = NtoC,
    OtoC = OtoC,
    PtoC = PtoC,
    StoC = StoC,
    moreRatios = moreRatios_ok,
    elementHeuristic = elementHeuristic_ok,
    DBE = DBE,
    DBE_ok = DBE_ok,
    parity = parity_val,
    parity_ok = parity_ok,
    valid = valid,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (return == "chemform") return(chemform[valid])
  if (return == "valid") return(valid)
  out
}
