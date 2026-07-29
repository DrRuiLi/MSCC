// MCP-only (formula enumeration) wrapper around vendored imslib `RealMassDecomposer`.
// This intentionally skips isotope distribution / scoring.

#include <Rcpp.h>

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <unordered_map>

#include <ims/weights.h>
#include <ims/decomp/realmassdecomposer.h>

using namespace Rcpp;

namespace {

// Build a stable cache key from element order + mono masses.
std::string make_cache_key(const CharacterVector& element_names, const NumericVector& mono_masses) {
  std::ostringstream oss;
  oss << std::setprecision(17);
  for (int i = 0; i < element_names.size(); ++i) {
    oss << as<std::string>(element_names[i]) << ':';
    oss << mono_masses[i] << ';';
  }
  return oss.str();
}

// Convert a count vector into a simple formula string like "C2H6O".
std::string counts_to_formula(const std::vector<std::string>& element_names, const std::vector<unsigned long long>& counts) {
  std::string out;
  out.reserve(16 + 8 * element_names.size());
  for (size_t i = 0; i < element_names.size(); ++i) {
    const auto c = counts[i];
    if (c == 0) continue;
    out += element_names[i];
    if (c != 1ULL) out += std::to_string(c);
  }
  return out;
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List mcp_decompose_mass(double mass,
                              double abs_error,
                              Rcpp::NumericVector mono_masses,
                              Rcpp::CharacterVector element_names,
                              Rcpp::IntegerVector min_counts,
                              Rcpp::IntegerVector max_counts) {
  if (mono_masses.size() == 0) {
    return List::create(_["formula"] = CharacterVector(),
                        _["exactmass"] = NumericVector());
  }
  const int n = mono_masses.size();

  if (element_names.size() != n || min_counts.size() != n || max_counts.size() != n) {
    stop("mcp_decompose_mass: size mismatch among inputs (mono_masses, element_names, min_counts, max_counts).");
  }

  for (int i = 0; i < n; ++i) {
    if (min_counts[i] > max_counts[i]) {
      stop("mcp_decompose_mass: min_counts must be <= max_counts for each element.");
    }
  }

  // Prepare element order + counts bounds.
  std::vector<std::string> el_names;
  el_names.reserve(n);
  for (int i = 0; i < n; ++i) {
    el_names.push_back(as<std::string>(element_names[i]));
  }

  std::vector<unsigned long long> min_c(n), max_c(n);
  for (int i = 0; i < n; ++i) {
    min_c[i] = static_cast<unsigned long long>(min_counts[i]);
    max_c[i] = static_cast<unsigned long long>(max_counts[i]);
  }

  // Cached decomposer for this alphabet.
  static std::unordered_map<std::string, std::unique_ptr<ims::RealMassDecomposer>> decomposer_cache;

  const double precision = 1.0e-5; // matches Rdisop default precision
  const std::string cache_key = make_cache_key(element_names, mono_masses);

  auto it = decomposer_cache.find(cache_key);
  if (it == decomposer_cache.end()) {
    std::vector<double> masses(n);
    for (int i = 0; i < n; ++i) masses[i] = mono_masses[i];

    ims::Weights weights(masses, precision);
    weights.divideByGCD();
    decomposer_cache.emplace(cache_key, std::make_unique<ims::RealMassDecomposer>(weights));
    it = decomposer_cache.find(cache_key);
  }

  // Enumerate all compositions within abs_error.
  auto decomps = it->second->getDecompositions(mass, abs_error);

  std::vector<std::string> formulas;
  std::vector<double> exactmasses;
  formulas.reserve(decomps.size());
  exactmasses.reserve(decomps.size());

  // Filter by per-element min/max counts and stringify.
  for (const auto& decomposition : decomps) {
    bool ok = true;
    for (int i = 0; i < n; ++i) {
      const unsigned long long c = static_cast<unsigned long long>(decomposition[i]);
      if (c < min_c[i] || c > max_c[i]) {
        ok = false;
        break;
      }
    }
    if (!ok) continue;

    double em = 0.0;
    std::vector<unsigned long long> counts(n);
    for (int i = 0; i < n; ++i) {
      const unsigned long long c = static_cast<unsigned long long>(decomposition[i]);
      counts[i] = c;
      em += static_cast<double>(c) * mono_masses[i];
    }

    formulas.push_back(counts_to_formula(el_names, counts));
    exactmasses.push_back(em);
  }

  return List::create(_["formula"] = formulas,
                       _["exactmass"] = exactmasses);
}

