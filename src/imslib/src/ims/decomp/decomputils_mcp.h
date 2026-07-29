#ifndef IMS_DECOMPUTILS_MCP_H
#define IMS_DECOMPUTILS_MCP_H

#include <utility>

namespace ims {
namespace DecompUtils {

// Minimal MCP-only utilities needed by `RealMassDecomposer`.
// We intentionally avoid the full imslib `decomputils.h` dependency chain
// (which pulls in alphabet/printing helpers not required for formula enumeration).

template <typename DecompositionWeights>
std::pair<typename DecompositionWeights::alphabet_mass_type,
          typename DecompositionWeights::alphabet_mass_type>
getMinMaxWeightsRoundingErrors(const DecompositionWeights& weights) {
  using alphabet_mass_type = typename DecompositionWeights::alphabet_mass_type;
  using size_type = typename DecompositionWeights::size_type;

  alphabet_mass_type precision = weights.getPrecision();
  alphabet_mass_type minNegativeError = 0;
  alphabet_mass_type maxPositiveError = 0;

  for (size_type i = 0; i < weights.size(); ++i) {
    alphabet_mass_type error =
      (precision * static_cast<alphabet_mass_type>(weights.getWeight(i)) -
       weights.getAlphabetMass(i)) / weights.getAlphabetMass(i);

    if (error < 0 && error < minNegativeError) {
      minNegativeError = error;
    } else if (error > 0 && error > maxPositiveError) {
      maxPositiveError = error;
    }
  }

  return std::make_pair(minNegativeError, maxPositiveError);
}

template <typename DecompositionWeights, typename DecompositionType>
typename DecompositionWeights::alphabet_mass_type
getParentMass(const DecompositionWeights& weights, const DecompositionType& decomposition) {
  using alphabet_mass_type = typename DecompositionWeights::alphabet_mass_type;

  alphabet_mass_type parentMass = 0;
  typename DecompositionType::size_type counter = 0;

  for (typename DecompositionType::const_iterator pos = decomposition.begin();
       pos != decomposition.end(); ++pos) {
    parentMass += (*pos) * weights.getAlphabetMass(counter++);
  }

  return parentMass;
}

} // namespace DecompUtils
} // namespace ims

#endif // IMS_DECOMPUTILS_MCP_H

