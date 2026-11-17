# Parse RXNMapper mapped reaction into canonical atom ids

**Deprecated:** temporarily unused. Prefer
[`rdkit_mcs_pairs`](https://drruili.github.io/MSCC/reference/rdkit_mcs_pairs.md)
/
[`rdkit_rxn_ids`](https://drruili.github.io/MSCC/reference/rdkit_rxn_ids.md)
(CLI) or a future atom-map API.

## Usage

``` r
RXNMapper_mapped_rxn_parse(rxn.result)
```

## Arguments

- rxn.result:

  Result from
  [`RXNMapper_map`](https://drruili.github.io/MSCC/reference/RXNMapper_map.md).

## Value

Nested list of map numbers named by canonical atom ids.
