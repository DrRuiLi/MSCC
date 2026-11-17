# Map reaction SMILES with RXNMapper

**Deprecated:** temporarily unused. Prefer
[`rdkit_mcs_pairs`](https://drruili.github.io/MSCC/reference/rdkit_mcs_pairs.md)
/
[`rdkit_rxn_ids`](https://drruili.github.io/MSCC/reference/rdkit_rxn_ids.md)
(CLI) or a future atom-map API.

## Usage

``` r
RXNMapper_map(from.smiles, to.smiles, RXNMapper = MSCC::get_RXNMapper())
```

## Arguments

- from.smiles:

  Substrate SMILES (character vector).

- to.smiles:

  Product SMILES (character vector).

- RXNMapper:

  Callable from
  [`get_RXNMapper`](https://drruili.github.io/MSCC/reference/get_RXNMapper.md).

## Value

RXNMapper result list, or `NA` if the reaction is too long.
