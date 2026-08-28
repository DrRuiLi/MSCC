# Molecular formula from SMILES

Molecular formula from SMILES

## Usage

``` r
get_smile_formula(smile)
```

## Arguments

- smile:

  Character SMILES string(s). Callers should unique the vector first
  when the same structure is repeated.

## Value

Character formula(s), formatted by
[`chemform_formate()`](https://drruili.github.io/MSCC/reference/chemform_formate.md).
