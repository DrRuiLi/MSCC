# Atom IDs from a ChemmineR SDF

Atom IDs from a ChemmineR SDF

## Usage

``` r
# S4 method for class 'SDF'
atom(object, element = "ANY")
```

## Arguments

- object:

  An `SDF` object.

- element:

  Element symbol(s) to keep, or `"ANY"` for all non-H atoms.

## Value

Character vector of atom IDs (vertex names).
