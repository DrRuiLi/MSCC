# Annotate Isotopologues in Mass Spectra

Annotates mass spectra with isotopologue information by matching peaks
to expected isotope patterns. This function assigns fragment groups and
isotope counts to observed peaks based on mass tolerance.

Requires a `fragment_group` column on `peak_assignment` or
`fragment_define`. For MSIPAtomMap workflows, prefer
[`MSIP::Spectra_annotate_MSIPAtomMap()`](https://rdrr.io/pkg/MSIP/man/Spectra_annotate_MSIPAtomMap.html)
after
[`MSIP::MSIPAtomMap_get_FG_map()`](https://rdrr.io/pkg/MSIP/man/MSIPAtomMap_get_FG_map.html).

## Usage

``` r
CFM_annotate_isotopologues(
  sp,
  cfmd,
  iso_ele = "[13]C",
  iso_count = 0,
  ppm = 20
)
```

## Arguments

- sp:

  A Spectra object containing experimental mass spectra data

- cfmd:

  A CFM_data object containing fragment definitions and peak assignments
  (raw CFM slots only). `fragment_group` must already be present as a
  column.

- iso_ele:

  Isotope element specification, e.g., "\[13\]C" for carbon-13 (default:
  "\[13\]C")

- iso_count:

  Maximum number of isotope incorporations to consider (default: 0)

- ppm:

  Mass tolerance in parts per million for isotope matching (default: 20)

## Value

A data frame containing annotated spectrum data with fragment groups,
isotope counts, intensity ratios, and summed intensities for each
fragment group

## Details

Annotate Isotopologues in Mass Spectra

## Functions

- `CFM_annotate_isotopologues()`: annotate isotopologues
