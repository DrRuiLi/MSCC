# Export a Spectra object to a CFM-ID peak-list text file

Writes energy0/1/2 blocks (CE 10/20/40) for
[`CFM_annotate`](https://drruili.github.io/MSCC/reference/CFM.md).

## Usage

``` r
export_Spectra_peak_list_for_cfm(sp, file)
```

## Arguments

- sp:

  A `Spectra` object with `collisionEnergy`.

- file:

  Output file path.

## Value

Invisibly, `file`.
