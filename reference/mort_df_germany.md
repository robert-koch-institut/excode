# German all-cause mortality data

A dataset containing the number of weekly deaths in Germany reported to
the German Federal Statistical Office
<https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Sterbefaelle-Lebenserwartung/sterbefallzahlen.html>.
The variables are as follows:

## Usage

``` r
mort_df_germany
```

## Format

A tibble with 418 rows and 7 variables:

- date:

  Week of death

- observed:

  Number of deaths

- timepoint:

  A numeric column counting sequential weeks from 0 to total number of
  weeks, representing time progression for time-trend modeling

- sin52:

  Computed sinus component used for harmonic modeling

- cos52:

  Computed cosine component used for harmonic modeling

- bckg_week:

  A logical (\`TRUE\`/\`FALSE\`) variable marking background mortality
  weeks, used for estimating model parameters.
