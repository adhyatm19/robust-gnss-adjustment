# Experiment A: Madrid real network (quiet_formal)

- network: `results/real_network/networks/quiet_formal`
- reference coordinates: `results/real_network/networks/quiet_formal/reference_coordinates.csv`

- stations: 4, baselines: 9, observations: 27, parameters: 9
- datum: station `MADR` held fixed; rank(A) = 9; redundancy r = 18; cond(N) = 28.9
- a-priori variance factor sigma0^2 = 1

## Global model test and method outcomes

| method | T = v'Pv/sigma0^2 | critical | passed | sigma0_hat^2 | flagged | iterations |
|---|---|---|---|---|---|---|
| WLS | 12.37 | 28.87 | True | 0.687 | - | 1 |
| DIA-reject | 12.37 | 28.87 | True | 0.687 | - | 1 |
| DIA-inflate | 12.37 | 28.87 | True | 0.687 | - | 1 |
| Huber | 12.37 | 28.87 | True | 0.687 | - | 2 |
| Hampel | 12.37 | 28.87 | True | 0.687 | - | 2 |
| SD-reject | 12.37 | 28.87 | True | 0.687 | - | 1 |

Statistical flags indicate inconsistency with the stochastic model; they do not by themselves prove a physical observation error.

## Coordinate differences to the reference (relative to the datum station)

| method      | station   |   diff_e_mm |   diff_n_mm |   diff_u_mm |   sigma_e_mm |   sigma_n_mm |   sigma_u_mm |
|:------------|:----------|------------:|------------:|------------:|-------------:|-------------:|-------------:|
| WLS         | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| WLS         | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| WLS         | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |
| DIA-reject  | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| DIA-reject  | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| DIA-reject  | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |
| DIA-inflate | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| DIA-inflate | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| DIA-inflate | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |
| Huber       | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| Huber       | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| Huber       | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |
| Hampel      | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| Hampel      | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| Hampel      | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |
| SD-reject   | CEBR      |       11.60 |       -3.95 |        5.97 |         1.82 |         0.85 |         2.02 |
| SD-reject   | VILL      |       -4.95 |       -0.17 |        4.25 |         1.80 |         0.84 |         2.03 |
| SD-reject   | YEBE      |        3.09 |       -3.62 |        7.72 |         2.22 |         1.03 |         2.52 |

## WLS residuals and local tests (ENU of from-station)

|   baseline_id | from   | to   |   length_m |   v_e_mm |   v_n_mm |   v_u_mm |   T_b |   T_b_critical |   local_redundancy |   sd_statistic |   sd_critical |   sd_dir_elevation_deg |
|--------------:|:-------|:-----|-----------:|---------:|---------:|---------:|------:|---------------:|-------------------:|---------------:|--------------:|-----------------------:|
|             1 | MADR   | CEBR |   1.04e+04 |    1.59  |   0.259  |    -1.02 | 0.408 |           16.3 |               2.15 |          0.638 |          4.03 |                   32.5 |
|             2 | MADR   | VILL |   2.53e+04 |   -3.18  |   1.16   |    -4.88 | 3.94  |           16.3 |               2.15 |          1.99  |          4.03 |                   55.1 |
|             3 | VILL   | YEBE |   7.38e+04 |   -2.08  |  -0.853  |     1.44 | 1.33  |           16.3 |               1.68 |          1.15  |          4.03 |                   32.6 |
|             4 | CEBR   | VILL |   3.53e+04 |   -0.349 |  -0.381  |    -1.41 | 0.396 |           16.3 |               2.12 |          0.629 |          4.03 |                   69.8 |
|             5 | MADR   | CEBR |   1.04e+04 |    3.24  |  -1.25   |    -1.73 | 2.53  |           16.3 |               2.15 |          1.59  |          4.03 |                   26.4 |
|             6 | MADR   | YEBE |   9.9e+04  |   -2.25  |   0.845  |     1.9  | 1.86  |           16.3 |               1.76 |          1.37  |          4.03 |                   37.7 |
|             7 | CEBR   | YEBE |   1.09e+05 |    4.44  |   0.0271 |    -3.28 | 4.18  |           16.3 |               1.71 |          2.04  |          4.03 |                   36.6 |
|             8 | CEBR   | VILL |   3.53e+04 |    0.609 |  -0.506  |     1.97 | 0.531 |           16.3 |               2.13 |          0.729 |          4.03 |                   68.3 |
|             9 | MADR   | VILL |   2.53e+04 |    0.805 |  -1.13   |     5.77 | 3.51  |           16.3 |               2.15 |          1.87  |          4.03 |                   76.4 |
