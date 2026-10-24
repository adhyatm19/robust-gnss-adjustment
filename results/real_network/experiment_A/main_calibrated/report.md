# Experiment A: Madrid real network (main_calibrated)

- network: `results/real_network/networks/main_calibrated`
- reference coordinates: `results/real_network/networks/main_calibrated/reference_coordinates.csv`

- stations: 4, baselines: 9, observations: 27, parameters: 9
- datum: station `MADR` held fixed; rank(A) = 9; redundancy r = 18; cond(N) = 28.8
- a-priori variance factor sigma0^2 = 1

## Global model test and method outcomes

| method | T = v'Pv/sigma0^2 | critical | passed | sigma0_hat^2 | flagged | iterations |
|---|---|---|---|---|---|---|
| WLS | 108.13 | 28.87 | False | 6.007 | - | 1 |
| DIA-reject | 10.22 | 16.92 | True | 1.136 | [1, 8, 9] | 4 |
| DIA-inflate | 11.92 | 28.87 | True | 0.662 | [1, 8, 9] | 4 |
| Huber | 64.51 | 28.87 | False | 3.584 | [8] | 15 |
| Hampel | 57.05 | 28.87 | False | 3.169 | [8] | 24 |
| SD-reject | 10.22 | 16.92 | True | 1.136 | [1, 8, 9] | 4 |

Statistical flags indicate inconsistency with the stochastic model; they do not by themselves prove a physical observation error.

## Coordinate differences to the reference (relative to the datum station)

| method      | station   |   diff_e_mm |   diff_n_mm |   diff_u_mm |   sigma_e_mm |   sigma_n_mm |   sigma_u_mm |
|:------------|:----------|------------:|------------:|------------:|-------------:|-------------:|-------------:|
| WLS         | CEBR      |       11.85 |       -2.47 |        5.91 |         1.61 |         0.74 |         1.76 |
| WLS         | VILL      |       -4.19 |        0.87 |        6.14 |         1.58 |         0.73 |         1.76 |
| WLS         | YEBE      |        0.78 |       -2.28 |        0.95 |         1.97 |         0.90 |         2.20 |
| DIA-reject  | CEBR      |       10.19 |       -2.68 |       13.26 |         2.16 |         0.98 |         2.36 |
| DIA-reject  | VILL      |       -5.56 |        2.53 |       -0.95 |         2.09 |         0.96 |         2.33 |
| DIA-reject  | YEBE      |       -0.31 |       -1.77 |        0.92 |         2.14 |         0.97 |         2.39 |
| DIA-inflate | CEBR      |       10.23 |       -2.67 |       13.13 |         2.15 |         0.97 |         2.35 |
| DIA-inflate | VILL      |       -5.54 |        2.50 |       -0.83 |         2.08 |         0.96 |         2.32 |
| DIA-inflate | YEBE      |       -0.29 |       -1.78 |        0.92 |         2.14 |         0.97 |         2.39 |
| Huber       | CEBR      |       11.47 |       -2.74 |        7.98 |         1.78 |         0.81 |         1.95 |
| Huber       | VILL      |       -4.42 |        1.43 |        3.89 |         1.75 |         0.81 |         1.95 |
| Huber       | YEBE      |        0.56 |       -2.17 |        0.85 |         2.00 |         0.91 |         2.24 |
| Hampel      | CEBR      |       11.40 |       -2.84 |        9.14 |         1.80 |         0.82 |         1.97 |
| Hampel      | VILL      |       -4.53 |        1.60 |        2.74 |         1.77 |         0.82 |         1.97 |
| Hampel      | YEBE      |        0.48 |       -2.15 |        0.83 |         2.01 |         0.92 |         2.25 |
| SD-reject   | CEBR      |       10.19 |       -2.68 |       13.26 |         2.16 |         0.98 |         2.36 |
| SD-reject   | VILL      |       -5.56 |        2.53 |       -0.95 |         2.09 |         0.96 |         2.33 |
| SD-reject   | YEBE      |       -0.31 |       -1.77 |        0.92 |         2.14 |         0.97 |         2.39 |

## WLS residuals and local tests (ENU of from-station)

|   baseline_id | from   | to   |   length_m |   v_e_mm |   v_n_mm |   v_u_mm |   T_b |   T_b_critical |   local_redundancy |   sd_statistic |   sd_critical |   sd_dir_elevation_deg |
|--------------:|:-------|:-----|-----------:|---------:|---------:|---------:|------:|---------------:|-------------------:|---------------:|--------------:|-----------------------:|
|             1 | MADR   | CEBR |   1.04e+04 |   -1.6   |   0.867  |     8.82 | 12.7  |           16.3 |               2.11 |           3.57 |          4.03 |                   78.3 |
|             2 | MADR   | VILL |   2.53e+04 |   -2.07  |  -0.0979 |     6.65 |  6.7  |           16.3 |               2.11 |           2.59 |          4.03 |                   72.3 |
|             3 | VILL   | YEBE |   7.38e+04 |   -2.82  |   3.1    |    -5.44 | 14.6  |           16.3 |               1.64 |           3.82 |          4.03 |                   52   |
|             4 | CEBR   | VILL |   3.53e+04 |    0.486 |  -1.38   |    16.5  | 37.5  |           16.3 |               2.14 |           6.12 |          4.03 |                   84.4 |
|             5 | MADR   | CEBR |   1.04e+04 |    3.66  |  -0.107  |    -5.86 |  5.75 |           16.3 |               2.17 |           2.4  |          4.03 |                   58.7 |
|             6 | MADR   | YEBE |   9.9e+04  |    2.9   |  -1.78   |    -1.55 |  5.39 |           16.3 |               1.78 |           2.32 |          4.03 |                   22   |
|             7 | CEBR   | YEBE |   1.09e+05 |    0.612 |  -1.51   |     6.84 |  8.31 |           16.3 |               1.73 |           2.88 |          4.03 |                   77.2 |
|             8 | CEBR   | VILL |   3.53e+04 |   -1.01  |   3.62   |   -19.3  | 54.3  |           16.3 |               2.13 |           7.37 |          4.03 |                   79.3 |
|             9 | MADR   | VILL |   2.53e+04 |   -2.62  |   1.25   |    -9.75 | 12.6  |           16.3 |               2.17 |           3.55 |          4.03 |                   74.1 |
