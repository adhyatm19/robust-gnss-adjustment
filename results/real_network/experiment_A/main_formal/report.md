# Experiment A: Madrid real network (main_formal)

- network: `results/real_network/networks/main_formal`
- reference coordinates: `results/real_network/networks/main_formal/reference_coordinates.csv`

- stations: 4, baselines: 9, observations: 27, parameters: 9
- datum: station `MADR` held fixed; rank(A) = 9; redundancy r = 18; cond(N) = 28.8
- a-priori variance factor sigma0^2 = 1

## Global model test and method outcomes

| method | T = v'Pv/sigma0^2 | critical | passed | sigma0_hat^2 | flagged | iterations |
|---|---|---|---|---|---|---|
| WLS | 79.29 | 28.87 | False | 4.405 | - | 1 |
| DIA-reject | 20.87 | 21.03 | True | 1.739 | [8, 9] | 3 |
| DIA-inflate | 21.85 | 28.87 | True | 1.214 | [8, 9] | 3 |
| Huber | 54.10 | 28.87 | False | 3.006 | [8] | 15 |
| Hampel | 53.72 | 28.87 | False | 2.985 | [8] | 20 |
| SD-reject | 20.87 | 21.03 | True | 1.739 | [8, 9] | 3 |

Statistical flags indicate inconsistency with the stochastic model; they do not by themselves prove a physical observation error.

## Coordinate differences to the reference (relative to the datum station)

| method      | station   |   diff_e_mm |   diff_n_mm |   diff_u_mm |   sigma_e_mm |   sigma_n_mm |   sigma_u_mm |
|:------------|:----------|------------:|------------:|------------:|-------------:|-------------:|-------------:|
| WLS         | CEBR      |       11.85 |       -2.47 |        5.91 |         1.88 |         0.86 |         2.06 |
| WLS         | VILL      |       -4.19 |        0.87 |        6.14 |         1.85 |         0.86 |         2.06 |
| WLS         | YEBE      |        0.78 |       -2.28 |        0.95 |         2.30 |         1.05 |         2.57 |
| DIA-reject  | CEBR      |       11.40 |       -2.87 |        7.65 |         2.03 |         0.93 |         2.22 |
| DIA-reject  | VILL      |       -5.01 |        2.43 |       -3.81 |         2.33 |         1.07 |         2.59 |
| DIA-reject  | YEBE      |        0.38 |       -1.86 |       -2.02 |         2.38 |         1.09 |         2.66 |
| DIA-inflate | CEBR      |       11.41 |       -2.86 |        7.63 |         2.03 |         0.93 |         2.22 |
| DIA-inflate | VILL      |       -5.00 |        2.41 |       -3.63 |         2.32 |         1.07 |         2.58 |
| DIA-inflate | YEBE      |        0.39 |       -1.87 |       -1.96 |         2.38 |         1.09 |         2.66 |
| Huber       | CEBR      |       11.66 |       -2.74 |        7.31 |         1.99 |         0.91 |         2.18 |
| Huber       | VILL      |       -4.26 |        1.26 |        4.52 |         1.97 |         0.91 |         2.19 |
| Huber       | YEBE      |        0.68 |       -2.23 |        0.85 |         2.31 |         1.06 |         2.59 |
| Hampel      | CEBR      |       11.66 |       -2.75 |        7.39 |         2.00 |         0.91 |         2.18 |
| Hampel      | VILL      |       -4.27 |        1.27 |        4.44 |         1.97 |         0.91 |         2.19 |
| Hampel      | YEBE      |        0.67 |       -2.23 |        0.85 |         2.31 |         1.06 |         2.59 |
| SD-reject   | CEBR      |       11.40 |       -2.87 |        7.65 |         2.03 |         0.93 |         2.22 |
| SD-reject   | VILL      |       -5.01 |        2.43 |       -3.81 |         2.33 |         1.07 |         2.59 |
| SD-reject   | YEBE      |        0.38 |       -1.86 |       -2.02 |         2.38 |         1.09 |         2.66 |

## WLS residuals and local tests (ENU of from-station)

|   baseline_id | from   | to   |   length_m |   v_e_mm |   v_n_mm |   v_u_mm |   T_b |   T_b_critical |   local_redundancy |   sd_statistic |   sd_critical |   sd_dir_elevation_deg |
|--------------:|:-------|:-----|-----------:|---------:|---------:|---------:|------:|---------------:|-------------------:|---------------:|--------------:|-----------------------:|
|             1 | MADR   | CEBR |   1.04e+04 |   -1.6   |   0.867  |     8.82 |  9.34 |           16.3 |               2.11 |           3.06 |          4.03 |                   78.3 |
|             2 | MADR   | VILL |   2.53e+04 |   -2.07  |  -0.0979 |     6.65 |  4.91 |           16.3 |               2.11 |           2.22 |          4.03 |                   72.3 |
|             3 | VILL   | YEBE |   7.38e+04 |   -2.82  |   3.1    |    -5.44 | 10.7  |           16.3 |               1.64 |           3.27 |          4.03 |                   52   |
|             4 | CEBR   | VILL |   3.53e+04 |    0.486 |  -1.38   |    16.5  | 27.5  |           16.3 |               2.14 |           5.24 |          4.03 |                   84.4 |
|             5 | MADR   | CEBR |   1.04e+04 |    3.66  |  -0.107  |    -5.86 |  4.22 |           16.3 |               2.17 |           2.05 |          4.03 |                   58.7 |
|             6 | MADR   | YEBE |   9.9e+04  |    2.9   |  -1.78   |    -1.55 |  3.95 |           16.3 |               1.78 |           1.99 |          4.03 |                   22   |
|             7 | CEBR   | YEBE |   1.09e+05 |    0.612 |  -1.51   |     6.84 |  6.09 |           16.3 |               1.73 |           2.47 |          4.03 |                   77.2 |
|             8 | CEBR   | VILL |   3.53e+04 |   -1.01  |   3.62   |   -19.3  | 39.8  |           16.3 |               2.13 |           6.31 |          4.03 |                   79.3 |
|             9 | MADR   | VILL |   2.53e+04 |   -2.62  |   1.25   |    -9.75 |  9.25 |           16.3 |               2.17 |           3.04 |          4.03 |                   74.1 |
