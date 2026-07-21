# Demo experiment report

> Synthetic baseline-level research prototype; not raw GNSS/RINEX processing.

- seed: `42`
- stations: 8, baselines: 16
- true outlier baselines: [3]

## Method comparison

| method      |    rmse_m |   mean_pos_error_m |   max_pos_error_m |   tp |   fp |   fn |   precision |   recall |   f1 |   fpr | localization_correct   |   iterations | converged   |   runtime_s |   condition_number |   variance_factor |   n_flagged |   n_true_outliers |
|:------------|----------:|-------------------:|------------------:|-----:|-----:|-----:|------------:|---------:|-----:|------:|:-----------------------|-------------:|:------------|------------:|-------------------:|------------------:|------------:|------------------:|
| WLS         | 0.0087649 |          0.013436  |          0.022697 |    0 |    0 |    1 |           0 |        0 |    0 |     0 | False                  |            1 | True        |   0.0010662 |             94.218 |           2.8009  |           0 |                 1 |
| DIA-reject  | 0.0044312 |          0.0071728 |          0.011212 |    1 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            2 | True        |   0.0036312 |             89.171 |           0.78514 |           1 |                 1 |
| DIA-inflate | 0.0044213 |          0.0071723 |          0.011084 |    1 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            2 | True        |   0.0037377 |             89.208 |           0.73442 |           1 |                 1 |
| Huber       | 0.0054676 |          0.008689  |          0.015519 |    1 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            9 | True        |   0.0059878 |             90.536 |           1.6563  |           1 |                 1 |
| Hampel      | 0.0048554 |          0.0077992 |          0.013989 |    1 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |           12 | True        |   0.0076338 |             89.98  |           1.3378  |           1 |                 1 |

## DIA audit trails

### DIA-reject (status: adapted)

| iter | T_global | crit | passed | identified | T_b | action |
|---|---|---|---|---|---|---|
| 0 | 75.63 | 40.11 | False | 3 | 56.78 | reject:baseline_3 |
| 1 | 18.84 | 36.42 | True | None |  | stopped:global_test_passed |

### DIA-inflate (status: adapted)

| iter | T_global | crit | passed | identified | T_b | action |
|---|---|---|---|---|---|---|
| 0 | 75.63 | 40.11 | False | 3 | 56.78 | inflate:baseline_3:x100 |
| 1 | 19.83 | 40.11 | True | None |  | stopped:global_test_passed |

### Huber: converged=True in 9 iterations; down-weighted baselines: [3]

### Hampel: converged=True in 12 iterations; down-weighted baselines: [3]
