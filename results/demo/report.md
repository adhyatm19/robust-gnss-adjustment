# Demo experiment report

> Synthetic baseline-level research prototype; not raw GNSS/RINEX processing.

- seed: `42`
- stations: 8, baselines: 16
- true outlier baselines: none

## Method comparison

| method      |    rmse_m |   mean_pos_error_m |   max_pos_error_m |   tp |   fp |   fn |   precision |   recall |   f1 |   fpr | localization_correct   |   iterations | converged   |   runtime_s |   condition_number |   variance_factor |   n_flagged |   n_true_outliers |
|:------------|----------:|-------------------:|------------------:|-----:|-----:|-----:|------------:|---------:|-----:|------:|:-----------------------|-------------:|:------------|------------:|-------------------:|------------------:|------------:|------------------:|
| WLS         | 0.0040851 |          0.0066533 |           0.01166 |    0 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            1 | True        |   0.0030222 |             94.218 |           0.75722 |           0 |                 0 |
| DIA-reject  | 0.0040851 |          0.0066533 |           0.01166 |    0 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            1 | True        |   0.0030273 |             94.218 |           0.75722 |           0 |                 0 |
| DIA-inflate | 0.0040851 |          0.0066533 |           0.01166 |    0 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            1 | True        |   0.0020733 |             94.218 |           0.75722 |           0 |                 0 |
| Huber       | 0.0040851 |          0.0066533 |           0.01166 |    0 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            2 | True        |   0.0019605 |             94.218 |           0.75722 |           0 |                 0 |
| Hampel      | 0.0040851 |          0.0066533 |           0.01166 |    0 |    0 |    0 |           1 |        1 |    1 |     0 | True                   |            2 | True        |   0.0019244 |             94.218 |           0.75722 |           0 |                 0 |

## DIA audit trails

### DIA-reject (status: clean)

| iter | T_global | crit | passed | identified | T_b | action |
|---|---|---|---|---|---|---|
| 0 | 20.44 | 40.11 | True | None |  | stopped:global_test_passed |

### DIA-inflate (status: clean)

| iter | T_global | crit | passed | identified | T_b | action |
|---|---|---|---|---|---|---|
| 0 | 20.44 | 40.11 | True | None |  | stopped:global_test_passed |

### Huber: converged=True in 2 iterations; down-weighted baselines: none

### Hampel: converged=True in 2 iterations; down-weighted baselines: none
