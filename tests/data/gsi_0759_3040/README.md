# GSI GEONET 0759 / 3040 test observations (real GNSS data)

Real RINEX 2.10 observation and GPS navigation files for two stations of
the Geospatial Information Authority of Japan (GSI) GEONET network,
redistributed unchanged from the RTKLIB source repository
(`https://github.com/tomojitakasu/RTKLIB`, tag `v2.4.3-b34`, commit
`180043ee24b6d2b168f98b64be15f69d50046b1a`, directory `test/data/rinex/`),
which is distributed under the BSD 2-clause licence.

| item | value |
|---|---|
| stations | 0759, 3040 (GSI GEONET) |
| date | 2005-04-02, 00:00:00 - 00:59:30 GPST (1 hour) |
| sampling | 30 s, 120 epochs per station |
| constellation / signals | GPS only; L1, C1, L2, P2 |
| receivers / antennas | TRIMBLE 5700 / TRM29659.00 (both stations) |
| antenna delta H/E/N | 0 / 0 / 0 (header) |
| approx. baseline length | 3.34 km |

Checksums are in `SHA256SUMS`. The header `APPROX POSITION XYZ` values are
receiver-approximate and differ from the processed baseline by about 0.17 m;
**no independent reference coordinates are available for this pair**, so it
is used only to verify the processing chain (format, conventions,
covariance reconstruction, antenna-setting precedence), not accuracy.
