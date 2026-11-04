| eta | método | parámetro | error relativo | residuo | residuo / delta | norma |
|---|---|---|---|---|---|---|
| 0.0001 | Pseudoinversa completa (128 σ) | k = 128, sin umbral | 6.287e+11 | 0.000456 | 0.9733 | 3.121e+12 |
| 0.0001 | Pseudoinversa numérica (scipy) | rango 89, umbral 2.84e-14 | 8.088e+07 | 0.0003225 | 0.6884 | 4.014e+08 |
| 0.0001 | SVD truncada · discrepancia | k = 38 | 0.09393 | 0.0003995 | 0.8527 | 4.942 |
| 0.0001 | SVD truncada · oráculo | k = 40 | 0.09307 | 0.0003702 | 0.7903 | 4.942 |
| 0.0001 | Tikhonov · discrepancia | λ/σ1² = 6.81e-06 | 0.09406 | 0.0004792 | 1.023 | 4.939 |
| 0.0001 | Tikhonov · oráculo | λ/σ1² = 3.16e-07 | 0.09137 | 0.0003608 | 0.7701 | 4.943 |
| 0.001 | Pseudoinversa completa (128 σ) | k = 128, sin umbral | 6.287e+12 | 0.004695 | 1.002 | 3.121e+13 |
| 0.001 | Pseudoinversa numérica (scipy) | rango 89, umbral 2.84e-14 | 8.088e+08 | 0.00251 | 0.5359 | 4.014e+09 |
| 0.001 | SVD truncada · discrepancia | k = 31 | 0.1139 | 0.004737 | 1.011 | 4.932 |
| 0.001 | SVD truncada · oráculo | k = 38 | 0.1029 | 0.003725 | 0.7952 | 4.946 |
| 0.001 | Tikhonov · discrepancia | λ/σ1² = 3.16e-04 | 0.1087 | 0.004834 | 1.032 | 4.929 |
| 0.001 | Tikhonov · oráculo | λ/σ1² = 2.15e-05 | 0.1008 | 0.00378 | 0.8068 | 4.938 |
| 0.01 | Pseudoinversa completa (128 σ) | k = 128, sin umbral | 6.287e+13 | 0.04753 | 1.015 | 3.121e+14 |
| 0.01 | Pseudoinversa numérica (scipy) | rango 89, umbral 2.84e-14 | 8.088e+09 | 0.02504 | 0.5344 | 4.014e+10 |
| 0.01 | SVD truncada · discrepancia | k = 23 | 0.1218 | 0.04193 | 0.8951 | 4.936 |
| 0.01 | SVD truncada · oráculo | k = 27 | 0.1185 | 0.04013 | 0.8567 | 4.938 |
| 0.01 | Tikhonov · discrepancia | λ/σ1² = 4.64e-03 | 0.1194 | 0.04844 | 1.034 | 4.904 |
| 0.01 | Tikhonov · oráculo | λ/σ1² = 3.16e-03 | 0.1191 | 0.04409 | 0.9411 | 4.914 |
