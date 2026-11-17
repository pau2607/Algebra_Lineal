| eta | fracción error² en [0.55,0.8] | ‖sesgo‖/‖x‖ | ‖ruido prop.‖/‖x‖ | res. pinv explícita/‖e‖ | res. TSVD k=rango/‖e‖ | ‖pinv‖ | ε‖A‖‖x_completa‖/‖e‖ | res. mínimo Tikhonov en malla/‖e‖ | σ_i > ‖e‖/√n | |u_iᵀAx| > ‖e‖/√n | último i informado | primer k con residuo TSVD creciente |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0001 | 0.983 | 0.09402 | 0.004603 | 0.6884 | 0.5344 | 4.014e+08 | 1.481 | 0.6957 | 50 | 40 | 41 | 94 |
| 0.001 | 0.9865 | 0.1079 | 0.007195 | 0.5359 | 0.5344 | 4.014e+09 | 1.481 | 0.6956 | 44 | 33 | 37 | 94 |
| 0.01 | 0.9612 | 0.1168 | 0.01892 | 0.5344 | 0.5344 | 4.014e+10 | 1.481 | 0.6956 | 37 | 23 | 26 | 94 |
