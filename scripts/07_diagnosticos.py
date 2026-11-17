"""Diagnósticos que respaldan las interpretaciones de docs/resultados_preliminares.md.

Para la semilla 2026 y cada eta:
  * dónde se concentra el error de Tikhonov (región del pulso);
  * separación del error en sesgo y ruido propagado;
  * residuos de la pseudoinversa (explícita vs. espectral) y residuo mínimo de la malla;
  * cuántas direcciones informan los datos (coeficientes de Picard sobre el nivel del ruido);
  * si el residuo calculado de la TSVD es monótono en k.
Además:
  * barrido de eta hasta 1e-10 (¿el error baja con menos ruido?);
  * pérdida de rango al formar A^T A y contraste con las ecuaciones normales regularizadas.
"""

import numpy as np
from comun import ETAS, KS, LAMBDA_REL, N, S, SEMILLA, guardar_tabla

from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import seleccion as SEL
from reconstruccion.metricas import error_relativo, residuo

EPS = np.finfo(float).eps
t = datos.malla(N)
A = datos.matriz_suavizado(t, S)
x_true = datos.senal_verdadera(t)
b_limpio = A @ x_true
z = datos.direccion_ruido(N, SEMILLA)
svd = M.SVD.de(A)
sigma = svd.sigma
lambdas = sigma[0] ** 2 * LAMBDA_REL
pulso = (t >= 0.55) & (t <= 0.8)
c_limpio = np.abs(svd.coeficientes(b_limpio))

# --------------------------------------------------------------------------
# (1) Diagnósticos por nivel de ruido
# --------------------------------------------------------------------------
print("(1) Diagnósticos por nivel de ruido (semilla 2026)")
filas = []
for eta in ETAS:
    e = datos.ruido(b_limpio, eta, z)
    b = b_limpio + e
    ne = np.linalg.norm(e)
    lam, _ = SEL.discrepancia_tikhonov(A, b, svd, lambdas, ne)
    x_l = M.tikhonov_espectral(svd, b, lam)
    r = x_l - x_true
    sesgo = M.tikhonov_espectral(svd, b_limpio, lam) - x_true
    ruido_prop = M.tikhonov_espectral(svd, e, lam)

    x_pinv, rango, _ = M.pseudoinversa_numerica(A, b)
    x_k89 = M.svd_truncada(svd, b, rango)
    x_completa = M.pseudoinversa_completa(svd, b)
    res_k = np.array([residuo(A, M.svd_truncada(svd, b, int(k)), b) for k in KS])
    subidas = np.flatnonzero(np.diff(res_k) > 0)
    nivel = ne / np.sqrt(N)
    informados = np.flatnonzero(c_limpio > nivel)

    fila = {
        "eta": eta,
        "fracción error² en [0.55,0.8]": float((r[pulso] ** 2).sum() / (r ** 2).sum()),
        "‖sesgo‖/‖x‖": float(np.linalg.norm(sesgo) / np.linalg.norm(x_true)),
        "‖ruido prop.‖/‖x‖": float(np.linalg.norm(ruido_prop) / np.linalg.norm(x_true)),
        "res. pinv explícita/‖e‖": residuo(A, x_pinv, b) / ne,
        "res. TSVD k=rango/‖e‖": residuo(A, x_k89, b) / ne,
        "‖pinv‖": float(np.linalg.norm(x_pinv)),
        "ε‖A‖‖x_completa‖/‖e‖": float(EPS * sigma[0] * np.linalg.norm(x_completa) / ne),
        "res. mínimo Tikhonov en malla/‖e‖": residuo(A, M.tikhonov_espectral(svd, b, lambdas[0]), b) / ne,
        "σ_i > ‖e‖/√n": int((sigma > nivel).sum()),
        "|u_iᵀAx| > ‖e‖/√n": int(informados.size),
        "último i informado": int(informados[-1] + 1),
        "primer k con residuo TSVD creciente": int(KS[subidas[0] + 1]) if subidas.size else None,
    }
    filas.append(fila)
    print("  " + ", ".join(f"{k} = {v:.3g}" if isinstance(v, float) else f"{k} = {v}" for k, v in fila.items()))
guardar_tabla(filas, "07_diagnosticos")

# --------------------------------------------------------------------------
# (2) Barrido de eta: ¿el error del oráculo baja con menos ruido?
# --------------------------------------------------------------------------
print("(2) Barrido de eta")
LAMBDA_REL_EXT = np.logspace(-20, 0, 121)   # malla ampliada solo para este barrido
lambdas_ext = sigma[0] ** 2 * LAMBDA_REL_EXT
filas = []
for eta in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 0.0]:
    b = b_limpio + datos.ruido(b_limpio, eta, z)
    err_l = [error_relativo(M.tikhonov_espectral(svd, b, lam), x_true) for lam in lambdas]
    err_e = [error_relativo(M.tikhonov_espectral(svd, b, lam), x_true) for lam in lambdas_ext]
    err_k = [error_relativo(M.svd_truncada(svd, b, int(k)), x_true) for k in KS]
    i, j, k = int(np.argmin(err_l)), int(np.argmin(err_e)), int(np.argmin(err_k))
    filas.append({
        "eta": eta,
        "Tikhonov oráculo (malla 1e-10)": float(err_l[i]), "λ/σ1² (malla 1e-10)": float(LAMBDA_REL[i]),
        "Tikhonov oráculo (malla 1e-20)": float(err_e[j]), "λ/σ1² (malla 1e-20)": float(LAMBDA_REL_EXT[j]),
        "TSVD oráculo": float(err_k[k]), "k oráculo": int(KS[k]),
    })
    print(f"  eta = {eta:g}: Tikhonov {err_l[i]:.4f} / {err_e[j]:.4f} (malla ampliada), TSVD {err_k[k]:.4f} (k = {KS[k]})")
guardar_tabla(filas, "07_barrido_eta")

# --------------------------------------------------------------------------
# (3) Ecuaciones normales
# --------------------------------------------------------------------------
print("(3) Ecuaciones normales")
b = b_limpio + datos.ruido(b_limpio, 1e-3, z)
dif = max(np.linalg.norm(M.tikhonov_normales(A, b, lam) - M.tikhonov_espectral(svd, b, lam))
          / np.linalg.norm(M.tikhonov_espectral(svd, b, lam)) for lam in lambdas)
filas = [
    {"comprobación": "rango numérico de A (matrix_rank)", "valor": int(np.linalg.matrix_rank(A))},
    {"comprobación": "rango numérico de AᵀA (matrix_rank)", "valor": int(np.linalg.matrix_rank(A.T @ A))},
    {"comprobación": "σ_i > √ε·σ1 (sobreviven al formar AᵀA)", "valor": int((sigma > np.sqrt(EPS) * sigma[0]).sum())},
    {"comprobación": "máx. dif. relativa normales regularizadas vs espectral (η=1e-3, 61 λ)", "valor": float(dif)},
]
for f in filas:
    print(f"  {f['comprobación']}: {f['valor']:.3g}" if isinstance(f["valor"], float) else f"  {f['comprobación']}: {f['valor']}")
guardar_tabla(filas, "07_ecuaciones_normales")
