"""Comprobaciones de robustez de los resultados preliminares.

(a) Repetición con 50 semillas: ¿la discrepancia queda cerca del oráculo de forma
    sistemática o solo con la semilla 2026?
(b) Sensibilidad a la cota de ruido: delta = c * ||e|| con c en {0.5, 1, 2}.
(c) Línea base sin ruido (e = 0): mejor error alcanzable en la malla; mide el sesgo
    que introduce la regularización aunque los datos sean exactos.
"""

import numpy as np
from comun import ETAS, KS, LAMBDA_REL, N, S, TAU, guardar_tabla

from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import seleccion as SEL
from reconstruccion.metricas import error_relativo

t = datos.malla(N)
A = datos.matriz_suavizado(t, S)
x_true = datos.senal_verdadera(t)
b_limpio = A @ x_true
svd = M.SVD.de(A)
lambdas = svd.sigma[0] ** 2 * LAMBDA_REL
SEMILLAS = range(1, 51)


def errores_tikhonov(b):
    return np.array([error_relativo(M.tikhonov_espectral(svd, b, lam), x_true) for lam in lambdas])


def errores_tsvd(b):
    return np.array([error_relativo(M.svd_truncada(svd, b, int(k)), x_true) for k in KS])


# (a) y (b)
print("(a)-(b) 50 semillas, tres niveles de ruido, tres cotas de ruido")
filas = []
for eta in ETAS:
    cocientes = {c: {"tik": [], "tsvd": []} for c in (0.5, 1.0, 2.0)}
    err_oraculo_tik = []
    for semilla in SEMILLAS:
        e = datos.ruido(b_limpio, eta, datos.direccion_ruido(N, semilla))
        b = b_limpio + e
        err_l = errores_tikhonov(b)
        err_k = errores_tsvd(b)
        err_oraculo_tik.append(err_l.min())
        for c in cocientes:
            delta = c * np.linalg.norm(e)
            _, i_d = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta, TAU)
            _, j_d = SEL.discrepancia_tsvd(A, b, svd, KS, delta, TAU)
            cocientes[c]["tik"].append(err_l[i_d] / err_l.min() if i_d is not None else np.nan)
            cocientes[c]["tsvd"].append(err_k[j_d] / err_k.min() if j_d is not None else np.nan)
    def mediana(v):
        return float(np.nanmedian(v)) if np.isfinite(v).any() else None

    def maximo(v):
        return float(np.nanmax(v)) if np.isfinite(v).any() else None

    for c, d in cocientes.items():
        tik, tsvd = np.array(d["tik"]), np.array(d["tsvd"])
        filas.append({
            "eta": eta, "delta / ||e||": c,
            "Tikhonov: mediana error_disc/error_orac": mediana(tik),
            "Tikhonov: máximo": maximo(tik),
            "Tikhonov: semillas no viables": int(np.isnan(tik).sum()),
            "TSVD: mediana error_disc/error_orac": mediana(tsvd),
            "TSVD: máximo": maximo(tsvd),
            "TSVD: semillas no viables": int(np.isnan(tsvd).sum()),
        })
    print(f"  eta = {eta:g}: error oráculo Tikhonov, mediana en 50 semillas = {np.median(err_oraculo_tik):.4f}")
guardar_tabla(filas, "03_robustez_semillas")

# (c) Línea base sin ruido
err_l0 = errores_tikhonov(b_limpio)
err_k0 = errores_tsvd(b_limpio)
i0, j0 = int(np.argmin(err_l0)), int(np.argmin(err_k0))
x_pinv0, rango0, _ = M.pseudoinversa_numerica(A, b_limpio)
filas0 = [
    {"método (e = 0)": "Tikhonov, mejor λ de la malla", "parámetro": f"λ/σ1² = {LAMBDA_REL[i0]:.1e}",
     "error relativo": float(err_l0[i0])},
    {"método (e = 0)": "SVD truncada, mejor k", "parámetro": f"k = {KS[j0]}", "error relativo": float(err_k0[j0])},
    {"método (e = 0)": "Pseudoinversa numérica (scipy)", "parámetro": f"rango {rango0}",
     "error relativo": error_relativo(x_pinv0, x_true)},
]
print("(c) sin ruido:")
for f in filas0:
    print(f"  {f['método (e = 0)']:32s} {f['parámetro']:18s} error = {f['error relativo']:.3g}")
guardar_tabla(filas0, "03_linea_base_sin_ruido")
