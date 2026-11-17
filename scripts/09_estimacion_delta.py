"""Estimación del nivel de ruido delta a partir de los datos (desarrollo §8.4).

Hasta ahora la discrepancia usaba delta = ||e|| exacto, que solo conocemos porque
generamos el ruido. Aquí delta se estima con la cola de los coeficientes de Picard
(`seleccion.estimar_delta`) y se mide:
  * qué tan cerca queda delta_hat de ||e|| (50 semillas, tres eta);
  * cuánto empeora la discrepancia por usar delta_hat en lugar de ||e||.
"""

import numpy as np
from comun import (C1, C2, C4, COLOR_ETA, ETAS, KS, LAMBDA_REL, N, S, SEMILLA, TAU, TENUE, TINTA, etiqueta_eta,
                   guardar, guardar_tabla, plt)

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
FRACCION = 0.5
I0 = int(round((1 - FRACCION) * N))

print("(1) Calidad del estimador y efecto en la discrepancia, 50 semillas")
filas = []
for eta in ETAS:
    cocientes, err_exacto, err_estimado, err_oraculo, err_tsvd_est, err_tsvd_orac, no_viables = [], [], [], [], [], [], 0
    for semilla in range(1, 51):
        e = datos.ruido(b_limpio, eta, datos.direccion_ruido(N, semilla))
        b = b_limpio + e
        delta_hat = SEL.estimar_delta(svd, b, FRACCION)
        cocientes.append(delta_hat / np.linalg.norm(e))
        err_l = np.array([error_relativo(M.tikhonov_espectral(svd, b, lam), x_true) for lam in lambdas])
        err_k = np.array([error_relativo(M.svd_truncada(svd, b, int(k)), x_true) for k in KS])
        _, i_ex = SEL.discrepancia_tikhonov(A, b, svd, lambdas, np.linalg.norm(e), TAU)
        _, i_es = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta_hat, TAU)
        _, j_es = SEL.discrepancia_tsvd(A, b, svd, KS, delta_hat, TAU)
        err_oraculo.append(err_l.min())
        err_exacto.append(err_l[i_ex])
        err_tsvd_orac.append(err_k.min())
        if i_es is None or j_es is None:
            no_viables += 1
            err_estimado.append(np.nan)
            err_tsvd_est.append(np.nan)
        else:
            err_estimado.append(err_l[i_es])
            err_tsvd_est.append(err_k[j_es])
    c = np.array(cocientes)
    filas.append({
        "eta": eta, "δ̂/‖e‖ mediana": float(np.median(c)), "δ̂/‖e‖ mín": float(c.min()), "δ̂/‖e‖ máx": float(c.max()),
        "Tikhonov: mediana error, δ exacto": float(np.median(err_exacto)),
        "Tikhonov: mediana error, δ̂ estimado": float(np.nanmedian(err_estimado)),
        "Tikhonov: mediana error oráculo": float(np.median(err_oraculo)),
        "TSVD: mediana error, δ̂ estimado": float(np.nanmedian(err_tsvd_est)),
        "TSVD: mediana error oráculo": float(np.median(err_tsvd_orac)),
        "semillas no viables con δ̂": no_viables,
    })
    print(f"  eta = {eta:g}: δ̂/‖e‖ en [{c.min():.2f}, {c.max():.2f}] (mediana {np.median(c):.3f}); "
          f"error Tikhonov: exacto {np.median(err_exacto):.4f}, estimado {np.nanmedian(err_estimado):.4f}, "
          f"oráculo {np.median(err_oraculo):.4f}; no viables {no_viables}")
guardar_tabla(filas, "09_estimacion_delta")

print("(2) Refinar lambda por bisección entre dos puntos de la malla (semilla 2026)")
filas = []
for eta in ETAS:
    e = datos.ruido(b_limpio, eta, datos.direccion_ruido(N, SEMILLA))
    b = b_limpio + e
    delta_hat = SEL.estimar_delta(svd, b, FRACCION)
    lam_m, _ = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta_hat, TAU)
    lam_r = SEL.discrepancia_tikhonov_refinada(A, b, svd, lambdas, delta_hat, TAU)
    filas.append({"eta": eta, "λ/σ1² malla": lam_m / svd.sigma[0] ** 2, "λ/σ1² refinado": lam_r / svd.sigma[0] ** 2,
                  "error malla": error_relativo(M.tikhonov_espectral(svd, b, lam_m), x_true),
                  "error refinado": error_relativo(M.tikhonov_espectral(svd, b, lam_r), x_true)})
    print(f"  eta = {eta:g}: λ/σ1² {lam_m / svd.sigma[0] ** 2:.3e} -> {lam_r / svd.sigma[0] ** 2:.3e}; "
          f"error {filas[-1]['error malla']:.4f} -> {filas[-1]['error refinado']:.4f}")
guardar_tabla(filas, "09_lambda_refinado")

# Figura 18: coeficientes de Picard con la cola usada y el nivel estimado (semilla 2026)
z = datos.direccion_ruido(N, SEMILLA)
fig, ax = plt.subplots(figsize=(7.5, 4.0))
i = np.arange(1, N + 1)
ax.semilogy(i, svd.sigma, color=TINTA, lw=2, label=r"$\sigma_i$")
ax.semilogy(i, np.abs(svd.coeficientes(b_limpio)) + 1e-300, "-", lw=1, color=TENUE, label=r"$|u_i^TAx_{true}|$ (sin ruido)")
for eta in ETAS:
    e = datos.ruido(b_limpio, eta, z)
    b = b_limpio + e
    delta_hat = SEL.estimar_delta(svd, b, FRACCION)
    ax.semilogy(i, np.abs(svd.coeficientes(b)), "o", ms=3, color=COLOR_ETA[eta], label=rf"$|u_i^Tb|$, {etiqueta_eta(eta)}")
    ax.axhline(delta_hat / np.sqrt(N), color=COLOR_ETA[eta], lw=1, ls="--")
    ax.axhline(np.linalg.norm(e) / np.sqrt(N), color=COLOR_ETA[eta], lw=1, ls=":")
ax.axvspan(I0 + 0.5, N + 0.5, color=C4, alpha=0.12, lw=0)
ax.text(I0 + 2, 1e-18, f"cola usada para estimar δ (i > {I0})", fontsize=8, color=TINTA)
ax.plot([], [], color=TENUE, ls="--", label=r"$\hat\delta/\sqrt{n}$ (estimado)")
ax.plot([], [], color=TENUE, ls=":", label=r"$\|e\|/\sqrt{n}$ (real)")
ax.set_ylim(1e-19, 10)
ax.set_xlabel("índice i")
ax.set_title("Estimar el nivel de ruido con la cola de los coeficientes de Picard")
ax.legend(loc="lower left", fontsize=7, ncol=2)
guardar(fig, "18_estimacion_delta.png")
