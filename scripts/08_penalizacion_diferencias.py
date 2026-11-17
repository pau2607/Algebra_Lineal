"""Ampliación: Tikhonov generalizado con penalización de diferencias (desarrollo §10).

Compara tres penalizaciones lambda ||L x||^2:
    L = I   (Tikhonov clásico),  L_1 primeras diferencias,  L_2 segundas diferencias,
sobre dos señales:
    * la señal del experimento (gaussiana + pulso): ¿mejoran los bordes?
    * una señal suave (seno + rampa, no nula en los extremos): ¿importa la penalización?
Para cada caso: lambda por discrepancia (sin x_true) y por oráculo, en los tres eta y
en 50 semillas. También comprueba la GSVD contra el sistema aumentado con QR.
"""

import numpy as np
from comun import (C1, C2, C3, COLOR_ETA, ETAS, N, S, SEMILLA, TAU, TENUE, TINTA, etiqueta_eta, guardar,
                   guardar_tabla, plt)

from reconstruccion import datos
from reconstruccion import general as G
from reconstruccion.metricas import error_relativo, residuo

t = datos.malla(N)
A = datos.matriz_suavizado(t, S)
z = datos.direccion_ruido(N, SEMILLA)
LAMBDAS = np.logspace(-14, 4, 109)           # malla absoluta (sigma_1(A) ≈ 1), 6 puntos por década
PENALIZACIONES = {"I": np.eye(N), "L1": G.primeras_diferencias(N), "L2": G.segundas_diferencias(N)}
COLOR_L = {"I": C2, "L1": C1, "L2": C3}
SENALES = {"pulso": datos.senal_verdadera(t), "suave": datos.senal_suave(t)}
PULSO = (t >= 0.55) & (t <= 0.8)

# --------------------------------------------------------------------------
# (1) Unicidad y GSVD
# --------------------------------------------------------------------------
print("(1) Unicidad y GSVD")
gsvds = {}
filas = []
for nombre, L in PENALIZACIONES.items():
    gsvds[nombre] = G.GSVD.de(A, L)
    g = gsvds[nombre]
    rec_A = np.linalg.norm(g.U @ np.diag(g.c) @ np.linalg.inv(g.X) - A) / np.linalg.norm(A)
    b = A @ SENALES["pulso"] + datos.ruido(A @ SENALES["pulso"], 1e-3, z)
    dif = max(np.linalg.norm(G.tikhonov_general_espectral(g, b, lam) - G.tikhonov_general_aumentado(A, L, b, lam))
              / np.linalg.norm(G.tikhonov_general_aumentado(A, L, b, lam)) for lam in LAMBDAS[::8])
    gam = g.gamma[np.isfinite(g.gamma)]
    filas.append({"L": nombre, "dim ker L": N - np.linalg.matrix_rank(L), "ker A ∩ ker L = {0}": G.nucleos_disjuntos(A, L),
                  "‖UCX⁻¹−A‖/‖A‖": float(rec_A), "máx dif espectral vs QR": float(dif),
                  "γ máx (finito)": float(gam.max()), "γ mín": float(gam.min())})
    print(f"  {nombre}: dim ker L = {filas[-1]['dim ker L']}, núcleos disjuntos = {G.nucleos_disjuntos(A, L)}, "
          f"GSVD reconstruye A con error {rec_A:.1e}, espectral vs QR {dif:.1e}")
guardar_tabla(filas, "08_gsvd")

# --------------------------------------------------------------------------
# (2) Comparación en la semilla 2026
# --------------------------------------------------------------------------
print("(2) Comparación, semilla 2026")
filas = []
rec = {}          # (senal, eta, L) -> reconstrucción por discrepancia
curvas = {}       # (senal, eta, L) -> errores en la malla
for senal, x_true in SENALES.items():
    b_limpio = A @ x_true
    for eta in ETAS:
        e = datos.ruido(b_limpio, eta, z)
        b = b_limpio + e
        delta = np.linalg.norm(e)
        for nombre, L in PENALIZACIONES.items():
            g = gsvds[nombre]
            X = [G.tikhonov_general_espectral(g, b, lam) for lam in LAMBDAS]
            err = np.array([error_relativo(xx, x_true) for xx in X])
            lam_d, i_d = G.discrepancia_general(A, L, b, LAMBDAS, delta, TAU, gsvd=g)
            i_o = int(np.argmin(err))
            rec[(senal, eta, nombre)] = X[i_d]
            curvas[(senal, eta, nombre)] = err
            r = X[i_d] - x_true
            filas.append({"señal": senal, "eta": eta, "L": nombre,
                          "λ discrepancia": float(lam_d), "error discrepancia": float(err[i_d]),
                          "λ oráculo": float(LAMBDAS[i_o]), "error oráculo": float(err[i_o]),
                          "error en [0.55,0.8] (disc.)": float(np.linalg.norm(r[PULSO]) / np.linalg.norm(x_true)),
                          "residuo/δ (disc.)": residuo(A, X[i_d], b) / delta})
            print(f"  {senal:6s} eta={eta:g} {nombre:2s}: discrepancia {err[i_d]:.4f} (λ={lam_d:.1e}) | "
                  f"oráculo {err[i_o]:.4f} (λ={LAMBDAS[i_o]:.1e})")
guardar_tabla(filas, "08_penalizaciones")

# --------------------------------------------------------------------------
# (3) Robustez: 50 semillas, cociente error(L) / error(I) con discrepancia y con oráculo
# --------------------------------------------------------------------------
print("(3) Robustez con 50 semillas")
filas = []
for senal, x_true in SENALES.items():
    b_limpio = A @ x_true
    for eta in ETAS:
        err_d = {k: [] for k in PENALIZACIONES}
        err_o = {k: [] for k in PENALIZACIONES}
        for semilla in range(1, 51):
            e = datos.ruido(b_limpio, eta, datos.direccion_ruido(N, semilla))
            b = b_limpio + e
            for nombre, L in PENALIZACIONES.items():
                g = gsvds[nombre]
                err = np.array([error_relativo(G.tikhonov_general_espectral(g, b, lam), x_true) for lam in LAMBDAS])
                _, i_d = G.discrepancia_general(A, L, b, LAMBDAS, np.linalg.norm(e), TAU, gsvd=g)
                err_d[nombre].append(err[i_d] if i_d is not None else np.nan)
                err_o[nombre].append(err.min())
        fila = {"señal": senal, "eta": eta}
        for nombre in PENALIZACIONES:
            fila[f"{nombre}: mediana error disc."] = float(np.nanmedian(err_d[nombre]))
            fila[f"{nombre}: mediana error orác."] = float(np.median(err_o[nombre]))
        filas.append(fila)
        print(f"  {senal:6s} eta={eta:g}: " + ", ".join(f"{k} {fila[f'{k}: mediana error disc.']:.4f}" for k in PENALIZACIONES))
guardar_tabla(filas, "08_robustez_penalizaciones")

# --------------------------------------------------------------------------
# Figuras
# --------------------------------------------------------------------------
eta = 1e-3

# 15. Reconstrucciones con las tres penalizaciones: pulso (y su ampliación) y señal suave (y su error)
fig, axs = plt.subplots(2, 2, figsize=(12, 7.2))
axs = axs.ravel()
x_true = SENALES["pulso"]
axs[0].plot(t, x_true, color=TINTA, lw=2.4, label="verdad")
for nombre in PENALIZACIONES:
    axs[0].plot(t, rec[("pulso", eta, nombre)], color=COLOR_L[nombre], lw=1.5, label=f"L = {nombre}")
axs[0].set_title(rf"Señal con pulso, {etiqueta_eta(eta)}, λ por discrepancia")
axs[0].legend(loc="upper right", fontsize=8)
axs[1].plot(t, x_true, color=TINTA, lw=2.4, label="verdad")
for nombre in PENALIZACIONES:
    axs[1].plot(t, rec[("pulso", eta, nombre)], color=COLOR_L[nombre], lw=1.5, marker="o", ms=3, label=f"L = {nombre}")
axs[1].set_xlim(0.5, 0.85)
axs[1].set_title("Ampliación del pulso: las tres penalizaciones casi coinciden")
x_true = SENALES["suave"]
axs[2].plot(t, x_true, color=TINTA, lw=2.4, label="verdad")
for nombre in PENALIZACIONES:
    axs[2].plot(t, rec[("suave", eta, nombre)], color=COLOR_L[nombre], lw=1.5, label=f"L = {nombre}")
axs[2].set_title(rf"Señal suave, {etiqueta_eta(eta)}, λ por discrepancia")
axs[2].legend(loc="upper right", fontsize=8)
for nombre in PENALIZACIONES:
    axs[3].semilogy(t, np.abs(rec[("suave", eta, nombre)] - x_true) + 1e-12, color=COLOR_L[nombre], lw=1.3,
                    label=f"L = {nombre}: error rel. {error_relativo(rec[('suave', eta, nombre)], x_true):.1e}")
axs[3].set_ylim(1e-6, 1)
axs[3].set_title("Señal suave: error punto a punto $|x_\\lambda - x_{true}|$")
axs[3].legend(loc="lower right", fontsize=8)
for ax in axs:
    ax.set_xlabel("t")
fig.tight_layout()
guardar(fig, "15_penalizaciones_reconstrucciones.png")

# 16. Error frente a lambda para las tres penalizaciones, ambas señales
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for ax, senal in zip(axs, SENALES):
    for nombre in PENALIZACIONES:
        ax.loglog(LAMBDAS, curvas[(senal, eta, nombre)], color=COLOR_L[nombre], label=f"L = {nombre}")
    ax.set_title(rf"Error relativo frente a $\lambda$, señal {senal}, {etiqueta_eta(eta)}")
    ax.set_xlabel(r"$\lambda$")
    ax.legend(loc="upper left", fontsize=8)
guardar(fig, "16_penalizaciones_error_lambda.png")

# 17. Factores de filtro generalizados en el lambda de discrepancia (pulso)
fig, ax = plt.subplots(figsize=(7.5, 3.8))
b = A @ SENALES["pulso"] + datos.ruido(A @ SENALES["pulso"], eta, z)
for nombre, L in PENALIZACIONES.items():
    g = gsvds[nombre]
    lam_d, _ = G.discrepancia_general(A, L, b, LAMBDAS, np.linalg.norm(datos.ruido(A @ SENALES["pulso"], eta, z)), TAU, gsvd=g)
    f = G.factores_filtro_general(g, lam_d)
    ax.semilogy(np.arange(1, N + 1), np.maximum(f, 1e-20), color=COLOR_L[nombre], label=f"L = {nombre}, λ = {lam_d:.1e}")
ax.axhline(0.5, color=TENUE, lw=1, ls=":")
ax.set_ylim(1e-12, 3)
ax.set_xlabel("índice i (direcciones generalizadas, ordenadas por $c_i$ decreciente)")
ax.set_title(rf"Factores de filtro $f_i=\gamma_i^2/(\gamma_i^2+\lambda)$ por discrepancia, pulso, {etiqueta_eta(eta)}")
ax.legend(loc="lower left", fontsize=8)
guardar(fig, "17_penalizaciones_filtros.png")
