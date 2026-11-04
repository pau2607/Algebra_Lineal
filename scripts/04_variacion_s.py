import numpy as np
from comun import COLOR_ETA, LAMBDA_REL, N, SEMILLA, TAU, TINTA, guardar, guardar_tabla, plt

from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import seleccion as SEL
from reconstruccion.metricas import error_relativo

ETA = 1e-3
SS = [0.02, 0.03, 0.05]
COLS = dict(zip(SS, COLOR_ETA.values()))
t = datos.malla(N)
x_true = datos.senal_verdadera(t)
z = datos.direccion_ruido(N, SEMILLA)

filas = []
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for s in SS:
    A = datos.matriz_suavizado(t, s)
    b_limpio = A @ x_true
    e = datos.ruido(b_limpio, ETA, z)
    b = b_limpio + e
    svd = M.SVD.de(A)
    lambdas = svd.sigma[0] ** 2 * LAMBDA_REL
    delta = np.linalg.norm(e)
    _, rango, _ = M.pseudoinversa_numerica(A, b)
    X = [M.tikhonov_espectral(svd, b, lam) for lam in lambdas]
    _, i_d = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta, TAU)
    i_o, err_o = SEL.oraculo(X, x_true)
    nivel = delta / np.sqrt(N)
    filas.append({
        "s": s, "kappa_2(A)": float(svd.sigma[0] / svd.sigma[-1]), "rango pinv": rango,
        "σ_i > ||e||/√n": int((svd.sigma > nivel).sum()),
        "λ/σ1² discrepancia": float(LAMBDA_REL[i_d]),
        "error discrepancia": error_relativo(X[i_d], x_true),
        "λ/σ1² oráculo": float(LAMBDA_REL[i_o]), "error oráculo": err_o,
    })
    axs[0].semilogy(np.arange(1, N + 1), svd.sigma, color=COLS[s], label=f"s = {s}")
    axs[1].plot(t, X[i_d], color=COLS[s], lw=1.6, label=f"s = {s} (discrepancia)")
    print(f"  s = {s}: rango pinv {rango}, error discrepancia {filas[-1]['error discrepancia']:.4f}, "
          f"oráculo {err_o:.4f}")

axs[0].axhline(nivel, color=TINTA, lw=1, ls="--")
axs[0].set_title(r"Espectro singular según s  (-- : $\|e\|/\sqrt{n}$)")
axs[0].set_xlabel("índice i")
axs[0].legend(loc="upper right")
axs[1].plot(t, x_true, color=TINTA, lw=2.4, label="verdad")
axs[1].set_title(r"Tikhonov por discrepancia, $\eta = 10^{-3}$")
axs[1].set_xlabel("t")
axs[1].legend(loc="upper right", fontsize=8)
guardar(fig, "09_variacion_s.png")
guardar_tabla(filas, "04_variacion_s")
