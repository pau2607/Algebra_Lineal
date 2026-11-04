"""Pasos 2-6 del experimento obligatorio (sección 9 del artículo).

Para s = 0.03 y eta en {1e-4, 1e-3, 1e-2}:
  * construye señal, matriz y datos perturbados;
  * compara pseudoinversa numérica, SVD truncada (todos los k) y Tikhonov
    (malla logarítmica lambda / sigma_1^2 en [1e-10, 1]);
  * mide error relativo, residuo y norma; muestra espectro y coeficientes de Picard;
  * elige parámetros con el principio de discrepancia SIN usar x_true y después
    los compara con el mejor valor oráculo.
"""

import numpy as np
from comun import (C1, C2, C3, C4, COLOR_ETA, ETAS, KS, LAMBDA_REL, N, REJILLA, S, SEMILLA, SUPERFICIE, TAU, TENUE,
                   TINTA, TINTA_2, entorno, etiqueta_eta, guardar, guardar_json, guardar_tabla, plt)

from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import seleccion as SEL
from reconstruccion.metricas import error_relativo, norma, residuo

t = datos.malla(N)
A = datos.matriz_suavizado(t, S)
x_true = datos.senal_verdadera(t)
b_limpio = A @ x_true
z = datos.direccion_ruido(N, SEMILLA)
svd = M.SVD.de(A)
sigma = svd.sigma
lambdas = sigma[0] ** 2 * LAMBDA_REL

print(f"Problema: n = {N}, s = {S}, semilla = {SEMILLA}")
print(f"  sigma_1 = {sigma[0]:.6g}, sigma_n = {sigma[-1]:.3g}, kappa_2(A) = {sigma[0] / sigma[-1]:.3g}")
print(f"  A simétrica: {np.allclose(A, A.T)}; filas suman 1: {np.allclose(A.sum(axis=1), 1)}")


dif_max = 0.0
b_prueba = b_limpio + datos.ruido(b_limpio, 1e-3, z)
for lam in lambdas[::5]:
    xe = M.tikhonov_espectral(svd, b_prueba, lam)
    xa = M.tikhonov_aumentado(A, b_prueba, lam)
    dif_max = max(dif_max, np.linalg.norm(xe - xa) / np.linalg.norm(xa))
print(f"  espectral vs aumentado (QR): máx. diferencia relativa = {dif_max:.2e}")

resultados = {}
filas_resumen = []
filas_curva = []

for eta in ETAS:
    e = datos.ruido(b_limpio, eta, z)
    b = b_limpio + e
    delta = np.linalg.norm(e)   # cota de ruido: en este experimento se toma delta = ||e||
    print(f"\neta = {eta:g}  (||e|| = {delta:.4g}, ||e||/||Ax|| = {delta / np.linalg.norm(b_limpio):.3g})")

   
    x_pinv, rango, umbral = M.pseudoinversa_numerica(A, b)
    x_pinv_c = M.pseudoinversa_completa(svd, b)

    
    X_k = [M.svd_truncada(svd, b, int(k)) for k in KS]
    err_k = np.array([error_relativo(x, x_true) for x in X_k])
    res_k = np.array([residuo(A, x, b) for x in X_k])
    nor_k = np.array([norma(x) for x in X_k])

    X_l = [M.tikhonov_espectral(svd, b, lam) for lam in lambdas]
    err_l = np.array([error_relativo(x, x_true) for x in X_l])
    res_l = np.array([residuo(A, x, b) for x in X_l])
    nor_l = np.array([norma(x) for x in X_l])

    lam_d, i_d = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta, TAU)
    k_d, j_d = SEL.discrepancia_tsvd(A, b, svd, KS, delta, TAU)
    
    i_o, _ = SEL.oraculo(X_l, x_true)
    j_o, _ = SEL.oraculo(X_k, x_true)

    resultados[eta] = dict(b=b, e=e, delta=delta, x_pinv=x_pinv, x_pinv_c=x_pinv_c, rango=rango, umbral=umbral,
                           err_k=err_k, res_k=res_k, nor_k=nor_k, err_l=err_l, res_l=res_l, nor_l=nor_l,
                           i_d=i_d, j_d=j_d, i_o=i_o, j_o=j_o,
                           x_tik_d=X_l[i_d], x_tik_o=X_l[i_o], x_tsvd_d=X_k[j_d], x_tsvd_o=X_k[j_o])

    def fila(metodo, parametro, x):
        return {"eta": eta, "método": metodo, "parámetro": parametro,
                "error relativo": error_relativo(x, x_true), "residuo": residuo(A, x, b),
                "residuo / delta": residuo(A, x, b) / delta, "norma": norma(x)}

    filas_resumen += [
        fila("Pseudoinversa completa (128 σ)", "k = 128, sin umbral", x_pinv_c),
        fila("Pseudoinversa numérica (scipy)", f"rango {rango}, umbral {umbral:.2e}", x_pinv),
        fila("SVD truncada · discrepancia", f"k = {k_d}", X_k[j_d]),
        fila("SVD truncada · oráculo", f"k = {KS[j_o]}", X_k[j_o]),
        fila("Tikhonov · discrepancia", f"λ/σ1² = {LAMBDA_REL[i_d]:.2e}", X_l[i_d]),
        fila("Tikhonov · oráculo", f"λ/σ1² = {LAMBDA_REL[i_o]:.2e}", X_l[i_o]),
    ]
    for lr, er, rs, nr in zip(LAMBDA_REL, err_l, res_l, nor_l):
        filas_curva.append({"eta": eta, "lambda/sigma1^2": lr, "error relativo": er, "residuo": rs, "norma": nr})

    print(f"  pinv numérica: rango {rango} (umbral {umbral:.2e}), error relativo = {error_relativo(x_pinv, x_true):.3g}")
    print(f"  pinv completa: error relativo = {error_relativo(x_pinv_c, x_true):.3g}")
    print(f"  TSVD   discrepancia k = {k_d:3d}  error = {err_k[j_d]:.4f} | oráculo k = {KS[j_o]:3d}  error = {err_k[j_o]:.4f}")
    print(f"  Tikh.  discrepancia λ/σ1² = {LAMBDA_REL[i_d]:.2e} error = {err_l[i_d]:.4f} | "
          f"oráculo λ/σ1² = {LAMBDA_REL[i_o]:.2e} error = {err_l[i_o]:.4f}")

guardar_tabla(filas_resumen, "02_resumen_metodos")
guardar_tabla(filas_curva, "02_curvas_tikhonov")
guardar_json({
    "n": N, "s": S, "semilla": SEMILLA, "etas": ETAS, "tau": TAU,
    "malla_lambda_rel": {"min": float(LAMBDA_REL[0]), "max": float(LAMBDA_REL[-1]), "puntos": len(LAMBDA_REL)},
    "sigma_1": float(sigma[0]), "sigma_n": float(sigma[-1]), "kappa_2": float(sigma[0] / sigma[-1]),
    "pinv_scipy": {"rtol": "max(m,n)*eps (por defecto)", "umbral_absoluto": float(resultados[ETAS[0]]["umbral"]),
                   "rango_numerico": int(resultados[ETAS[0]]["rango"])},
    "numpy_matrix_rank_por_defecto": int(np.linalg.matrix_rank(A)),
    "valores_singulares": [float(v) for v in sigma],
    "entorno": entorno(),
}, "configuracion.json")


# Figuras


# Señal, datos limpios y perturbados
fig, axs = plt.subplots(1, 3, figsize=(14, 3.6))
axs[0].plot(t, x_true, color=TINTA, label=r"señal $x_{true}$")
axs[0].plot(t, b_limpio, color=C1, label=r"datos limpios $Ax_{true}$")
axs[0].set_title(f"Señal y datos limpios (s = {S})")
axs[0].set_xlabel("t")
axs[0].legend(loc="upper right")
axs[1].plot(t, b_limpio, color=C1, lw=2, label=r"$Ax_{true}$")
axs[1].plot(t, resultados[1e-2]["b"], color=C3, lw=1.2, marker="o", ms=3, label=r"$b$, $\eta = 10^{-2}$")
axs[1].set_xlim(0.45, 0.95)
axs[1].set_title(r"Datos perturbados (ampliación), $\eta = 10^{-2}$")
axs[1].set_xlabel("t")
axs[1].legend(loc="upper right")
for eta in ETAS[::-1]:
    axs[2].semilogy(t, np.abs(resultados[eta]["e"]), color=COLOR_ETA[eta], lw=1.2, label=etiqueta_eta(eta))
axs[2].set_title(r"Magnitud del ruido $|e_j|$ (misma dirección $z$)")
axs[2].set_xlabel("t")
axs[2].legend(loc="lower right", ncol=3, fontsize=8)
guardar(fig, "02_senal_y_datos.png")

# Espectro singular y coeficientes de Picard
fig, ax = plt.subplots(figsize=(7.5, 4.2))
i = np.arange(1, N + 1)
ax.semilogy(i, sigma, color=TINTA, lw=2, label=r"$\sigma_i$")
for eta in ETAS:
    c = np.abs(svd.coeficientes(resultados[eta]["b"]))
    ax.semilogy(i, c, "o", ms=3, color=COLOR_ETA[eta], label=rf"$|u_i^Tb|$, {etiqueta_eta(eta)}")
ax.semilogy(i, np.abs(svd.coeficientes(b_limpio)), "-", lw=1, color=TENUE, label=r"$|u_i^TAx_{true}|$ (sin ruido)")
ax.axhline(resultados[ETAS[0]]["umbral"], color=C4, lw=1.2, ls="--",
           label=f"umbral de pinv (scipy) ≈ {resultados[ETAS[0]]['umbral']:.1e}")
ax.set_xlabel("índice i")
ax.set_title("Espectro singular y coeficientes de Picard")
ax.set_ylim(1e-19, 10)
ax.legend(loc="lower left", ncol=1, frameon=True, facecolor=SUPERFICIE, edgecolor=REJILLA, framealpha=1)
guardar(fig, "03_espectro_picard.png")

# Curvas de Tikhonov: error, residuo y norma frente a lambda
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
for eta in ETAS:
    r = resultados[eta]
    col = COLOR_ETA[eta]
    axs[0].loglog(LAMBDA_REL, r["err_l"], color=col, label=etiqueta_eta(eta))
    axs[0].plot(LAMBDA_REL[r["i_o"]], r["err_l"][r["i_o"]], "o", color=col, ms=8, mfc="white", mew=2)
    axs[0].plot(LAMBDA_REL[r["i_d"]], r["err_l"][r["i_d"]], "s", color=col, ms=8)
    axs[1].loglog(LAMBDA_REL, r["res_l"], color=col, label=etiqueta_eta(eta))
    axs[1].axhline(TAU * r["delta"], color=col, lw=1, ls="--")
    axs[1].plot(LAMBDA_REL[r["i_d"]], r["res_l"][r["i_d"]], "s", color=col, ms=8)
    axs[2].loglog(LAMBDA_REL, r["nor_l"], color=col, label=etiqueta_eta(eta))
axs[0].set_title("Error relativo de reconstrucción")
axs[1].set_title(r"Residuo $\|Ax_\lambda - b\|$  (-- : $\tau\delta$)")
axs[2].set_title(r"Norma $\|x_\lambda\|$")
for ax in axs:
    ax.set_xlabel(r"$\lambda / \sigma_1^2$")
axs[0].plot([], [], "o", color=TENUE, mfc="white", mew=2, label="oráculo (usa $x_{true}$)")
axs[0].plot([], [], "s", color=TENUE, label="discrepancia (sin $x_{true}$)")
axs[0].legend(loc="upper right", fontsize=8)
axs[2].legend(loc="upper right", fontsize=8)
guardar(fig, "04_curvas_tikhonov.png")

# Curvas de SVD truncada
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
for eta in ETAS:
    r = resultados[eta]
    col = COLOR_ETA[eta]
    axs[0].semilogy(KS, r["err_k"], color=col, label=etiqueta_eta(eta))
    axs[0].plot(KS[r["j_o"]], r["err_k"][r["j_o"]], "o", color=col, ms=8, mfc="white", mew=2)
    axs[0].plot(KS[r["j_d"]], r["err_k"][r["j_d"]], "s", color=col, ms=8)
    axs[1].semilogy(KS, r["res_k"], color=col, label=etiqueta_eta(eta))
    axs[1].axhline(TAU * r["delta"], color=col, lw=1, ls="--")
axs[0].axvline(resultados[ETAS[0]]["rango"], color=TENUE, lw=1, ls=":")
axs[0].text(resultados[ETAS[0]]["rango"] + 2, 10, "rango numérico\nde pinv", fontsize=8, color=TINTA_2)
axs[0].set_title("SVD truncada: error relativo")
axs[1].set_title(r"SVD truncada: residuo  (-- : $\tau\delta$)")
for ax in axs:
    ax.set_xlabel("k (número de valores singulares conservados)")
axs[0].plot([], [], "o", color=TENUE, mfc="white", mew=2, label="oráculo")
axs[0].plot([], [], "s", color=TENUE, label="discrepancia")
axs[0].legend(loc="upper left", fontsize=8)
guardar(fig, "05_curvas_tsvd.png")

# Reconstrucciones elegidas por discrepancia
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
for ax, eta in zip(axs, ETAS):
    r = resultados[eta]
    ax.plot(t, x_true, color=TINTA, lw=2.4, label="verdad")
    ax.plot(t, r["x_tsvd_d"], color=C1, lw=1.6, label=f"SVD truncada, k = {KS[r['j_d']]}")
    ax.plot(t, r["x_tik_d"], color=C2, lw=1.6, label=f"Tikhonov, λ/σ1² = {LAMBDA_REL[r['i_d']]:.1e}")
    ax.set_title(f"{etiqueta_eta(eta)} — parámetros por discrepancia")
    ax.set_xlabel("t")
    ax.legend(loc="upper right", fontsize=8)
axs[0].set_ylim(-0.3, 1.3)
guardar(fig, "06_reconstrucciones.png")

# La pseudoinversa numérica: amplificación del ruido
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [2, 1]})
x_p = resultados[ETAS[0]]["x_pinv"]
axs[0].plot(t, x_p, color=C1, lw=1.0, label=f"pinv numérica, {etiqueta_eta(ETAS[0])}")
axs[0].plot(t, x_true, color=TINTA, lw=2, label=r"verdad ($\max|x| = 1$)")
axs[0].set_title(f"Pseudoinversa numérica (rango {resultados[ETAS[0]]['rango']} de {N}), el ruido más pequeño")
axs[0].set_xlabel("t")
axs[0].legend(loc="lower left", fontsize=8)
normas = [np.linalg.norm(resultados[eta]["x_pinv"]) for eta in ETAS]
barras = axs[1].bar([etiqueta_eta(eta) for eta in ETAS], normas, 0.6, color=[COLOR_ETA[eta] for eta in ETAS])
axs[1].axhline(np.linalg.norm(x_true), color=TINTA, lw=1.5, ls="--")
axs[1].text(-0.4, np.linalg.norm(x_true) * 2.5, r"$\|x_{true}\| \approx 5$", fontsize=8, color=TINTA_2)
for barra, v in zip(barras, normas):
    axs[1].text(barra.get_x() + barra.get_width() / 2, v * 1.6, f"{v:.1e}", ha="center", fontsize=8, color=TINTA)
axs[1].set_yscale("log")
axs[1].set_ylim(1, 1e12)
axs[1].set_title(r"$\|x_{pinv}\|$ crece en proporción a $\eta$")
axs[1].grid(axis="x", visible=False)
guardar(fig, "06b_pseudoinversa.png")

# Curva L (exploratoria: solo se dibuja, no se usa para elegir)
fig, ax = plt.subplots(figsize=(6, 4.4))
for eta in ETAS:
    r = resultados[eta]
    ax.loglog(r["res_l"], r["nor_l"], color=COLOR_ETA[eta], label=etiqueta_eta(eta))
    ax.plot(r["res_l"][r["i_d"]], r["nor_l"][r["i_d"]], "s", color=COLOR_ETA[eta], ms=7)
    ax.plot(r["res_l"][r["i_o"]], r["nor_l"][r["i_o"]], "o", color=COLOR_ETA[eta], ms=7, mfc="white", mew=2)
ax.set_xlabel(r"$\|Ax_\lambda - b\|$")
ax.set_ylabel(r"$\|x_\lambda\|$")
ax.set_title("Curva L de Tikhonov (exploratoria)")
ax.plot([], [], "o", color=TENUE, mfc="white", mew=2, label="oráculo")
ax.plot([], [], "s", color=TENUE, label="discrepancia")
ax.legend(loc="upper right", fontsize=8)
guardar(fig, "07_curva_L.png")

# Factores de filtro
fig, ax = plt.subplots(figsize=(7.5, 3.8))
eta = 1e-3
r = resultados[eta]
lam_d = lambdas[r["i_d"]]
ax.semilogy(i, M.factores_filtro(sigma, lam_d), color=C2, label=rf"Tikhonov $f_i(\lambda)$, λ/σ1² = {LAMBDA_REL[r['i_d']]:.1e}")
ax.semilogy(i, (i <= KS[r["j_d"]]).astype(float) + 1e-300, color=C1, drawstyle="steps-mid",
            label=f"SVD truncada (corte en k = {KS[r['j_d']]})")
ax.set_ylim(1e-20, 3)
ax.set_xlabel("índice i")
ax.set_title(rf"Factores de filtro elegidos por discrepancia, {etiqueta_eta(eta)}")
ax.legend(loc="lower left", fontsize=8)
guardar(fig, "08_factores_filtro.png")
