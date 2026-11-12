"""Ejemplo práctico 1: imagen desenfocada y con ruido.

Modelo:  B = A X A^T + E,  con A = datos.matriz_suavizado (la misma construcción
del experimento principal) aplicada en las dos direcciones de la imagen.
Se reconstruye con pseudoinversa numérica, SVD truncada y Tikhonov, eligiendo
el parámetro por discrepancia (sin usar la imagen original).

Uso:
    python scripts/05_imagen.py                      # gato + texto
    python scripts/05_imagen.py --imagen foto.jpg    # agrega una foto propia
"""

import argparse

import numpy as np
from PIL import Image
from comun import (C1, C2, C3, LAMBDA_REL, RAIZ, REJILLA, SEMILLA, SUPERFICIE, TAU, TENUE,
                   TINTA, TINTA_2, guardar, guardar_tabla, plt)

from reconstruccion import datos
from reconstruccion import separable as SEP

N_IMG = 256
S_PX = 0.03 * 127          # ancho del suavizado en píxeles: el mismo del experimento principal (≈ 3.8)
S_IMG = S_PX / (N_IMG - 1)  # el mismo ancho expresado en la malla t_j = j/(n-1) de 256 puntos
ETAS_IMG = [1e-3, 1e-2]
COLOR_ETA_IMG = {1e-3: C2, 1e-2: C3}
DIR_IMG = RAIZ / "resultados" / "imagen"
DIR_IMG.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------------
# Imágenes de prueba
# --------------------------------------------------------------------------------

def cargar_imagen(ruta, n=N_IMG):
    """Escala de grises, recorte cuadrado centrado, n x n, valores en [0, 1]."""
    im = Image.open(ruta).convert("L")
    w, h = im.size
    lado = min(w, h)
    im = im.crop(((w - lado) // 2, (h - lado) // 2, (w + lado) // 2, (h + lado) // 2))
    im = im.resize((n, n), Image.LANCZOS)
    return np.asarray(im, dtype=float) / 255.0


def imagen_texto(n=N_IMG):
    """Texto negro sobre blanco, generado con matplotlib (determinista)."""
    fig = plt.figure(figsize=(1, 1), dpi=n, facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_facecolor("white")
    ax.text(0.5, 0.80, "ÁLGEBRA", ha="center", va="center", fontsize=8.5, weight="bold", color="black")
    ax.text(0.5, 0.58, "LINEAL", ha="center", va="center", fontsize=8.5, weight="bold", color="black")
    ax.text(0.5, 0.36, "APLICADA", ha="center", va="center", fontsize=6.5, color="black")
    ax.text(0.5, 0.16, "b = Ax + e    UP 2026-2", ha="center", va="center", fontsize=3.6, color="black")
    fig.canvas.draw()
    rgb = np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(float).mean(axis=2) / 255.0
    plt.close(fig)
    return rgb


def guardar_png(X, nombre):
    Image.fromarray(np.round(np.clip(X, 0, 1) * 255).astype(np.uint8), mode="L").save(DIR_IMG / nombre)


parser = argparse.ArgumentParser()
parser.add_argument("--imagen", help="ruta a una foto propia (opcional)")
args, _ = parser.parse_known_args()

imagenes = {
    "gato": cargar_imagen(RAIZ / "recursos" / "gato_256.png"),
    "texto": imagen_texto(),
}
if args.imagen:
    imagenes["propia"] = cargar_imagen(args.imagen)

# --------------------------------------------------------------------------------
# Operador y comprobaciones
# --------------------------------------------------------------------------------

t = datos.malla(N_IMG)
A = datos.matriz_suavizado(t, S_IMG)
op = SEP.OperadorSeparable.de(A, A)
lambdas = op.sigma1 ** 2 * LAMBDA_REL
rango = op.rango_numerico()
KS = np.unique(np.round(np.logspace(0, np.log10(op.S.size), 150)).astype(int))

print(f"Imagen {N_IMG} x {N_IMG}: x tiene {N_IMG ** 2} incógnitas; la matriz de Kronecker sería "
      f"{N_IMG ** 2} x {N_IMG ** 2} y no se forma.")
print(f"  suavizado gaussiano de {S_PX:.2f} px (s = {S_IMG:.4f}) en cada dirección")
print(f"  S_max = {op.S.max():.4g}, S_min = {op.S.min():.3g}, rango numérico (umbral de pinv) = {rango} de {op.S.size}")

B_prueba = op.aplicar(imagenes["gato"]) + SEP.ruido_matricial(op.aplicar(imagenes["gato"]), 1e-2, SEMILLA)
lam_prueba = lambdas[30]
dif = abs(op.residuo(op.tikhonov(B_prueba, lam_prueba), B_prueba) - op.residuo_tikhonov_espectral(B_prueba, lam_prueba))
print(f"  residuo directo vs. fórmula espectral: diferencia {dif:.1e}")

# --------------------------------------------------------------------------------
# Experimento
# --------------------------------------------------------------------------------

resultados = {}
filas = []
for nombre, X in imagenes.items():
    B_limpio = op.aplicar(X)
    for eta in ETAS_IMG:
        E = SEP.ruido_matricial(B_limpio, eta, SEMILLA)
        B = B_limpio + E
        delta = float(np.linalg.norm(E))

        X_pinv = op.svd_truncada(B, rango)
        X_pc = op.pseudoinversa_completa(B)
        k_d, _ = SEP.discrepancia_tsvd(op, B, KS, delta, TAU)
        X_tsvd = op.svd_truncada(B, k_d)
        lam_d, i_d = SEP.discrepancia_tikhonov(op, B, lambdas, delta, TAU)
        X_tik = op.tikhonov(B, lam_d)
        err_l = np.array([np.linalg.norm(op.tikhonov(B, lam) - X) / np.linalg.norm(X) for lam in lambdas])
        i_o = int(np.argmin(err_l))
        X_tik_o = op.tikhonov(B, lambdas[i_o])

        resultados[(nombre, eta)] = dict(X=X, B=B, E=E, pinv=X_pinv, tsvd=X_tsvd, tik=X_tik, k_d=k_d,
                                         i_d=i_d, i_o=i_o, err_l=err_l)

        def fila(metodo, parametro, Xr):
            return {"imagen": nombre, "eta": eta, "método": metodo, "parámetro": parametro,
                    "error relativo": float(np.linalg.norm(Xr - X) / np.linalg.norm(X)),
                    "PSNR (dB)": SEP.psnr(Xr, X), "residuo / delta": op.residuo(Xr, B) / delta}

        filas += [
            fila("Observada (sin reconstruir)", "-", B),
            fila("Pseudoinversa completa", f"{op.S.size} direcciones", X_pc),
            fila("Pseudoinversa numérica", f"rango {rango}", X_pinv),
            fila("SVD truncada · discrepancia", f"k = {k_d}", X_tsvd),
            fila("Tikhonov · discrepancia", f"λ/σ1² = {LAMBDA_REL[i_d]:.2e}", X_tik),
            fila("Tikhonov · oráculo", f"λ/σ1² = {LAMBDA_REL[i_o]:.2e}", X_tik_o),
        ]
        print(f"\n{nombre}, eta = {eta:g}")
        for f in filas[-6:]:
            print(f"  {f['método']:<30} {f['parámetro']:<22} error = {f['error relativo']:.3g}  "
                  f"PSNR = {f['PSNR (dB)']:7.1f} dB")

guardar_tabla(filas, "05_imagen")

# --------------------------------------------------------------------------------
# Figuras
# --------------------------------------------------------------------------------


def mostrar(ax, Xr, X, titulo):
    ax.imshow(np.clip(Xr, 0, 1), cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    if Xr is X:
        ax.set_title(titulo, fontsize=10)
    else:
        err = np.linalg.norm(Xr - X) / np.linalg.norm(X)
        ax.set_title(f"{titulo}\nerror {err:.2g} · PSNR {SEP.psnr(Xr, X):.1f} dB", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for lado in ax.spines.values():
        lado.set_visible(False)


for numero, nombre in [("10", "gato"), ("11", "texto")] + ([("13b", "propia")] if args.imagen else []):
    fig, axs = plt.subplots(len(ETAS_IMG), 5, figsize=(15, 7.4))
    for fila_ax, eta in zip(axs, ETAS_IMG):
        r = resultados[(nombre, eta)]
        mostrar(fila_ax[0], r["X"], r["X"], "original (desconocida)")
        mostrar(fila_ax[1], r["B"], r["X"], rf"observada, $\eta=10^{{{int(np.log10(eta))}}}$")
        mostrar(fila_ax[2], r["pinv"], r["X"], f"pseudoinversa (rango {rango})")
        mostrar(fila_ax[3], r["tsvd"], r["X"], f"SVD truncada, k = {r['k_d']}")
        mostrar(fila_ax[4], r["tik"], r["X"], rf"Tikhonov, $\lambda/\sigma_1^2$ = {LAMBDA_REL[r['i_d']]:.1e}")
    fig.suptitle(f"Desenfoque gaussiano de {S_PX:.1f} px + ruido; parámetros elegidos por discrepancia "
                 "(sin mirar la original)", color=TINTA, fontsize=11)
    fig.tight_layout(h_pad=2.0)
    guardar(fig, f"{numero}_imagen_{nombre}.png")

# Espectro (Picard en 2D) y curvas de error
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
orden = op.orden_espectral
S_ord = op.S.ravel()[orden]
idx = np.unique(np.round(np.logspace(0, np.log10(S_ord.size), 1500)).astype(int)) - 1
axs[0].loglog(idx + 1, S_ord[idx], color=TINTA, lw=2, label=r"$S_{ij}=\sigma_i\sigma_j$ (ordenados)")
for eta in ETAS_IMG:
    r = resultados[("gato", eta)]
    C_ord = np.abs(op.coeficientes(r["B"])).ravel()[orden]
    axs[0].loglog(idx + 1, C_ord[idx], ".", ms=2, color=COLOR_ETA_IMG[eta],
                  label=rf"$|C_{{ij}}|$ (gato), $\eta=10^{{{int(np.log10(eta))}}}$")
    axs[0].axhline(np.linalg.norm(r["E"]) / N_IMG, color=COLOR_ETA_IMG[eta], lw=1, ls="--")
axs[0].plot([], [], "--", color=TENUE, lw=1, label=r"nivel del ruido $\|E\|_F/n$")
axs[0].set_ylim(1e-20, 1e3)
axs[0].set_xlabel("posición en el orden de $S_{ij}$ (de 65 536)")
axs[0].set_title("Espectro del desenfoque 2D y coeficientes de Picard")
axs[0].legend(loc="lower left", fontsize=8, frameon=True, facecolor=SUPERFICIE, edgecolor=REJILLA)

estilos = {"gato": "-", "texto": ":"}
for (nombre, eta), r in resultados.items():
    if nombre not in estilos:
        continue
    col = COLOR_ETA_IMG[eta]
    axs[1].loglog(LAMBDA_REL, r["err_l"], estilos[nombre], color=col,
                  label=rf"{nombre}, $\eta=10^{{{int(np.log10(eta))}}}$")
    axs[1].plot(LAMBDA_REL[r["i_o"]], r["err_l"][r["i_o"]], "o", color=col, ms=7, mfc="white", mew=2)
    axs[1].plot(LAMBDA_REL[r["i_d"]], r["err_l"][r["i_d"]], "s", color=col, ms=7)
axs[1].plot([], [], "o", color=TENUE, mfc="white", mew=2, label="oráculo")
axs[1].plot([], [], "s", color=TENUE, label="discrepancia")
axs[1].set_ylim(0.03, 10)
axs[1].set_xlabel(r"$\lambda/\sigma_1^2$")
axs[1].set_title("Tikhonov: error relativo de la imagen reconstruida")
axs[1].legend(loc="upper left", fontsize=8, ncol=2)
guardar(fig, "12_imagen_espectro_y_error.png")

# Ampliaciones para el video: ojo del gato y última línea del texto
eta = ETAS_IMG[0]
recortes = [("gato", (slice(70, 140), slice(40, 140)), "ojo del gato"),
            ("texto", (slice(200, 240), slice(20, 236)), "línea pequeña del texto")]
fig, axs = plt.subplots(2, 3, figsize=(13, 5.2), gridspec_kw={"height_ratios": [70 / 100, 40 / 216 * 2.2]})
for fila_ax, (nombre, rec, rotulo) in zip(axs, recortes):
    r = resultados[(nombre, eta)]
    for ax, Xr, tit in zip(fila_ax, [r["X"], r["B"], r["tik"]], ["original", "observada", "Tikhonov (discrepancia)"]):
        ax.imshow(np.clip(Xr[rec], 0, 1), cmap="gray", vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        ax.set_title(f"{tit} — {rotulo}", fontsize=9, color=TINTA_2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
fig.suptitle(rf"Ampliación, $\eta=10^{{{int(np.log10(eta))}}}$", color=TINTA)
fig.tight_layout()
guardar(fig, "13_imagen_ampliacion.png")

# Imágenes sueltas para las diapositivas del video
for nombre in imagenes:
    for eta in ETAS_IMG:
        r = resultados[(nombre, eta)]
        sufijo = f"eta1e{int(np.log10(eta))}"
        guardar_png(r["B"], f"{nombre}_observada_{sufijo}.png")
        guardar_png(r["pinv"], f"{nombre}_pseudoinversa_{sufijo}.png")
        guardar_png(r["tik"], f"{nombre}_tikhonov_{sufijo}.png")
    guardar_png(imagenes[nombre], f"{nombre}_original.png")
print("  imágenes sueltas -> resultados/imagen/")
