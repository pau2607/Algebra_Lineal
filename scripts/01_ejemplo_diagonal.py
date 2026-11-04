import numpy as np
from comun import C1, C2, C3, TINTA, guardar, guardar_tabla, plt

from reconstruccion import metodos as M
from reconstruccion.metricas import residuo

A = np.diag([1.0, 1e-3])
x_true = np.array([1.0, 1.0])
e = np.array([0.0, 0.01])
b = A @ x_true + e
svd = M.SVD.de(A)
lam = 1e-4

sols = {
    "Pseudoinversa": M.pseudoinversa_completa(svd, b),
    "SVD truncada (k=1)": M.svd_truncada(svd, b, 1),
    "Tikhonov (λ=1e-4)": M.tikhonov_espectral(svd, b, lam),
}

# (x, error, residuo) a mano
mano = {
    "Pseudoinversa": (np.array([1.0, 11.0]), 10.0, 0.0),
    "SVD truncada (k=1)": (np.array([1.0, 0.0]), 1.0, 0.011),
    "Tikhonov (λ=1e-4)": (np.array([1 / 1.0001, 0.000011 / 0.000101]), 0.8911, None),
}

print("Ejemplo diagonal")
filas = []
for nombre, x in sols.items():
    err = np.linalg.norm(x - x_true)
    res = residuo(A, x, b)
    x_m, err_m, res_m = mano[nombre]

    assert np.allclose(x, x_m, atol=1e-10), nombre
    assert abs(err - err_m) < 5e-4, nombre
    if res_m is not None:
        assert abs(res - res_m) < 1e-12, nombre

    filas.append({
        "método": nombre,
        "x1": float(x[0]),
        "x2": float(x[1]),
        "error ||x - x_true||": float(err),
        "residuo ||Ax - b||": float(res),
        "norma ||x||": float(np.linalg.norm(x)),
    })
    print(f"  {nombre:22s} x = ({x[0]:.5f}, {x[1]:.5f})  error = {err:.4f}  residuo = {res:.3g}")

f = M.factores_filtro(svd.sigma, lam)
print(f"  factores de filtro f_i(1e-4) = {np.array2string(f, precision=6)}")
print(f"  amplificación del ruido en la dirección 2: 1/sigma_2 = {1 / svd.sigma[1]:.0f}")
print("  Todas las cifras coinciden con los cálculos a mano.")
guardar_tabla(filas, "01_ejemplo_diagonal")

fig, ax = plt.subplots(figsize=(6.4, 3.6))
pos = np.arange(2)
ancho = 0.2
series = [
    ("Verdad", x_true, TINTA),
    ("Pseudoinversa (residuo 0)", sols["Pseudoinversa"], C2),
    ("SVD truncada k=1", sols["SVD truncada (k=1)"], C1),
    ("Tikhonov λ=1e-4", sols["Tikhonov (λ=1e-4)"], C3),
]
for j, (nom, vals, color) in enumerate(series):
    barras = ax.bar(pos + (j - 1.5) * ancho, vals, ancho * 0.9, color=color, label=nom)
    for bar, v in zip(barras, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.15, f"{v:.3g}",
                ha="center", va="bottom", fontsize=8, color=TINTA)

ax.set_xticks(pos, ["Componente 1", "Componente 2"])
ax.set_ylabel("valor reconstruido")
ax.set_title("Residuo cero no garantiza una buena reconstrucción")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
ax.set_ylim(0, 12.5)
guardar(fig, "01_ejemplo_diagonal.png")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
ax.set_ylim(0, 12.5)
guardar(fig, "01_ejemplo_diagonal.png")
