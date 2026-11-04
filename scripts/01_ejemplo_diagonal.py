"""Paso 1 del experimento: ejemplo diagonal de dos coordenadas (secciones 3.1, 4 y 5).

A = diag(1, 1e-3), x_true = (1, 1), e = (0, 0.01), b = (1, 0.011).
Se calculan pseudoinversa, SVD truncada con k = 1 y Tikhonov con lambda = 1e-4,
y se comparan con los valores obtenidos a mano.
"""

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

soluciones = {
    "Pseudoinversa": M.pseudoinversa_completa(svd, b),
    "SVD truncada (k=1)": M.svd_truncada(svd, b, 1),
    "Tikhonov (λ=1e-4)": M.tikhonov_espectral(svd, b, lam),
}


a_mano = {
    "Pseudoinversa": (np.array([1.0, 11.0]), 10.0, 0.0),
    "SVD truncada (k=1)": (np.array([1.0, 0.0]), 1.0, 0.011),
    "Tikhonov (λ=1e-4)": (np.array([1 / 1.0001, 0.000011 / 0.000101]), 0.8911, None),
}

print("Ejemplo diagonal")
filas = []
for nombre, x in soluciones.items():
    err = np.linalg.norm(x - x_true)
    res = residuo(A, x, b)
    x_ref, err_ref, res_ref = a_mano[nombre]
    assert np.allclose(x, x_ref, atol=1e-10), nombre
    assert abs(err - err_ref) < 5e-4, nombre
    if res_ref is not None:
        assert abs(res - res_ref) < 1e-12, nombre
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
ancho = 0.26
posiciones = np.arange(2)
series = [("Verdad", x_true, TINTA), ("Pseudoinversa (residuo 0)", soluciones["Pseudoinversa"], C2),
          ("SVD truncada k=1", soluciones["SVD truncada (k=1)"], C1),
          ("Tikhonov λ=1e-4", soluciones["Tikhonov (λ=1e-4)"], C3)]
ancho = 0.2
for j, (nombre, valores, color) in enumerate(series):
    barras = ax.bar(posiciones + (j - 1.5) * ancho, valores, ancho * 0.9, color=color, label=nombre)
    for barra, v in zip(barras, valores):
        ax.text(barra.get_x() + barra.get_width() / 2, v + 0.15, f"{v:.3g}", ha="center",
                va="bottom", fontsize=8, color=TINTA)
ax.set_xticks(posiciones, ["Componente 1", "Componente 2"])
ax.set_ylabel("valor reconstruido")
ax.set_title("Residuo cero no garantiza una buena reconstrucción")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
ax.set_ylim(0, 12.5)
guardar(fig, "01_ejemplo_diagonal.png")
