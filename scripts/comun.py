"""Configuración compartida por los scripts: parámetros del experimento, rutas y estilo de figuras."""

import json
import platform
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402

# --------------------------------------------------------------------------
# Parámetros del experimento obligatorio (sección 9 del artículo)
# --------------------------------------------------------------------------
N = 128                       # número de puntos
S = 0.03                      # ancho del suavizado gaussiano
ETAS = [1e-4, 1e-3, 1e-2]     # niveles de ruido relativo
SEMILLA = 2026                # semilla del vector gaussiano z (la misma dirección para los tres eta)
TAU = 1.1                     # factor de seguridad del principio de discrepancia
LAMBDA_REL = np.logspace(-10, 0, 61)  # malla de lambda / sigma_1^2 (6 puntos por década)
KS = np.arange(1, N + 1)      # valores de k para la SVD truncada

DIR_FIG = RAIZ / "resultados" / "figuras"
DIR_TAB = RAIZ / "resultados" / "tablas"
DIR_FIG.mkdir(parents=True, exist_ok=True)
DIR_TAB.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Estilo de figuras: superficie clara, tinta neutra, colores categóricos en orden fijo
# --------------------------------------------------------------------------
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
TENUE = "#898781"
REJILLA = "#e1e0d9"
EJE = "#c3c2b7"
SUPERFICIE = "#fcfcfb"
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # azul, naranja, aqua, amarillo
COLOR_ETA = {1e-4: C1, 1e-3: C2, 1e-2: C3}

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE,
    "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE,
    "figure.dpi": 110,
    "savefig.dpi": 160,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.edgecolor": EJE,
    "axes.labelcolor": TINTA_2,
    "axes.titlecolor": TINTA,
    "xtick.color": TENUE,
    "ytick.color": TENUE,
    "xtick.labelcolor": TINTA_2,
    "ytick.labelcolor": TINTA_2,
    "axes.grid": True,
    "grid.color": REJILLA,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "legend.labelcolor": TINTA_2,
})


def etiqueta_eta(eta: float) -> str:
    return rf"$\eta = 10^{{{int(np.log10(eta))}}}$"


def guardar(fig, nombre: str):
    ruta = DIR_FIG / nombre
    fig.savefig(ruta, bbox_inches="tight")
    plt.close(fig)
    print(f"  figura -> {ruta.relative_to(RAIZ)}")


def guardar_tabla(filas, nombre: str, columnas=None):
    """Guarda una lista de diccionarios como CSV y como tabla Markdown."""
    import csv

    columnas = columnas or list(filas[0].keys())
    with open(DIR_TAB / f"{nombre}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columnas)
        w.writeheader()
        w.writerows(filas)

    def fmt(v):
        if isinstance(v, float):
            return f"{v:.4g}"
        return "—" if v is None else str(v)

    lineas = ["| " + " | ".join(columnas) + " |", "|" + "---|" * len(columnas)]
    lineas += ["| " + " | ".join(fmt(fila[c]) for c in columnas) + " |" for fila in filas]
    (DIR_TAB / f"{nombre}.md").write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"  tabla  -> resultados/tablas/{nombre}.csv / .md")


def entorno():
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }


def guardar_json(obj, nombre: str):
    ruta = RAIZ / "resultados" / nombre
    ruta.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  json   -> {ruta.relative_to(RAIZ)}")
