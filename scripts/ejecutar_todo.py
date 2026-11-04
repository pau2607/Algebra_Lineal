import runpy
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))

SCRIPTS = [
    "01_ejemplo_diagonal.py",
    "02_experimento_principal.py",
    "03_robustez.py",
    "04_variacion_s.py",
]

for nombre in SCRIPTS:
    print(f"\n=== {nombre} ===")
    runpy.run_path(str(AQUI / nombre), run_name="__main__")

print("\nListo: resultados en resultados/figuras y resultados/tablas")
