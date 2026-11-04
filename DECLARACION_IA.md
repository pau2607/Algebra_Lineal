# Declaración de uso de inteligencia artificial

> Revisar y completar antes de entregar. La declaración debe describir el uso real; si el grupo modifica o reescribe partes, debe indicarlo.

| Campo | Detalle |
|---|---|
| Herramienta | Claude (Anthropic), modelo Claude Opus 5.5, en claude.ai |
| Fecha de uso | 4 de noviembre de 2026 |
| Finalidad | Generar una primera versión del código del experimento, de las pruebas automáticas y de los borradores de los documentos de la entrega parcial (desarrollo matemático, resultados, guion, organización y README) |
| Instrucción principal | "Tenemos que hacer la entrega parcial de este proyecto", adjuntando el artículo de referencia del docente y las rúbricas del trabajo |
| Partes en las que se usó | `reconstruccion/`, `scripts/`, `tests/`, `docs/`, `README.md` |
| Partes sin IA | *(completar: p. ej. distribución de responsabilidades, revisión de las demostraciones, decisiones sobre la ampliación)* |

## Procedimientos de verificación

Realizados al generar el material:

1. Las cifras del ejemplo diagonal se compararon con los cálculos a mano del artículo (aserciones en `scripts/01_ejemplo_diagonal.py` y pruebas `test_diagonal_*`).
2. La fórmula espectral de Tikhonov se contrastó con dos cálculos independientes: el sistema aumentado resuelto con QR y las ecuaciones normales.
3. Se comprobaron numéricamente las propiedades demostradas: gradiente nulo en $x_\lambda$, cota de estabilidad, límites en $\lambda$, monotonía del residuo y factores de filtro en $(0,1)$.
4. Se verificó que las funciones de selección no reciben $x_{\mathrm{true}}$.
5. Los números citados en `docs/resultados_preliminares.md` se cotejaron con las tablas generadas por los scripts.
6. Las figuras se revisaron visualmente.

Por realizar por el grupo *(completar y marcar)*:

- [ ] Cada integrante rehízo a mano al menos una demostración de `docs/desarrollo_matematico.md` (indicar cuál).
- [ ] Se ejecutó `python scripts/ejecutar_todo.py` y `python -m pytest -q` en el computador de al menos un integrante y se obtuvieron los mismos resultados.
- [ ] Se revisaron las referencias bibliográficas en sus fuentes originales.
- [ ] Se leyó y entendió cada función del paquete `reconstruccion/` (cada responsable la puede explicar).
