# Declaración de uso de inteligencia artificial

| Campo | Detalle |
|---|---|
| Herramienta | Claude (Anthropic), modelo Claude Opus 5.5, en claude.ai |
| Finalidad | En la elaboración de este trabajo se utilizó como apoyo en las siguientes tareas: comprensión del tema, investigación y búsqueda de fuentes, mejora de la redacción, aplicación de formato Markdown y elaboración de gráficos. Las fuentes identificadas con ayuda de la herramienta fueron revisadas y verificadas por los autores, y los gráficos fueron validados contra los datos originales. Las ideas, el análisis, las conclusiones y la responsabilidad final sobre el contenido son de los autores, quienes revisaron y editaron todo el material generado con apoyo de la IA. |

## Ampliación posterior a la entrega parcial (ejemplo práctico de imagen y audio)

| Campo | Detalle |
|---|---|
| Herramienta | Claude (Anthropic), modelo Claude Opus 5.5, noviembre de 2026 |
| Finalidad | Implementar el ejemplo práctico pedido en la retroalimentación del docente (reconstrucción de una imagen y de un audio) reutilizando el paquete `reconstruccion/`. |
| Instrucción relevante | Continuar el repositorio entregado agregando un ejemplo práctico con una imagen o un sonido, con el mismo método y la misma regla de selección del parámetro. |
| Partes asistidas | `reconstruccion/separable.py`, `reconstruccion/circulante.py`, `scripts/05_imagen.py`, `scripts/06_audio.py`, `tests/test_aplicaciones.py`, sección 9 de `docs/desarrollo_matematico.md`, sección 7 de `docs/resultados.md`, versión 2 de `docs/guion_video.md` y los cambios del README. |
| Verificación | 22 pruebas automáticas forman las matrices grandes (Kronecker y circulante) en casos pequeños y comprueban que las fórmulas rápidas dan los mismos valores singulares y las mismas soluciones que `metodos.py` y que el sistema aumentado resuelto con QR. Se revisaron visualmente todas las figuras. Los DOI de las referencias nuevas [9]–[11] se comprobaron en las páginas de las editoriales. **Pendiente: revisión de las demostraciones de la sección 9 por los tres integrantes.** |

## Entrega final (ampliación, estimación del ruido, revisión matemática y material del video)

| Campo | Detalle |
|---|---|
| Herramienta | Claude (Anthropic), modelo Claude Fable 5.1, en claude.ai, 17 de noviembre de 2026 |
| Finalidad | Completar la entrega final: la ampliación elegida por el grupo (penalización de diferencias $\lambda\lVert Lx\rVert^2$), la estimación del nivel de ruido a partir de los datos, la incorporación de las correcciones de una revisión matemática hecha sobre la entrega parcial, y el guion y las diapositivas del video. |
| Partes asistidas | `reconstruccion/general.py`, `scripts/07_diagnosticos.py`, `scripts/08_penalizacion_diferencias.py`, `scripts/09_estimacion_delta.py`, `tests/test_general.py`, secciones 3, 4, 5, 6, 7, 8.4, 10 y 11 de `docs/desarrollo_matematico.md`, `docs/resultados.md`, `docs/guion_video.md` (versión 3), `docs/diapositivas/`, README y `docs/organizacion.md`. |
| Revisión matemática | Después de la entrega parcial, un segundo agente de IA sin acceso a la redacción revisó las demostraciones y recalculó cada cifra citada (incluida una SVD con 80 dígitos de precisión). Confirmó las demostraciones y encontró tres errores de interpretación que corregimos en esta entrega: el $\kappa_2(A)\approx2\times10^{16}$ no es el número de condición verdadero sino un artefacto del redondeo ($A$ es invertible en aritmética exacta pero numéricamente singular); el residuo de la pseudoinversa completa se explica por el redondeo en la SVD, no en el producto $Ax$; y el error sí baja al reducir el ruido, aunque despacio. También corrigió detalles menores (meseta de $\sigma_i$ desde $i=98$, monotonía numérica del residuo de la TSVD, ruido esférico y no gaussiano). Todas las cifras que antes no eran reproducibles se incorporaron a `scripts/07_diagnosticos.py`. |
| Verificación | 61 pruebas automáticas. La GSVD calculada reconstruye $A$ y $L$ con error $10^{-14}$ y la fórmula espectral generalizada coincide con el sistema aumentado resuelto por QR. El límite $\lambda\to\infty$ con $L_1$ se comprobó contra el valor analítico ($\operatorname{media}(b)$). Los resultados se regeneraron desde cero en un entorno limpio. **Pendiente del grupo:** cada integrante rehace a mano la Proposición 5 y el Teorema 3 (sección 10 del desarrollo) y revisa las respuestas de la sección 11 antes de la ronda de preguntas. |
