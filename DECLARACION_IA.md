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
| Partes asistidas | `reconstruccion/separable.py`, `reconstruccion/circulante.py`, `scripts/05_imagen.py`, `scripts/06_audio.py`, `tests/test_aplicaciones.py`, sección 9 de `docs/desarrollo_matematico.md`, sección 7 de `docs/resultados_preliminares.md`, versión 2 de `docs/guion_video.md` y los cambios del README. |
| Verificación | 22 pruebas automáticas forman las matrices grandes (Kronecker y circulante) en casos pequeños y comprueban que las fórmulas rápidas dan los mismos valores singulares y las mismas soluciones que `metodos.py` y que el sistema aumentado resuelto con QR. Se revisaron visualmente todas las figuras. Los DOI de las referencias nuevas [9]–[11] se comprobaron en las páginas de las editoriales. **Pendiente: revisión de las demostraciones de la sección 9 por los tres integrantes.** |
