# Guion del video final (versión 2, con ejemplo práctico)

Duración objetivo: 9:50 (tolerancia 9–11 min, incluye identificación y créditos).

**Cambios respecto de la versión 1 (retroalimentación de la entrega parcial):**
- Se agrega una demostración con datos que se pueden ver y oír: una imagen y un audio (bloque 9, 1:40).
- Para hacer espacio, se comprimen los bloques del caso 1D. Estabilidad y sesgo pasa al bloque de Tikhonov, y la interpretación se acorta.
- El problema se abre con la imagen de texto borrosa como gancho, y al final se muestra la misma imagen reconstruida.

La asignación de bloques nuevos es una propuesta. Cada integrante habla alrededor de 3 a 3.5 minutos.

| # | Tiempo | Quién | Contenido | Recurso en pantalla |
|---|---|---|---|---|
| 0 | 0:00–0:20 | Paula | Tema, título, curso e integrantes | Portada |
| 1 | 0:20–1:10 | Paula | **Problema.** Gancho: "¿alguien puede leer esto?" (texto desenfocado). Ese es nuestro problema: $b=Ax_{\mathrm{true}}+e$. Qué es $A$ (suavizado gaussiano, $s=0.03$), la señal (zona suave + pulso), el ruido $\eta$. Pregunta guía: ¿cómo reconstruir sin amplificar lo que no sabemos? | `resultados/imagen/texto_observada_eta1e-3.png`, luego `02_senal_y_datos.png` |
| 2 | 1:10–2:00 | Paula | **Ejemplo diagonal.** $x^\dagger=(1,11)$ con residuo 0 y error 10; $1/\sigma_2=1000$. "Un residuo pequeño no certifica una buena reconstrucción". | `01_ejemplo_diagonal.png` + pizarra |
| 3 | 2:00–2:50 | Matías | **SVD y pseudoinversa.** $x^\dagger=\sum (u_i^Tb/\sigma_i)v_i$, solución de mínima norma, término $A^\dagger e$ y cota $\lVert e\rVert/\sigma_r$. Por qué no usamos ecuaciones normales ($\kappa^2$). | Desarrollo §3–4 |
| 4 | 2:50–3:50 | Matías | **Tikhonov.** Función objetivo, unicidad ($z^T(A^TA+\lambda I)z>0$), fórmula espectral, factores de filtro frente al corte de la TSVD. Cota $\lVert e\rVert/2\sqrt\lambda$ y el balance entre sesgo y ruido. | §6–7 + `08_factores_filtro.png` |
| 5 | 3:50–4:40 | Carlos | **Espectro de nuestro $A$.** Gráfico de Picard; rango numérico 89; la pinv da errores de $10^8$–$10^{10}$ con un residuo menor que el ruido. | `03_espectro_picard.png`, `06b_pseudoinversa.png` |
| 6 | 4:40–5:30 | Carlos | **Elegir el parámetro.** Curvas de error, residuo y norma; oráculo frente a discrepancia (qué usa cada uno). | `04_curvas_tikhonov.png` |
| 7 | 5:30–6:20 | Carlos | **Reconstrucciones 1D** por discrepancia para los tres $\eta$; robustez con 50 semillas; qué pasa si $\delta$ está mal estimado. | `06_reconstrucciones.png`, tabla 03 |
| 8 | 6:20–7:00 | Paula | **Interpretación 1D.** El error no baja de ≈ 9 % porque domina el sesgo en los bordes del pulso; sin ruido el mínimo es 5.4 %: información perdida. | tabla 03 sin ruido |
| 9a | 7:00–7:50 | Paula | **Demostración: imagen.** La misma $A$ en las dos direcciones: $B=AXA^T+E$. La matriz completa tendría 65 536 × 65 536, pero su SVD sale de la SVD de $A$ (Kronecker). Mostrar en orden: observada → pseudoinversa (ruido puro) → Tikhonov. El texto vuelve a leerse. | `11_imagen_texto.png`, `13_imagen_ampliacion.png` (o las imágenes de `resultados/imagen/` a pantalla completa) |
| 9b | 7:50–8:40 | Matías | **Demostración: audio.** Reproducir ~5 s de cada uno: original → observado (apagado) → Tikhonov. Explicar con el espectro: Tikhonov devuelve los armónicos hasta ≈ 1 kHz; más arriba la información se perdió. Avisar antes de reproducir la pseudoinversa (o solo mostrarla). | `resultados/audio/melodia_*.wav`, `14_audio.png` (panel central) |
| 10 | 8:40–9:30 | Matías | **Conclusiones y limitaciones.** Truncamiento y Tikhonov controlan la amplificación a costa de sesgo; Tikhonov rinde mejor en la imagen; la selección necesita una buena estimación del ruido; la regularización no recupera lo que la observación destruyó (el pulso, los bordes de las letras, los agudos). Cierre con el texto ya legible. | `resultados/imagen/texto_tikhonov_eta1e-3.png` + diapositiva final |
| 11 | 9:30–9:50 | Carlos | Referencias principales y enlace al repositorio | Créditos |

## Notas de producción

- Los audios están a 8000 Hz: en el video conviene reproducirlos desde el archivo y no grabarlos con micrófono.
- La imagen y el audio que se muestran se eligen sin mirar la original (discrepancia); conviene decirlo, porque es lo que haría alguien con datos reales.
- Pendiente: ensayar con cronómetro, sobre todo el bloque 9 (el audio toma tiempo real de reproducción).
