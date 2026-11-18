# Organización del trabajo

## 1. Integrantes y responsabilidades

| Integrante | Responsabilidad principal | Archivos a su cargo | Bloque del video |
|---|---|---|---|
| Paula Romero | Planteamiento, construcción del problema, SVD, pseudoinversa y TSVD | `reconstruccion/datos.py`, `scripts/01_ejemplo_diagonal.py`, parte de pseudoinversa y TSVD de `reconstruccion/metodos.py`, §1 a §5 del desarrollo | 0:00-3:30. Problema, ejemplo diagonal, SVD y pseudoinversa, espectro de A y pinv numérica |
| Carlos Ledo | Tikhonov, estabilidad, sesgo y robustez | parte de Tikhonov de `reconstruccion/metodos.py`, `tests/`, `scripts/03_robustez.py`, `scripts/04_variacion_s.py`, §6 y §7 del desarrollo | 3:30-6:35. Tikhonov y factores de filtro, estabilidad y sesgo, robustez y variación de s |
| Matias Desmonteix | Selección del parámetro, experimento principal e interpretación | `reconstruccion/seleccion.py`, `scripts/02_experimento_principal.py`, `docs/resultados.md`, guion, §8 del desarrollo | 6:35-10:00. Curvas y selección del parámetro, reconstrucciones, interpretación, conclusiones y créditos |

Cada uno habla alrededor de 3 minutos de contenido técnico. Todos tenemos que poder explicar el razonamiento completo, porque la ronda de preguntas es individual.

## 2. Registro de contribuciones (entrega parcial)

| Integrante | Aporte concreto en esta entrega | Evidencia (archivos y commits) |
|---|---|---|
| Paula Romero | Planteé el problema y construí los datos de prueba: implementé la matriz de suavizado gaussiano con $s=0.03$, la señal (zona suave más pulso) y el ruido con semilla registrada para los tres valores de $\eta$. Resolví a mano el ejemplo diagonal y lo reproduje en código. Implementé la pseudoinversa numérica, con su rango y su umbral, y la SVD truncada para todos los $k$. Analicé el espectro de $A$ y los coeficientes de Picard. Redacté las secciones 1 a 5 del desarrollo. | `reconstruccion/datos.py`, `scripts/01_ejemplo_diagonal.py`, `reconstruccion/metodos.py` (pseudoinversa y TSVD), desarrollo §1 a §5. Commits: completar |
| Carlos Ledo | Implementé Tikhonov con la fórmula espectral y lo contrasté con el sistema aumentado resuelto por QR. Escribí las pruebas de `tests/` para verificar los métodos. Hice el estudio de robustez con 50 semillas y la sensibilidad a $\delta$, y la variación del suavizado $s$. Redacté las secciones 6 y 7 del desarrollo (Tikhonov, estabilidad y sesgo). | `reconstruccion/metodos.py` (Tikhonov), `tests/`, `scripts/03_robustez.py`, `scripts/04_variacion_s.py`, desarrollo §6 y §7. Commits: completar |
| Matias Desmonteix | Implementé la selección del parámetro: principio de discrepancia sin usar $x_{\mathrm{true}}$, comparado con el oráculo, y el gráfico de la curva L. Armé el experimento principal que genera las curvas de error, residuo y norma, las reconstrucciones y las tablas. Redacté los resultados preliminares, la sección 8 del desarrollo y el guion del video. | `reconstruccion/seleccion.py`, `scripts/02_experimento_principal.py`, `docs/resultados.md`, guion, desarrollo §8. Commits: completar |

## 3. Registro de contribuciones (entrega final)

> Completar con lo que cada uno hizo de verdad entre la entrega parcial y la final (la rúbrica valora el contenido del aporte y la capacidad de explicarlo). Propuesta según la asignación de la sección 1:

| Integrante | Aporte en la entrega final | Evidencia |
|---|---|---|
| Paula Romero | README final, estimación del nivel de ruido $\hat\delta$ con la cola de Picard (sección 8.4 del desarrollo, script 09), etiqueta `entrega-final`, bloques 1, 2 y 8 del video | `scripts/09_estimacion_delta.py`, `README.md`, desarrollo §8.4. Commits: completar |
| Carlos Ledo | Ejemplos prácticos de imagen y audio (retroalimentación del docente), ampliación con penalización de diferencias (GSVD, Proposición 5 y Teorema 3), pruebas, revisión de reproducibilidad, bloques 5 a 7 del video | `reconstruccion/separable.py`, `reconstruccion/circulante.py`, `reconstruccion/general.py`, `scripts/05`–`08`, `tests/test_general.py`, desarrollo §9–10. Commits: completar |
| Matías Desmonteix | Incorporación de la revisión matemática en el desarrollo (secciones 3 a 8), resultados finales e interpretación, preguntas de defensa (sección 11), guion y edición del video, bloques 3, 4, 9b y 10 | `docs/resultados.md`, desarrollo §3–8 y §11, `docs/guion_video.md`, video. Commits: completar |

## 4. Estado final

| Elemento | Estado |
|---|---|
| Ejemplo diagonal reproducido y verificado | Completo |
| Señal, matriz y ruido ($s=0.03$, tres $\eta$, semilla registrada) | Completo |
| Pseudoinversa numérica (con rango y umbral), TSVD para todos los $k$, Tikhonov en malla logarítmica | Completo |
| Error relativo, residuo, norma, espectro y coeficientes de Picard | Completo |
| Discrepancia sin $x_{\mathrm{true}}$ comparada con oráculo | Completo |
| Contraste fórmula espectral vs. sistema aumentado (QR) y ecuaciones normales | Completo |
| Robustez con 50 semillas y sensibilidad a $\delta$ | Completo |
| Variación del suavizado $s$ y barrido de $\eta$ hasta $10^{-10}$ | Completo |
| Ejemplo práctico de imagen (separable, Kronecker) con pruebas | Completo |
| Ejemplo práctico de audio (circulante, Fourier) con pruebas | Completo |
| Correcciones de la revisión matemática de la entrega parcial | Completo (ver `DECLARACION_IA.md`) |
| Ampliación: penalización $\lambda\lVert Lx\rVert^2$ con $L_1$ y $L_2$, GSVD, 50 semillas, dos señales | Completo |
| Estimación de $\delta$ desde los datos y refinamiento de $\lambda$ por bisección | Completo |
| Preguntas de defensa con respuestas | Completo (desarrollo §11) |
| Curva L | Solo gráfico (no se usa para elegir) |
| Guion y diapositivas del video | Completo (`docs/guion_video.md`, `docs/diapositivas/`) |
| Grabación y edición del video | Completo |

