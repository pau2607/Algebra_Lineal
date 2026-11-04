# Organización del trabajo

## 1. Integrantes y responsabilidades

| Integrante | Responsabilidad principal | Archivos a su cargo | Bloque del video |
|---|---|---|---|
| Paula Romero | Planteamiento, construcción del problema, SVD, pseudoinversa y TSVD | `reconstruccion/datos.py`, `scripts/01_ejemplo_diagonal.py`, parte de pseudoinversa y TSVD de `reconstruccion/metodos.py`, §1 a §5 del desarrollo | 0:00-3:30. Problema, ejemplo diagonal, SVD y pseudoinversa, espectro de A y pinv numérica |
| Carlos Ledo | Tikhonov, estabilidad, sesgo y robustez | parte de Tikhonov de `reconstruccion/metodos.py`, `tests/`, `scripts/03_robustez.py`, `scripts/04_variacion_s.py`, §6 y §7 del desarrollo | 3:30-6:35. Tikhonov y factores de filtro, estabilidad y sesgo, robustez y variación de s |
| Matias Desmonteix | Selección del parámetro, experimento principal e interpretación | `reconstruccion/seleccion.py`, `scripts/02_experimento_principal.py`, `docs/resultados_preliminares.md`, guion, §8 del desarrollo | 6:35-10:00. Curvas y selección del parámetro, reconstrucciones, interpretación, conclusiones y créditos |

Cada uno habla alrededor de 3 minutos de contenido técnico. Todos tenemos que poder explicar el razonamiento completo, porque la ronda de preguntas es individual.

## 2. Registro de contribuciones (entrega parcial)

| Integrante | Aporte concreto en esta entrega | Evidencia (archivos y commits) |
|---|---|---|
| Paula Romero | Planteé el problema y construí los datos de prueba: implementé la matriz de suavizado gaussiano con $s=0.03$, la señal (zona suave más pulso) y el ruido con semilla registrada para los tres valores de $\eta$. Resolví a mano el ejemplo diagonal y lo reproduje en código. Implementé la pseudoinversa numérica, con su rango y su umbral, y la SVD truncada para todos los $k$. Analicé el espectro de $A$ y los coeficientes de Picard. Redacté las secciones 1 a 5 del desarrollo. | `reconstruccion/datos.py`, `scripts/01_ejemplo_diagonal.py`, `reconstruccion/metodos.py` (pseudoinversa y TSVD), desarrollo §1 a §5. Commits: completar |
| Carlos Ledo | Implementé Tikhonov con la fórmula espectral y lo contrasté con el sistema aumentado resuelto por QR. Escribí las pruebas de `tests/` para verificar los métodos. Hice el estudio de robustez con 50 semillas y la sensibilidad a $\delta$, y la variación del suavizado $s$. Redacté las secciones 6 y 7 del desarrollo (Tikhonov, estabilidad y sesgo). | `reconstruccion/metodos.py` (Tikhonov), `tests/`, `scripts/03_robustez.py`, `scripts/04_variacion_s.py`, desarrollo §6 y §7. Commits: completar |
| Matias Desmonteix | Implementé la selección del parámetro: principio de discrepancia sin usar $x_{\mathrm{true}}$, comparado con el oráculo, y el gráfico de la curva L. Armé el experimento principal que genera las curvas de error, residuo y norma, las reconstrucciones y las tablas. Redacté los resultados preliminares, la sección 8 del desarrollo y el guion del video. | `reconstruccion/seleccion.py`, `scripts/02_experimento_principal.py`, `docs/resultados_preliminares.md`, guion, desarrollo §8. Commits: completar |

## 3. Estado actual

| Elemento | Estado |
|---|---|
| Ejemplo diagonal reproducido y verificado | Completo |
| Señal, matriz y ruido ($s=0.03$, tres $\eta$, semilla registrada) | Completo |
| Pseudoinversa numérica (con rango y umbral), TSVD para todos los $k$, Tikhonov en malla logarítmica | Completo |
| Error relativo, residuo, norma, espectro y coeficientes de Picard | Completo |
| Discrepancia sin $x_{\mathrm{true}}$ comparada con oráculo | Completo |
| Contraste fórmula espectral vs. sistema aumentado con QR | Completo |
| Robustez con 50 semillas y sensibilidad a $\delta$ | Preliminar |
| Variación del suavizado $s$ | Exploratorio |
| Curva L | Parcial, solo gráfico, sin selección por esquina |
| Ampliación (penalización $\lVert Lx\rVert^2$ o GCV) | Pendiente |
| Video | Pendiente, guion preliminar listo |

## 4. Actividades pendientes

1. Incorporar la retroalimentación del docente a esta entrega (todos).
2. Ampliación elegida por el grupo (una de las dos):
   - penalización $\lambda\lVert Lx\rVert^2$ con $L$ de primeras diferencias, comprobando $\ker A\cap\ker L=\{0\}$ (la lidera Carlos);
   - selección por curva L (esquina por máxima curvatura) o validación cruzada generalizada, comparada con la discrepancia (la lidera Matias).
3. Refinar la malla de $\lambda$ cerca del valor elegido y ampliarla hacia abajo para el caso sin ruido (Carlos).
4. Estudiar una estimación de $\delta$ a partir de los datos (por ejemplo, con los coeficientes de Picard planos) en lugar de usar $\lVert e\rVert$ exacto (Paula).
5. Preparar las diapositivas de su propio bloque, ensayar con cronómetro, grabar y editar el video (9 a 11 min). La edición la hace Matias.
6. Completar el README final (Paula), el registro de contribuciones y la declaración de IA (todos), y poner la etiqueta `entrega-final` (Paula).
7. Preparación individual de la ronda de preguntas: cada uno repasa las seis preguntas de defensa (§9 del desarrollo) y un ejemplo de variación (cambiar $s$, $\eta$, $\tau$ o $\delta$).

## 5. Cronograma hasta la entrega final

Las fechas son referenciales y hay que ajustarlas a las publicadas en Blackboard.

| Semana | Fechas | Actividad | Responsables |
|---|---|---|---|
| 1 | 9-15 nov 2026 | Retroalimentación de la parcial; elegir y empezar la ampliación | Todos, y Carlos y Matias con la ampliación |
| 2 | 16-22 nov 2026 | Terminar la ampliación con sus pruebas; estimación de $\delta$; actualizar resultados y documentos | Carlos, Matias y Paula |
| 3 | 23-29 nov 2026 | Diapositivas, guion definitivo, ensayo cronometrado y grabación | Todos |
| 4 | 30 nov-6 dic 2026 | Edición del video (Matias); revisión de reproducibilidad desde cero (Carlos); README final y etiqueta `entrega-final` (Paula) | Todos |
| Final | Según Blackboard | Ronda de preguntas individual | Todos |
