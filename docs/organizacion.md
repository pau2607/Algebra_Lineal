# Organización del trabajo

> **Completar antes de entregar:** nombres reales, fechas oficiales de Blackboard y la descripción de lo que cada integrante hizo realmente. La rúbrica valora el contenido del aporte y la capacidad de explicarlo, no el número de commits.

## 1. Integrantes y responsabilidades

| Integrante | Responsabilidad principal | Archivos a su cargo | Bloque del video |
|---|---|---|---|
| Integrante A | Planteamiento, construcción del problema y ejemplo diagonal | `reconstruccion/datos.py`, `scripts/01_ejemplo_diagonal.py`, §1–2 y §4.1 del desarrollo | Problema y ejemplo diagonal |
| Integrante B | Métodos: pseudoinversa, SVD truncada, Tikhonov y su verificación | `reconstruccion/metodos.py`, `tests/`, §3, §5 y §6 del desarrollo | SVD, pseudoinversa y Tikhonov |
| Integrante C | Estabilidad, sesgo y selección del parámetro | `reconstruccion/seleccion.py`, `scripts/03_robustez.py`, §7–8 del desarrollo | Estabilidad, Picard y pinv numérica |
| Integrante D | Experimento principal, figuras, tablas e interpretación | `scripts/02_experimento_principal.py`, `scripts/04_variacion_s.py`, `docs/resultados_preliminares.md`, guion | Curvas, reconstrucciones y conclusiones |

Todos los integrantes deben poder explicar el razonamiento completo (ronda de preguntas individual).

## 2. Registro de contribuciones (entrega parcial)

| Integrante | Aporte concreto en esta entrega | Evidencia (archivos / commits) |
|---|---|---|
| Integrante A | *(describir: p. ej. "implementó la matriz de suavizado y comprobó que conserva constantes; redactó §1 y el ejemplo diagonal")* | |
| Integrante B | | |
| Integrante C | | |
| Integrante D | | |

## 3. Estado actual

| Elemento | Estado |
|---|---|
| Ejemplo diagonal reproducido y verificado | ✔ completo |
| Señal, matriz y ruido ($s=0.03$, tres $\eta$, semilla registrada) | ✔ completo |
| Pseudoinversa numérica (con rango y umbral), TSVD para todos los $k$, Tikhonov en malla logarítmica | ✔ completo |
| Error relativo, residuo, norma, espectro y coeficientes de Picard | ✔ completo |
| Discrepancia sin $x_{\mathrm{true}}$ comparada con oráculo | ✔ completo |
| Contraste fórmula espectral vs. sistema aumentado con QR | ✔ completo |
| Robustez con 50 semillas y sensibilidad a $\delta$ | ✔ preliminar |
| Variación del suavizado $s$ | ✔ exploratorio |
| Curva L | ◐ solo gráfico, sin selección por esquina |
| Ampliación (penalización $\lVert Lx\rVert^2$ o GCV) | ☐ pendiente |
| Video | ☐ pendiente (guion preliminar listo) |

## 4. Actividades pendientes

1. Incorporar la retroalimentación del docente a esta entrega.
2. Ampliación elegida por el grupo (una de):
   - penalización $\lambda\lVert Lx\rVert^2$ con $L$ de primeras diferencias, comprobando $\ker A\cap\ker L=\{0\}$;
   - selección por curva L (esquina por máxima curvatura) o validación cruzada generalizada, comparada con la discrepancia.
3. Refinar la malla de $\lambda$ cerca del valor elegido y ampliarla hacia abajo para el caso sin ruido.
4. Estudiar una estimación de $\delta$ a partir de los datos (p. ej. con los coeficientes de Picard planos) en lugar de usar $\lVert e\rVert$ exacto.
5. Preparar diapositivas, ensayar con cronómetro, grabar y editar el video (9–11 min).
6. Completar README final, registro de contribuciones, declaración de IA y etiqueta `entrega-final`.
7. Preparación individual de la ronda de preguntas: cada integrante repasa las seis preguntas de defensa (§9 del desarrollo) y un ejemplo de "variación" (cambiar $s$, $\eta$, $\tau$ o $\delta$).

## 5. Cronograma hasta la entrega final

Las fechas son referenciales; ajustar a las publicadas en Blackboard.

| Semana | Fechas | Actividad | Responsables |
|---|---|---|---|
| 1 | 9–15 nov 2026 | Retroalimentación de la parcial; elegir y empezar la ampliación | Todos / B y C |
| 2 | 16–22 nov 2026 | Terminar la ampliación con sus pruebas; actualizar resultados y documentos | B, C, D |
| 3 | 23–29 nov 2026 | Diapositivas, guion definitivo, ensayo cronometrado y grabación | Todos |
| 4 | 30 nov–6 dic 2026 | Edición del video; revisión de reproducibilidad desde cero; etiqueta `entrega-final` | A y D (revisión: todos) |
| — | Según Blackboard | Ronda de preguntas individual | Todos |
