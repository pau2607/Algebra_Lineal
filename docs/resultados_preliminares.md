# Resultados preliminares e interpretación

Todos los números de este documento salen de `python scripts/ejecutar_todo.py` (semilla 2026, $n=128$, $s=0.03$, $\tau=1.1$, malla de 61 valores $\lambda/\sigma_1^2\in[10^{-10},1]$). Las tablas completas están en `resultados/tablas/` y las figuras en `resultados/figuras/`.

## 1. Comprobaciones realizadas

| Comprobación | Cómo | Resultado |
|---|---|---|
| Ejemplo diagonal coincide con el cálculo a mano | `scripts/01_ejemplo_diagonal.py` (aserciones) y `tests/` | $x^\dagger=(1,11)$, $x_1=(1,0)$, $x_\lambda=(0.99990,\,0.10891)$; errores 10, 1 y 0.8911 ✔ |
| Filas de $A$ suman 1 y $A\neq A^T$ | `test_matriz_conserva_constantes_y_no_es_simetrica` | ✔ |
| Ruido relativo exactamente $\eta$ | `test_ruido_tiene_nivel_relativo_exacto` | ✔ (error relativo $<10^{-12}$) |
| Fórmula espectral = sistema aumentado (QR) = ecuaciones normales | `test_tikhonov_tres_formas_coinciden` y script 02 | diferencia relativa máx. $3\times10^{-10}$ en el problema real ✔ |
| $x_\lambda$ anula el gradiente de la función objetivo | `test_tikhonov_es_punto_estacionario` | ✔ |
| Cota $\lVert x_\lambda(b+e)-x_\lambda(b)\rVert\le\lVert e\rVert/2\sqrt\lambda$ | `test_cota_de_estabilidad` (60 perturbaciones) | ✔ |
| Límites $\lambda\to0^+$ y $\lambda\to\infty$ | `test_limites_de_lambda` | ✔ |
| Monotonía del residuo (base de la discrepancia) | `test_residuo_*` | ✔ |
| La discrepancia no usa $x_{\mathrm{true}}$ | `test_discrepancia_no_recibe_la_verdad` | ✔ |

Total: 23 pruebas, todas aprobadas (`python -m pytest -q`). Los ejemplos prácticos de la sección 7 agregan otras 22 (`tests/test_aplicaciones.py`): 45 en total.

## 2. Espectro del operador

- $\sigma_1\approx1.001$ y los $\sigma_i$ decaen de forma aproximadamente gaussiana hasta $\sim10^{-14}$ en $i\approx90$; desde $i\approx95$ hay una meseta en $\sim10^{-16}$ que es redondeo de la SVD.
- El $\kappa_2(A)$ calculado es $2.2\times10^{16}$: $A$ es numéricamente singular. `scipy.linalg.pinv` (umbral $2.84\times10^{-14}$) y `numpy.linalg.matrix_rank` dan rango numérico 89.
- En el gráfico de Picard (fig. 03), los $|u_i^Tb|$ siguen a $\sigma_i$ hasta que chocan con el nivel del ruido, $\approx\lVert e\rVert/\sqrt n$ ($4\times10^{-5}$, $4\times10^{-4}$ y $4\times10^{-3}$), y a partir de ahí quedan planos. Solo $\approx44$ valores singulares superan el nivel de ruido para $\eta=10^{-3}$ (tabla 04): esas son las direcciones que los datos realmente informan.

## 3. Comparación de métodos (semilla 2026)

Fragmento de `resultados/tablas/02_resumen_metodos.md` (error = error relativo):

| $\eta$ | Pinv numérica (rango 89) | TSVD discrepancia | TSVD oráculo | Tikhonov discrepancia | Tikhonov oráculo |
|---|---|---|---|---|---|
| $10^{-4}$ | $8.1\times10^{7}$ | 0.0939 ($k=38$) | 0.0931 ($k=40$) | 0.0941 ($\lambda/\sigma_1^2=6.8\times10^{-6}$) | 0.0914 ($3.2\times10^{-7}$) |
| $10^{-3}$ | $8.1\times10^{8}$ | 0.1139 ($k=31$) | 0.1029 ($k=38$) | 0.1087 ($3.2\times10^{-4}$) | 0.1008 ($2.2\times10^{-5}$) |
| $10^{-2}$ | $8.1\times10^{9}$ | 0.1218 ($k=23$) | 0.1185 ($k=27$) | 0.1194 ($4.6\times10^{-3}$) | 0.1191 ($3.2\times10^{-3}$) |

### Interpretación

1. La pseudoinversa es inservible aunque ajuste los datos. Su residuo está entre $0.54$ y $0.69\,\lVert e\rVert$, es decir, *más pequeño que el propio ruido*: está reproduciendo ruido. Su norma es $4\times10^{8}$ (la verdad tiene $\lVert x_{\mathrm{true}}\rVert\approx4.96$). El error crece exactamente $\times10$ cuando $\eta$ crece $\times10$ porque domina el término lineal $A^\dagger e$ y usamos la misma dirección $z$ para los tres niveles.
2. Cuidado con el residuo de la pseudoinversa completa (128 $\sigma_i$). En teoría debería ser casi cero, pero el calculado vale $\approx\lVert e\rVert$. Con $\lVert x\rVert\approx3\times10^{12}$, el redondeo al calcular $Ax$ es del orden de $\varepsilon\,\lVert A\rVert\,\lVert x\rVert\approx7\times10^{-4}$, comparable con $\lVert e\rVert$. Ese residuo no es confiable.
3. La regularización reduce el error de $10^{8}$–$10^{10}$ a 9–12 %. Las reconstrucciones (fig. 06) recuperan bien la zona suave y la posición del pulso, con oscilaciones cerca de sus bordes.
4. La discrepancia sobrerregulariza un poco. Elige $\lambda$ mayores (y $k$ menores) que el oráculo: el residuo del oráculo está en $0.77$–$0.94\,\lVert e\rVert$, por debajo del umbral $1.1\,\lVert e\rVert$. El costo en error es pequeño (ver sección 4).
5. Tikhonov y TSVD dan resultados parecidos. Coinciden en lo esencial porque los $\sigma_i$ caen muy rápido: el filtro suave y el corte abrupto descartan prácticamente las mismas direcciones (fig. 08).

## 4. Robustez (50 semillas) y sensibilidad a $\delta$

Cociente error(discrepancia) / error(oráculo), de `resultados/tablas/03_robustez_semillas.md`:

| $\eta$ | $\delta/\lVert e\rVert$ | Tikhonov: mediana (máx.) | TSVD: mediana (máx.) |
|---|---|---|---|
| $10^{-4}$ | 1 | 1.016 (1.07) | 1.013 (1.05) |
| $10^{-3}$ | 1 | 1.081 (1.15) | 1.135 (1.20) |
| $10^{-2}$ | 1 | 1.005 (1.04) | 1.041 (1.12) |
| todos | 2 | 1.06–1.12 (≤ 1.20) | 1.09–1.27 (≤ 1.35) |
| todos | 0.5 | no viable en 50/50 semillas | $\sim10^{8}$–$10^{10}$ (no viable en 13/50) |

- Con $\delta$ correcto, la mediana del error de la regla queda entre 0.5 % y 14 % por encima del oráculo, sin consultar la verdad.
- Subestimar el ruido es mucho peor que sobrestimarlo. Con $\delta=0.5\lVert e\rVert$ el umbral $0.55\lVert e\rVert$ está por debajo del residuo mínimo de Tikhonov en la malla ($\approx0.70\lVert e\rVert$ en $\lambda/\sigma_1^2=10^{-10}$), así que la regla no encuentra $\lambda$ y lo informa. En la TSVD el umbral solo se alcanza con $k$ muy grandes, que reproducen ruido: error de $10^{8}$ veces el oráculo. Con $\delta=2\lVert e\rVert$ se sobrerregulariza, pero la mediana del error sube solo entre 6 % y 27 % sobre el oráculo (máximo 35 %).

## 5. ¿Por qué el error no baja de ≈9 % aunque el ruido sea muy pequeño?

- En los tres niveles de ruido, entre el 96 % y el 99 % del error cuadrático de Tikhonov (con $\lambda$ por discrepancia) se concentra en $t\in[0.55,0.8]$, alrededor del pulso, que es solo 32 de los 128 puntos.
- Separando $x_\lambda-x_{\mathrm{true}}$ en sesgo y ruido propagado (sección 7.2 del desarrollo), para $\eta=10^{-3}$ el sesgo aporta $\approx0.108$ del error relativo y el ruido solo $\approx0.007$.
- Sin ruido ($e=0$, tabla 03), el mejor error posible es 5.4 % (TSVD, $k=91$) y 7.8 % (Tikhonov en el extremo inferior de la malla). Las discontinuidades del pulso necesitan componentes $v_i$ de alta frecuencia cuyos $\sigma_i$ están por debajo de la precisión de máquina: esa información no está en los datos de forma utilizable. Es el límite que el artículo anticipa: la regularización no recupera lo que la observación destruyó.

## 6. Exploración: más suavizado, menos información

| $s$ | rango pinv | $\sigma_i$ sobre el ruido | error Tikhonov (discrepancia) | error oráculo |
|---|---|---|---|---|
| 0.02 | 128 | 65 | 0.085 | 0.079 |
| 0.03 | 89 | 44 | 0.109 | 0.101 |
| 0.05 | 56 | 27 | 0.138 | 0.126 |

Al aumentar $s$ los $\sigma_i$ caen más rápido, menos direcciones superan el nivel del ruido y el error mínimo alcanzable sube (fig. 09).

## 7. Ejemplos prácticos: imagen y sonido (después de la entrega parcial)

Generados con `python scripts/05_imagen.py` y `python scripts/06_audio.py`, con la misma semilla (2026), $\tau=1.1$, $\delta=\lVert E\rVert_F$ y la misma malla de $\lambda/\sigma_1^2$. El suavizado tiene el ancho del experimento principal: 3.81 píxeles o muestras. El método está en la sección 9 del desarrollo. Las pruebas de `tests/test_aplicaciones.py` (22, todas aprobadas) verifican las fórmulas contra las matrices grandes formadas explícitamente en casos pequeños.

### 7.1 Imagen de $256\times256$ (65 536 incógnitas)

Tablas completas en `resultados/tablas/05_imagen.csv`. PSNR = $10\log_{10}(1/\text{ECM})$ en dB; más alto es mejor.

| Imagen, $\eta$ | Observada | Pinv numérica (rango 18 817) | SVD truncada · discrepancia | Tikhonov · discrepancia | Tikhonov · oráculo |
|---|---|---|---|---|---|
| gato, $10^{-3}$ | 0.117 · 22.3 dB | $5.1\times10^{6}$ | 0.074 · 26.3 dB ($k=2876$) | **0.070 · 26.8 dB** | 0.067 · 27.2 dB |
| gato, $10^{-2}$ | 0.118 · 22.3 dB | $5.1\times10^{7}$ | 0.092 · 24.4 dB ($k=1586$) | **0.081 · 25.6 dB** | 0.078 · 25.8 dB |
| texto, $10^{-3}$ | 0.152 · 16.6 dB | $5.1\times10^{6}$ | 0.105 · 19.9 dB ($k=3338$) | **0.102 · 20.1 dB** | 0.093 · 20.9 dB |
| texto, $10^{-2}$ | 0.152 · 16.6 dB | $5.1\times10^{7}$ | 0.141 · 17.3 dB ($k=1708$) | **0.123 · 18.4 dB** | 0.114 · 19.2 dB |

(error relativo · PSNR)

Interpretación:

1. **Se ve lo mismo que en la señal 1D, ahora a simple vista.** La pseudoinversa da ruido puro (figuras 10 y 11), con error de $10^{6}$ a $10^{7}$ (pinv numérica) y de $10^{28}$ a $10^{29}$ (sin umbral). Tikhonov con discrepancia devuelve una imagen nítida sin mirar la original: los ojos y el pelaje del gato reaparecen, y el texto "APLICADA" vuelve a ser legible (figura 13, ampliación).
2. **El error de la imagen observada casi no depende de $\eta$** (0.117 frente a 0.118). Lo que la degrada es el desenfoque, no el ruido. El ruido decide, en cambio, *cuánto* se puede deshacer el desenfoque: con $\eta=10^{-3}$ se conservan $k\approx2900$ direcciones y con $\eta=10^{-2}$ solo $\approx1600$. Es el gráfico de Picard de la figura 12: los $|C_{ij}|$ siguen a los $S_{ij}$ hasta chocar con el nivel del ruido.
3. **Tikhonov supera a la SVD truncada en las cuatro combinaciones**, y la ventaja crece con el ruido (texto, $\eta=10^{-2}$: 18.4 frente a 17.3 dB). En la imagen el corte abrupto de la TSVD deja ondulaciones visibles alrededor de las letras (efecto de Gibbs). En 1D los dos métodos empataban porque los $\sigma_i$ caen muy rápido. Aquí hay muchos más $S_{ij}$ de tamaño intermedio (cada nivel de $S$ se alcanza con muchas parejas $(i,j)$; curva negra de la figura 12), y nuestra explicación es que por eso importa más cómo se atenúan esas direcciones intermedias.
4. **La discrepancia queda cerca del oráculo**: la pérdida es de 0.2 a 0.8 dB, aunque elige un $\lambda$ entre 3 y 15 veces mayor que el oráculo (sobrerregulariza, igual que en 1D). En la SVD truncada, $k$ se busca en una malla logarítmica de 150 valores entre 1 y 65 536, no en todos los $k$.
5. **El residuo de la pinv numérica es $0.844\,\delta$ en los cuatro casos.** No es casualidad: esa reconstrucción deja fuera las mismas $65\,536-18\,817$ direcciones, el residuo es la parte del ruido en esas direcciones y usamos la misma dirección de ruido $Z$. Por eso el cociente es igual en los cuatro casos (medido: 0.8439) y está muy cerca del valor esperado para ruido isotrópico, $\sqrt{1-18\,817/65\,536}\approx0.8443$. El "residuo" de la pseudoinversa completa ($10^{15}\delta$) es redondeo, por el mismo motivo que la observación 2 de la sección 3: $\lVert X\rVert\sim10^{31}$.
6. **El texto es más difícil que la foto** (error de 0.10 a 0.12 frente a 0.07 a 0.08), igual que el pulso en 1D: los bordes nítidos de las letras necesitan direcciones de alta frecuencia con $S_{ij}$ diminutos.

### 7.2 Audio de 2 s (16 384 muestras a 8000 Hz)

Tabla en `resultados/tablas/06_audio.csv`; audios en `resultados/audio/`. SNR $=10\log_{10}(\lVert x\rVert^2/\lVert x_{\text{rec}}-x\rVert^2)$.

| $\eta$ | Observado | Pseudoinversa completa | Tikhonov · discrepancia | Tikhonov · oráculo |
|---|---|---|---|---|
| $10^{-3}$ | 0.592 · 4.6 dB | $1.9\times10^{13}$ | **0.266 · 11.5 dB** | 0.246 · 12.2 dB |
| $10^{-2}$ | 0.592 · 4.6 dB | $1.9\times10^{14}$ | **0.329 · 9.7 dB** | 0.296 · 10.6 dB |

(error relativo · SNR)

Interpretación (figura 14):

1. **El suavizado apaga los agudos.** El sonido observado conserva el tono de cada nota (los graves), pero pierde el "brillo" de los armónicos. En la ventana de 20 ms de la figura, la amplitud cae de 0.72 a 0.21 y desaparecen los picos (en todo el audio, el pico baja de 0.80 a 0.51).
2. **Tikhonov devuelve los armónicos hasta cerca de 1000 Hz** y la forma de onda vuelve a parecerse a la original (panel izquierdo). Por encima de ~1200 Hz el factor de filtro corta: ahí $\hat h_k$ está por debajo del ruido y esa información ya no está en los datos (panel central). Es la versión audible del límite de la sección 5: *la regularización no recupera lo que la observación destruyó*. Por eso el error no baja de 25 a 30 %.
3. **La pseudoinversa es inaudible como música:** divide entre $\hat h_k\sim10^{-17}$ (redondeo), y el resultado es ruido blanco de amplitud $10^{13}$. El archivo `..._pseudoinversa_..._VOLUMEN_REDUCIDO.wav` está reescalado para no dañar los parlantes.
4. **El ruido decide cuánto se recupera.** Con $\eta=10^{-3}$ la discrepancia permite un $\lambda$ 46 veces menor que con $\eta=10^{-2}$ y la SNR sube 1.8 dB.

## 8. Limitaciones de estos resultados

- $\delta=\lVert e\rVert$ es exacto porque generamos el ruido; en datos reales habría que estimarlo, y la sección 4 muestra que subestimarlo es peligroso.
- Los ejemplos de imagen y audio usan una sola semilla de ruido; no repetimos para ellos el estudio de robustez con 50 semillas.
- En la imagen, el desenfoque es gaussiano y separable con bordes renormalizados; otros desenfoques reales (movimiento en diagonal, desenfoque circular de una lente) no son separables y requeriría otra forma de calcular la SVD. En el audio, la convolución es circular (sección 9.2 del desarrollo).
- El ruido es gaussiano blanco; con ruido correlacionado o errores en $A$ las conclusiones pueden cambiar.
- La malla de $\lambda$ es finita (6 puntos por década): el $\lambda$ elegido tiene resolución de $10^{1/6}\approx1.47$. En el caso sin ruido el óptimo de Tikhonov cae en el borde de la malla.
- Un buen resultado con una señal sintética no garantiza el mismo desempeño con datos reales.
