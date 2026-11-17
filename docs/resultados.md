# Resultados e interpretación

Todos los números de este documento salen de `python scripts/ejecutar_todo.py`: semilla 2026, $n=128$, $s=0.03$, $\tau=1.1$ y una malla de 61 valores $\lambda/\sigma_1^2\in[10^{-10},1]$, salvo donde se indica otra cosa. Las tablas completas están en `resultados/tablas/` y las figuras en `resultados/figuras/`. Los diagnósticos de las secciones 2, 3 y 5 vienen de `scripts/07_diagnosticos.py` (tablas `07_*`); los ejemplos prácticos, de `scripts/05_imagen.py` y `scripts/06_audio.py`; la ampliación, de `scripts/08_penalizacion_diferencias.py`; la estimación de $\delta$, de `scripts/09_estimacion_delta.py`.

## 1. Comprobaciones realizadas

| Comprobación | Cómo | Resultado |
|---|---|---|
| Ejemplo diagonal coincide con el cálculo a mano | `scripts/01_ejemplo_diagonal.py` (aserciones) y `tests/` | $x^\dagger=(1,11)$, $x_1=(1,0)$, $x_\lambda=(0.99990,\,0.10891)$; errores 10, 1 y 0.8911 ✔ |
| Filas de $A$ suman 1 y $A\neq A^T$ | `test_matriz_conserva_constantes_y_no_es_simetrica` | ✔ |
| Ruido relativo exactamente $\eta$ | `test_ruido_tiene_nivel_relativo_exacto` | ✔ (error relativo $<10^{-12}$) |
| Fórmula espectral = sistema aumentado (QR) | script 02, problema real, 13 valores de $\lambda$ | diferencia relativa máx. $3\times10^{-10}$ ✔ |
| Fórmula espectral = ecuaciones normales regularizadas | `test_tikhonov_tres_formas_coinciden` (matriz de prueba) y script 07 (problema real, 61 valores de $\lambda$) | diferencia relativa máx. $4.5\times10^{-7}$ en el problema real ✔ |
| $x_\lambda$ anula el gradiente de la función objetivo | `test_tikhonov_es_punto_estacionario` | ✔ |
| Cota $\lVert x_\lambda(b+e)-x_\lambda(b)\rVert\le\lVert e\rVert/2\sqrt\lambda$ | `test_cota_de_estabilidad` (60 perturbaciones) | ✔ |
| Límites $\lambda\to0^+$ y $\lambda\to\infty$ | `test_limites_de_lambda` | ✔ |
| Monotonía del residuo en $\lambda$ y en $k$ | `test_residuo_*` (matriz de prueba bien escalada) | ✔; con nuestra $A$ la del residuo en $k$ falla numéricamente desde $k=94$ (ver §2) |
| La discrepancia no usa $x_{\mathrm{true}}$ | `test_discrepancia_no_recibe_la_verdad` | ✔ |

| GSVD reconstruye $A$ y $L$; fórmula generalizada = sistema aumentado | `tests/test_general.py` | error $10^{-14}$; diferencia relativa $\le7\times10^{-8}$ ✔ |
| Kronecker y circulantes coinciden con la matriz grande formada explícitamente | `tests/test_aplicaciones.py` (casos pequeños) | ✔ |
| $\hat\delta$ no recibe $e$ y aproxima $\lVert e\rVert$ | `test_estimar_delta_no_usa_e_y_aproxima_su_norma` | ✔ |

Total: 61 pruebas, todas aprobadas (`python -m pytest -q`): 25 del experimento principal, 22 de imagen y audio, 14 de la ampliación.

## 2. Espectro del operador

- Forma del espectro. $\sigma_1\approx1.001$. Los $\sigma_i$ decaen de forma aproximadamente gaussiana: $\log\sigma_i$ es casi cuadrático en $i$. Llegan a $\sim4\times10^{-14}$ en $i=89$.
- Meseta de redondeo. Desde $i=98$, los $\sigma_i$ calculados forman una meseta en $\approx7\times10^{-17}$. Es nivel de redondeo: en aritmética exacta $A$ es invertible (ver desarrollo §3), pero sus $\sigma_i$ verdaderos más pequeños son muchísimo menores que $\varepsilon_{\text{máq}}\sigma_1$.
- El $\kappa_2(A)$ calculado no es el verdadero. Vale $\approx2\times10^{16}$, pero solo dice que $A$ es numéricamente singular ($\kappa_2\gtrsim1/\varepsilon_{\text{máq}}$). Por eso la tabla 04 lo marca como no fiable para $s=0.03$ y $s=0.05$.
- Rango numérico. `scipy.linalg.pinv` (umbral $2.84\times10^{-14}$) y `numpy.linalg.matrix_rank` dan 89. Formar $A^TA$ baja el rango numérico a 63.
- Residuo de la TSVD. El residuo calculado deja de ser monótono en $k$ desde $k=94$, cuando entran direcciones con $\sigma_i$ al nivel del redondeo.
- Gráfico de Picard (fig. 03). Los $|u_i^Tb|$ decaen junto con los $|u_i^TAx_{\mathrm{true}}|=\sigma_i|v_i^Tx_{\mathrm{true}}|$, que bajan al menos tan rápido como $\sigma_i$. Cuando chocan con el nivel del ruido, $\approx\lVert e\rVert/\sqrt n$ ($4\times10^{-5}$, $4\times10^{-4}$ y $4\times10^{-3}$), quedan planos.
- Direcciones informadas para $\eta=10^{-3}$ (tabla `07_diagnosticos`):

  | Recuento | Valor |
  |---|---|
  | $\sigma_i$ que superan el nivel del ruido | 44 |
  | Coeficientes limpios $|u_i^TAx_{\mathrm{true}}|$ que lo superan | 33 |
  | Índice del último de esos coeficientes | 37 |

  Lo que cuenta son los coeficientes, no los $\sigma_i$: los datos informan realmente unas 35–40 direcciones. Eso coincide con los $k$ elegidos por la TSVD (31 por discrepancia y 38 por oráculo).

## 3. Comparación de métodos (semilla 2026)

Fragmento de `resultados/tablas/02_resumen_metodos.md` (error = error relativo):

| $\eta$ | Pinv numérica (rango 89) | TSVD discrepancia | TSVD oráculo | Tikhonov discrepancia | Tikhonov oráculo |
|---|---|---|---|---|---|
| $10^{-4}$ | $8.1\times10^{7}$ | 0.0939 ($k=38$) | 0.0931 ($k=40$) | 0.0941 ($\lambda/\sigma_1^2=6.8\times10^{-6}$) | 0.0914 ($3.2\times10^{-7}$) |
| $10^{-3}$ | $8.1\times10^{8}$ | 0.1139 ($k=31$) | 0.1029 ($k=38$) | 0.1087 ($3.2\times10^{-4}$) | 0.1008 ($2.2\times10^{-5}$) |
| $10^{-2}$ | $8.1\times10^{9}$ | 0.1218 ($k=23$) | 0.1185 ($k=27$) | 0.1194 ($4.6\times10^{-3}$) | 0.1191 ($3.2\times10^{-3}$) |

### Interpretación

1. La pseudoinversa es inservible aunque ajuste los datos.
   - Calculada como TSVD con $k=89$, su residuo es $0.534\,\lVert e\rVert$ en los tres niveles. Esto concuerda con el valor esperado $\sqrt{39/128}\approx0.55$: ajusta la parte del ruido que cae en las 89 direcciones conservadas, es decir, *reproduce ruido*.
   - Al formar `pinv(A)` explícitamente, el residuo sube a $0.688\,\lVert e\rVert$ para $\eta=10^{-4}$, por redondeo de la matriz explícita, que no depende de $\eta$.
   - Su norma es $4\times10^{8}$, $4\times10^{9}$ y $4\times10^{10}$, mientras que $\lVert x_{\mathrm{true}}\rVert\approx4.96$.
   - El error crece prácticamente $\times10$ cuando $\eta$ crece $\times10$. Domina el término lineal $A^\dagger e$, y usamos la misma dirección $z$ para los tres niveles.
2. Cuidado con el residuo de la pseudoinversa completa (128 $\sigma_i$).
   - En aritmética exacta sería casi cero; el calculado vale $\approx\lVert e\rVert$.
   - La causa: $\lVert x\rVert\approx3\times10^{12}$, y cualquier error relativo de orden $\varepsilon_{\text{máq}}$ en los factores de la SVD o al formar $x$ produce un residuo de orden $\varepsilon_{\text{máq}}\lVert A\rVert\lVert x\rVert\approx1.5\,\lVert e\rVert$ en los tres niveles.
   - Ese residuo calculado no es confiable.
3. La regularización reduce el error de $10^{8}$–$10^{10}$ a 9–12 %. Las reconstrucciones (fig. 06) recuperan bien la zona suave y la posición del pulso, con oscilaciones cerca de sus bordes.
4. La discrepancia sobrerregulariza un poco. Elige $\lambda$ mayores (y $k$ menores) que el oráculo. El residuo del oráculo está en $0.77$–$0.94\,\lVert e\rVert$, por debajo del umbral $1.1\,\lVert e\rVert$. El costo en error es pequeño (sección 4).
5. Tikhonov y TSVD dan resultados parecidos. Los $\sigma_i$ caen muy rápido, así que el filtro suave y el corte abrupto descartan prácticamente las mismas direcciones (fig. 08).

## 4. Robustez (50 semillas) y sensibilidad a $\delta$

Cociente error(discrepancia) / error(oráculo), de `resultados/tablas/03_robustez_semillas.md`:

| $\eta$ | $\delta/\lVert e\rVert$ | Tikhonov: mediana (máx.) | TSVD: mediana (máx.) |
|---|---|---|---|
| $10^{-4}$ | 1 | 1.016 (1.07) | 1.013 (1.05) |
| $10^{-3}$ | 1 | 1.081 (1.15) | 1.135 (1.20) |
| $10^{-2}$ | 1 | 1.005 (1.04) | 1.041 (1.12) |
| todos | 2 | 1.06–1.12 (≤ 1.20) | 1.09–1.27 (≤ 1.35) |
| todos | 0.5 | no viable en la malla (50/50 semillas) | $\sim10^{8}$–$10^{10}$ (no viable en 13/50) |

- Con $\delta$ correcto, la mediana del error de la regla queda entre 0.5 % y 14 % por encima del oráculo, sin consultar la verdad.
- Subestimar el ruido es mucho peor que sobrestimarlo.
  - Con $\delta=0.5\lVert e\rVert$, el umbral es $0.55\lVert e\rVert$. Para la semilla 2026 queda por debajo del residuo mínimo de Tikhonov en la malla, que es $0.696\lVert e\rVert$ en $\lambda/\sigma_1^2=10^{-10}$. La regla no encuentra $\lambda$ y lo informa.
  - Esta "no viabilidad" depende del extremo inferior de la malla. Como $A$ es invertible, con $\lambda$ aún menores el umbral se alcanzaría, pero con un $\lambda$ diminuto que amplifica el ruido. Eso es lo que le pasa a la TSVD: solo alcanza el umbral con $k$ muy grandes, que reproducen ruido, y el error resulta $10^{8}$ veces el del oráculo.
  - Las 13 semillas no viables de la TSVD también son un efecto numérico: el residuo calculado deja de ser monótono para $k\ge94$.
  - Con $\delta=2\lVert e\rVert$ se sobrerregulariza, pero la mediana del error sube solo entre 6 % y 27 % sobre el oráculo (máximo 35 %).

## 5. Por qué el error baja tan despacio al reducir el ruido

El error se concentra en el pulso. En los tres niveles de ruido, entre el 96 % y el 99 % del error cuadrático de Tikhonov (con $\lambda$ por discrepancia) está en $t\in[0.55,0.8]$, alrededor del pulso. Esa región es solo 32 de los 128 puntos.

Domina el sesgo, no el ruido. Separando $x_\lambda-x_{\mathrm{true}}$ en sesgo y ruido propagado (desarrollo, sección 7.2), para $\eta=10^{-3}$ se obtiene $\lVert\text{sesgo}\rVert/\lVert x_{\mathrm{true}}\rVert=0.108$ y $\lVert\text{ruido propagado}\rVert/\lVert x_{\mathrm{true}}\rVert=0.007$. Son vectores, así que sus normas no se suman, pero está claro cuál domina.

El error sí baja con menos ruido, pero lentamente (tabla `07_barrido_eta`):

| $\eta$ | Tikhonov oráculo (malla hasta $10^{-20}$) | TSVD oráculo ($k$) |
|---|---|---|
| $10^{-2}$ | 0.119 | 0.119 (27) |
| $10^{-3}$ | 0.101 | 0.103 (38) |
| $10^{-4}$ | 0.091 | 0.093 (40) |
| $10^{-6}$ | 0.077 | 0.078 (53) |
| $10^{-8}$ | 0.068 | 0.070 (65) |
| $10^{-10}$ | 0.065 | 0.065 (70) |
| $0$ | 0.061 | 0.054 (91) |

- Dividir el ruido por $10^{8}$ solo baja el error de 12 % a 6.5 %.
- La razón: los $\sigma_i$ decaen como una gaussiana en $i$, así que reducir el ruido apenas suma direcciones útiles. Pasar de $\eta=10^{-2}$ a $10^{-10}$ lleva el $k$ óptimo solo de 27 a 70.
- Las discontinuidades del pulso necesitan componentes $v_i$ de alta frecuencia, cuyos $\sigma_i$ quedan por debajo del nivel del ruido.
- Sin ruido, el límite lo pone la doble precisión, no la observación: $A$ es invertible en aritmética exacta. Las direcciones con $\sigma_i$ bajo el redondeo no se pueden usar.
- En la malla principal (hasta $10^{-10}$), el óptimo de Tikhonov para $\eta\le10^{-6}$ cae en el borde de la malla (7.8 %). Por eso la tabla usa una malla ampliada para este barrido.
- Es lo que anticipa el artículo: la regularización no recupera lo que el proceso de observación, junto con el ruido o el redondeo, dejó fuera del alcance de los datos.

## 6. Exploración: más suavizado, menos información

| $s$ | rango pinv | $\sigma_i$ sobre el ruido | error Tikhonov (discrepancia) | error oráculo |
|---|---|---|---|---|
| 0.02 | 128 | 65 | 0.085 | 0.079 |
| 0.03 | 89 | 44 | 0.109 | 0.101 |
| 0.05 | 56 | 27 | 0.138 | 0.126 |

Al aumentar $s$ los $\sigma_i$ caen más rápido, menos direcciones superan el nivel del ruido y el error mínimo alcanzable sube (fig. 09). Solo para $s=0.02$ el $\sigma_{\min}$ calculado está por encima del redondeo, y su $\kappa_2(A)\approx2.8\times10^{13}$ es fiable.

## 7. Ejemplos prácticos: imagen y sonido

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

## 8. Ampliación: penalización de diferencias (entrega final)

Tablas `08_gsvd`, `08_penalizaciones` y `08_robustez_penalizaciones`; figuras 15 a 17. Método en la sección 10 del desarrollo. Tres penalizaciones, $L=I$ (Tikhonov clásico), $L_1$ (primeras diferencias) y $L_2$ (segundas), sobre dos señales: la del experimento (gaussiana + pulso) y una suave, $x(t)=0.5+0.4\sin(2\pi t)+0.2t$. Malla absoluta de $\lambda$ en $[10^{-14},10^{4}]$; $\lambda$ por discrepancia con $\delta=\lVert e\rVert$.

Comprobaciones previas: $\ker A\cap\ker L=\{0\}$ para $L_1$ y $L_2$ (dimensiones de $\ker L$: 1 y 2); la GSVD calculada reconstruye $A$ con error $\le1.3\times10^{-14}$ y la fórmula espectral generalizada coincide con el sistema aumentado (QR) hasta $7\times10^{-8}$.

Mediana del error relativo en 50 semillas (discrepancia / oráculo):

| Señal | $\eta$ | $L=I$ | $L_1$ | $L_2$ |
|---|---|---|---|---|
| pulso | $10^{-4}$ | 0.0941 / 0.0927 | 0.0943 / 0.0927 | 0.0944 / 0.0926 |
| pulso | $10^{-3}$ | 0.1081 / 0.0996 | 0.1104 / 0.0994 | 0.1119 / 0.0992 |
| pulso | $10^{-2}$ | 0.1181 / 0.1174 | 0.1227 / 0.1170 | 0.1242 / 0.1171 |
| suave | $10^{-4}$ | 0.0044 / 0.0040 | 0.0025 / 0.0020 | 0.0002 / 0.0001 |
| suave | $10^{-3}$ | 0.0095 / 0.0094 | 0.0047 / 0.0039 | 0.0009 / 0.0005 |
| suave | $10^{-2}$ | 0.0235 / 0.0228 | 0.0124 / 0.0090 | 0.0064 / 0.0033 |

Interpretación:

1. Para el pulso, la penalización no cambia nada. El oráculo da el mismo error con las tres (diferencias menores del 0.5 %) y la discrepancia, errores dentro del 5 % entre sí, con $L=I$ ligeramente mejor. En la figura 15 (arriba a la derecha) las tres reconstrucciones se superponen punto a punto.
2. La razón está en los factores de filtro (figura 17): en el $\lambda$ elegido, las tres penalizaciones dejan pasar las mismas $\approx30$ direcciones generalizadas ($f_i>0.5$ para 31, 30 y 29 direcciones). El corte lo impone la caída de los valores singulares de $A$ frente al ruido; la penalización solo cambia la forma de la transición. Un salto no es "simple" para ninguna penalización cuadrática.
3. Para la señal suave, la penalización lo cambia todo. Con $\eta=10^{-3}$ el error baja de $0.95\,\%$ ($L=I$) a $0.47\,\%$ ($L_1$) y $0.09\,\%$ ($L_2$): diez veces. En la figura 15 (abajo a la derecha) el error punto a punto de $L=I$ oscila alrededor de $5\times10^{-3}$ y el de $L_2$ alrededor de $3\times10^{-4}$. La penalización $\lVert x\rVert^2$ empuja hacia cero una señal que no es pequeña (vale entre 0.25 y 0.95); $\lVert L_2x\rVert^2$ solo castiga la curvatura, que esta señal casi no tiene.
4. Los $\lambda$ elegidos son muy distintos ($10^{-4}$ para $I$, $10^{-2}$ para $L_1$, $10^{2}$ para $L_2$ con la señal suave): cada penalización tiene su propia escala, porque $\lVert Lx\rVert$ y $\lVert x\rVert$ miden cosas distintas. La discrepancia lo resuelve sola, porque mira el residuo y no $\lambda$.
5. Con $\lambda\to\infty$ y $L_1$, la reconstrucción tiende a la constante $\operatorname{media}(b)$, y con $L=I$, a cero (figura 16, extremo derecho): el error de $L=I$ se va a $1$, el de $L_1$ a $0.7$ (pulso).

Conclusión: la penalización es un modelo a priori de la señal; ayuda cuando describe lo que la señal es y no cuando no. El límite del pulso no está en el método sino en lo que la observación dejó en los datos. Para recuperar bordes haría falta una penalización no cuadrática (variación total), que sale del álgebra lineal.

## 9. Estimar el nivel de ruido a partir de los datos (entrega final)

Tabla `09_estimacion_delta`, figura 18. Método en la sección 8.4 del desarrollo: $\hat\delta=\sqrt{n\cdot\operatorname{media}_{i>64}(u_i^Tb)^2}$, usando la mitad final de los coeficientes de Picard, donde ya no hay señal.

| $\eta$ | $\hat\delta/\lVert e\rVert$ (mín, mediana, máx) | Tikhonov · disc. con $\lVert e\rVert$ | Tikhonov · disc. con $\hat\delta$ | oráculo | TSVD · disc. con $\hat\delta$ |
|---|---|---|---|---|---|
| $10^{-4}$ | 0.89, 0.99, 1.15 | 0.0941 | 0.0941 | 0.0927 | 0.0939 |
| $10^{-3}$ | 0.89, 0.99, 1.15 | 0.1081 | 0.1079 | 0.0996 | 0.1138 |
| $10^{-2}$ | 0.89, 0.99, 1.15 | 0.1181 | 0.1180 | 0.1174 | 0.1220 |

(medianas en 50 semillas; el cociente $\hat\delta/\lVert e\rVert$ no depende de $\eta$ porque la cola es ruido puro en los tres niveles)

- El estimador queda dentro de $\pm15\,\%$ del valor real, como predice su desviación teórica ($1/\sqrt{2\cdot64}\approx9\,\%$), y no fue inviable en ninguna semilla.
- La discrepancia con $\hat\delta$ da el mismo error que con $\lVert e\rVert$ exacto. Toda la cadena (SVD, estimación del ruido, elección de $\lambda$, reconstrucción) corre sin conocer $e$ ni $x_{\mathrm{true}}$; la verdad solo se usa para calificar.
- El refinamiento de $\lambda$ por bisección (tabla `09_lambda_refinado`) mueve $\lambda$ entre 1 % y 5 % y cambia el error en menos de $10^{-4}$: la malla ya era suficiente.

## 10. Limitaciones de estos resultados

- El ruido tiene norma fija y dirección uniforme en la esfera (isotrópico, "blanco"). Con ruido correlacionado o errores en $A$, las conclusiones pueden cambiar, y el estimador $\hat\delta$ supone ruido isotrópico.
- Los ejemplos de imagen y audio usan una sola semilla de ruido; no repetimos para ellos el estudio de robustez con 50 semillas.
- En la imagen, el desenfoque es gaussiano y separable con bordes renormalizados; otros desenfoques reales (movimiento en diagonal, desenfoque circular de una lente) no son separables y requerirían otra forma de calcular la SVD. En el audio, la convolución es circular (sección 9.2 del desarrollo).
- La malla de $\lambda$ es finita, con 6 puntos por década. Su extremo inferior influye en la viabilidad de la discrepancia y en el óptimo cuando hay poco ruido; por eso la sección 5 usa una malla ampliada.
- Varias cantidades calculadas están afectadas por el redondeo: el $\kappa_2$, la meseta de $\sigma_i$, el residuo de la pinv completa y la monotonía del residuo de la TSVD para $k\ge94$. Las hemos señalado donde aparecen.
- La ampliación compara penalizaciones cuadráticas; no probamos penalizaciones no cuadráticas (variación total), que son las que recuperan bordes.
- Un buen resultado con una señal sintética no garantiza el mismo desempeño con datos reales.
