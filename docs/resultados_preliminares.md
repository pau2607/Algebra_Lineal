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

Total: 23 pruebas, todas aprobadas (`python -m pytest -q`).

## 2. Espectro del operador

- $\sigma_1\approx1.001$ y los $\sigma_i$ decaen de forma aproximadamente gaussiana hasta $\sim10^{-14}$ en $i\approx90$; desde $i\approx95$ hay una meseta en $\sim10^{-16}$ que es redondeo de la SVD.
- El $\kappa_2(A)$ calculado es $2.2\times10^{16}$: $A$ es numéricamente singular. `scipy.linalg.pinv` (umbral $2.84\times10^{-14}$) y `numpy.linalg.matrix_rank` dan rango numérico **89**.
- En el gráfico de Picard (fig. 03), los $|u_i^Tb|$ siguen a $\sigma_i$ hasta que chocan con el nivel del ruido, $\approx\lVert e\rVert/\sqrt n$ ($4\times10^{-5}$, $4\times10^{-4}$ y $4\times10^{-3}$), y a partir de ahí quedan planos. Solo $\approx44$ valores singulares superan el nivel de ruido para $\eta=10^{-3}$ (tabla 04): esas son las direcciones que los datos realmente informan.

## 3. Comparación de métodos (semilla 2026)

Fragmento de `resultados/tablas/02_resumen_metodos.md` (error = error relativo):

| $\eta$ | Pinv numérica (rango 89) | TSVD discrepancia | TSVD oráculo | Tikhonov discrepancia | Tikhonov oráculo |
|---|---|---|---|---|---|
| $10^{-4}$ | $8.1\times10^{7}$ | 0.0939 ($k=38$) | 0.0931 ($k=40$) | 0.0941 ($\lambda/\sigma_1^2=6.8\times10^{-6}$) | 0.0914 ($3.2\times10^{-7}$) |
| $10^{-3}$ | $8.1\times10^{8}$ | 0.1139 ($k=31$) | 0.1029 ($k=38$) | 0.1087 ($3.2\times10^{-4}$) | 0.1008 ($2.2\times10^{-5}$) |
| $10^{-2}$ | $8.1\times10^{9}$ | 0.1218 ($k=23$) | 0.1185 ($k=27$) | 0.1194 ($4.6\times10^{-3}$) | 0.1191 ($3.2\times10^{-3}$) |

### Interpretación

1. **La pseudoinversa es inservible aunque ajuste los datos.** Su residuo está entre $0.54$ y $0.69\,\lVert e\rVert$, es decir, *más pequeño que el propio ruido*: está reproduciendo ruido. Su norma es $4\times10^{8}$ (la verdad tiene $\lVert x_{\mathrm{true}}\rVert\approx4.96$). El error crece exactamente $\times10$ cuando $\eta$ crece $\times10$ porque domina el término lineal $A^\dagger e$ y usamos la misma dirección $z$ para los tres niveles.
2. **Cuidado con el residuo de la pseudoinversa completa (128 $\sigma_i$).** En teoría debería ser casi cero, pero el calculado vale $\approx\lVert e\rVert$. Con $\lVert x\rVert\approx3\times10^{12}$, el redondeo al calcular $Ax$ es del orden de $\varepsilon\,\lVert A\rVert\,\lVert x\rVert\approx7\times10^{-4}$, comparable con $\lVert e\rVert$. Ese residuo no es confiable.
3. **La regularización reduce el error de $10^{8}$–$10^{10}$ a 9–12 %.** Las reconstrucciones (fig. 06) recuperan bien la zona suave y la posición del pulso, con oscilaciones cerca de sus bordes.
4. **La discrepancia sobrerregulariza un poco.** Elige $\lambda$ mayores (y $k$ menores) que el oráculo: el residuo del oráculo está en $0.77$–$0.94\,\lVert e\rVert$, por debajo del umbral $1.1\,\lVert e\rVert$. El costo en error es pequeño (ver sección 4).
5. **Tikhonov y TSVD dan resultados parecidos.** Coinciden en lo esencial porque los $\sigma_i$ caen muy rápido: el filtro suave y el corte abrupto descartan prácticamente las mismas direcciones (fig. 08).

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
- **Subestimar el ruido es mucho peor que sobrestimarlo.** Con $\delta=0.5\lVert e\rVert$ el umbral $0.55\lVert e\rVert$ está por debajo del residuo mínimo de Tikhonov en la malla ($\approx0.70\lVert e\rVert$ en $\lambda/\sigma_1^2=10^{-10}$), así que la regla no encuentra $\lambda$ y lo informa. En la TSVD el umbral solo se alcanza con $k$ muy grandes, que reproducen ruido: error de $10^{8}$ veces el oráculo. Con $\delta=2\lVert e\rVert$ se sobrerregulariza, pero la mediana del error sube solo entre 6 % y 27 % sobre el oráculo (máximo 35 %).

## 5. ¿Por qué el error no baja de ≈9 % aunque el ruido sea muy pequeño?

- En los tres niveles de ruido, entre el **96 % y el 99 %** del error cuadrático de Tikhonov (con $\lambda$ por discrepancia) se concentra en $t\in[0.55,0.8]$, alrededor del pulso, que es solo 32 de los 128 puntos.
- Separando $x_\lambda-x_{\mathrm{true}}$ en sesgo y ruido propagado (sección 7.2 del desarrollo), para $\eta=10^{-3}$ el sesgo aporta $\approx0.108$ del error relativo y el ruido solo $\approx0.007$.
- **Sin ruido** ($e=0$, tabla 03), el mejor error posible es 5.4 % (TSVD, $k=91$) y 7.8 % (Tikhonov en el extremo inferior de la malla). Las discontinuidades del pulso necesitan componentes $v_i$ de alta frecuencia cuyos $\sigma_i$ están por debajo de la precisión de máquina: esa información no está en los datos de forma utilizable. Es el límite que el artículo anticipa: la regularización no recupera lo que la observación destruyó.

## 6. Exploración: más suavizado, menos información

| $s$ | rango pinv | $\sigma_i$ sobre el ruido | error Tikhonov (discrepancia) | error oráculo |
|---|---|---|---|---|
| 0.02 | 128 | 65 | 0.085 | 0.079 |
| 0.03 | 89 | 44 | 0.109 | 0.101 |
| 0.05 | 56 | 27 | 0.138 | 0.126 |

Al aumentar $s$ los $\sigma_i$ caen más rápido, menos direcciones superan el nivel del ruido y el error mínimo alcanzable sube (fig. 09).

## 7. Limitaciones de estos resultados

- $\delta=\lVert e\rVert$ es exacto porque generamos el ruido; en datos reales habría que estimarlo, y la sección 4 muestra que subestimarlo es peligroso.
- El ruido es gaussiano blanco; con ruido correlacionado o errores en $A$ las conclusiones pueden cambiar.
- La malla de $\lambda$ es finita (6 puntos por década): el $\lambda$ elegido tiene resolución de $10^{1/6}\approx1.47$. En el caso sin ruido el óptimo de Tikhonov cae en el borde de la malla.
- Un buen resultado con una señal sintética no garantiza el mismo desempeño con datos reales.
