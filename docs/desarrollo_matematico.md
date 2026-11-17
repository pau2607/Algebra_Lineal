# Desarrollo matemático

**Proyecto:** Reconstrucción de señales con ruido pseudoinversa, valores singulares y regularización de Tikhonov
**Curso:** Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico
**Estado:** entrega final + ejemplos prácticos (imagen y audio) agregados después de la retroalimentación

Notación: $\lVert\cdot\rVert$ es la norma euclídea; $A^{\dagger}$, la pseudoinversa de Moore–Penrose; $\kappa_2(A)=\sigma_1/\sigma_n$, el número de condición en norma 2 (para columnas independientes; si $A$ es singular se toma $\kappa_2=\infty$). En Tikhonov, $\lambda$ es directamente el coeficiente de la penalización (no $\lambda^2$), igual que en el artículo de referencia [1].

## 1. Problema, objetivo y alcance

### 1.1 Problema

Una señal desconocida $x_{\mathrm{true}}\in\mathbb R^{n}$ se observa a través de un proceso lineal de suavizado y con ruido de medición:

$$
b = A x_{\mathrm{true}} + e,\qquad A\in\mathbb R^{m\times n}.
$$

Queremos recuperar $x_{\mathrm{true}}$ a partir de $b$ y $A$. La dificultad es que $A$ atenúa fuertemente los detalles finos (los bordes del pulso), de modo que su contribución a $b$ es del mismo orden que el ruido o menor.

### 1.2 Problema concreto que delimitamos

- Malla $t_j=j/(n-1)$ con $n=128$.
- Señal $x_j=\exp[-80(t_j-0.3)^2]+0.6\,\mathbf 1_{[0.6,\,0.75]}(t_j)$: una zona suave y un pulso con dos discontinuidades.
- Suavizado gaussiano con filas normalizadas, $s=0.03$:

$$
\tilde a_{ij}=\exp\!\left(-\frac{(t_i-t_j)^2}{2s^2}\right),\qquad a_{ij}=\frac{\tilde a_{ij}}{\sum_{\ell}\tilde a_{i\ell}} .
$$

- Ruido $e=\eta\,\lVert Ax_{\mathrm{true}}\rVert\, z/\lVert z\rVert$ con $z$ gaussiano de semilla 2026 y $\eta\in\{10^{-4},10^{-3},10^{-2}\}$.

### 1.3 Objetivo verificable

1. Reproducir el ejemplo diagonal de dos coordenadas y comprobar que el código coincide con los cálculos a mano (las soluciones $x$ con tolerancia $10^{-10}$; los errores, que el artículo da redondeados, con $5\times10^{-4}$).
2. Para cada $\eta$, comparar pseudoinversa numérica, SVD truncada y Tikhonov midiendo error relativo $\lVert x-x_{\mathrm{true}}\rVert/\lVert x_{\mathrm{true}}\rVert$, residuo $\lVert Ax-b\rVert$ y norma $\lVert x\rVert$.
3. Elegir el parámetro con una regla que no usa $x_{\mathrm{true}}$ (principio de discrepancia) y cuantificar cuánto se aleja del mejor parámetro posible (oráculo), con 50 semillas.
4. Explicar con la SVD por qué la pseudoinversa falla y qué parte del error no puede eliminarse regularizando.

### 1.4 Alcance

Incluye: señal unidimensional con las tres reconstrucciones, SVD calculada con NumPy (no se programa desde cero), principio de discrepancia con $\delta$ conocido y con $\delta$ estimado de los datos (sección 8.4), ejemplos prácticos con una imagen y un sonido (sección 9) y, como ampliación, la penalización de diferencias $\lambda\lVert Lx\rVert^2$ (sección 10). La curva L solo se grafica.

Después de la entrega parcial, y siguiendo la retroalimentación del docente, agregamos dos ejemplos prácticos con el mismo método (sección 9): una imagen de $256\times256$ desenfocada (fotografía y texto) y un audio de 2 s con los agudos apagados. Se pueden ver y escuchar en `resultados/imagen/` y `resultados/audio/`.

## 2. Conceptos del curso que se utilizan

| Concepto | Dónde aparece en el proyecto |
|---|---|
| Mínimos cuadrados y ecuaciones normales | Planteamiento del ajuste; por qué no conviene usar $A^TA$ |
| Núcleo, imagen y rango | Solución de mínima norma; componentes que ningún método recupera |
| Bases ortonormales y proyección ortogonal | Coordenadas en $\{u_i\}$ y $\{v_i\}$; residuo fuera de la imagen |
| Formas cuadráticas definidas positivas | Unicidad de la solución de Tikhonov |
| SVD y pseudoinversa | Fórmulas de las tres reconstrucciones, filtrado espectral |
| Número de condición | Amplificación del ruido |
| Producto de Kronecker y vectorización | Desenfoque de imágenes como $b=Ax+e$; SVD de $A_h\otimes A_v$ (sección 9.1) |
| Números complejos, producto interno complejo y Fourier | Las matrices circulantes se diagonalizan con la DFT; filtrado del audio por frecuencias (sección 9.2) |

## 3. Mínimos cuadrados y pseudoinversa

El ajuste básico es

$$
\min_{x}\ \lVert Ax-b\rVert^2 .
$$

Sus ecuaciones normales son $A^TAx=A^Tb$. No las usamos para calcular [8]: para columnas independientes $\kappa_2(A^TA)=\kappa_2(A)^2$, porque los valores singulares de $A^TA$ son los $\sigma_i^2$.

¿Es $A$ singular? En aritmética exacta, no: $\tilde A$ es la matriz de un núcleo gaussiano evaluado en puntos distintos, que es simétrica definida positiva, y $A=D^{-1}\tilde A$ con $D=\operatorname{diag}(\sum_\ell\tilde a_{i\ell})$ diagonal positiva; luego $A$ es invertible y $\ker A=\{0\}$. Pero sus $\sigma_i$ verdaderos más pequeños están muy por debajo de $\varepsilon_{\text{máq}}\,\sigma_1\approx2.2\times10^{-16}$, así que en doble precisión $A$ es numéricamente singular. El $\sigma_{\min}$ que devuelve la SVD ($\approx5\times10^{-17}$) es nivel de redondeo, y el $\kappa_2(A)\approx2\times10^{16}$ calculado no es el número de condición verdadero: solo indica que $\kappa_2(A)\gtrsim1/\varepsilon_{\text{máq}}$ (otros algoritmos de SVD dan valores distintos del mismo orden o mayores). Por eso no reportamos $\kappa_2(A)^2$.

Qué se pierde al formar $A^TA$. En doble precisión se pierde todo $\sigma_i^2$ menor que $\varepsilon_{\text{máq}}\,\sigma_1^2$, es decir, todo $\sigma_i\lesssim\sqrt{\varepsilon_{\text{máq}}}\,\sigma_1\approx1.5\times10^{-8}$. En nuestra $A$ solo 67 valores singulares superan ese nivel, y el rango numérico de la $A^TA$ calculada es 63, frente a 89 para $A$ (tabla `07_ecuaciones_normales`). Trabajar con la SVD de $A$ (o con QR) evita esa pérdida. Con regularización la situación mejora: $\kappa_2(A^TA+\lambda I)\approx\sigma_1^2/\lambda\le10^{10}$ en nuestra malla, y las ecuaciones normales regularizadas difieren de la fórmula espectral en $4.5\times10^{-7}$ como máximo. La advertencia vale sobre todo para mínimos cuadrados sin regularizar.

### 3.1 Solución de mínima norma

Sea la SVD compacta de rango $r$:

$$
A=U_r\Sigma_rV_r^T=\sum_{i=1}^{r}\sigma_iu_iv_i^T,\qquad \sigma_1\ge\cdots\ge\sigma_r>0 .
$$

**Proposición 1.** La solución de mínimos cuadrados de norma mínima es

$$
x^{\dagger}=A^{\dagger}b=\sum_{i=1}^{r}\frac{u_i^Tb}{\sigma_i}\,v_i .
$$

*Demostración.* Escribimos $x=\sum_{i=1}^r z_iv_i+x_0$ con $x_0\in\ker A$ (los $v_i$ con $i\le r$ forman una base ortonormal de $(\ker A)^{\perp}$). Entonces $Ax=\sum_i\sigma_iz_iu_i$. Completamos $\{u_i\}$ a una base ortonormal de $\mathbb R^m$ y separamos $b=\sum_{i\le r}(u_i^Tb)u_i+b_{\perp}$, donde $b_\perp$ es ortogonal a la imagen de $A$. Por Pitágoras,

$$
\lVert Ax-b\rVert^2=\sum_{i=1}^{r}\big(\sigma_iz_i-u_i^Tb\big)^2+\lVert b_\perp\rVert^2 .
$$

Cada sumando se anula con $z_i=u_i^Tb/\sigma_i$; $\lVert b_\perp\rVert^2$ no depende de $x$. Como $x_0$ no cambia el residuo y $\lVert x\rVert^2=\sum z_i^2+\lVert x_0\rVert^2$, la norma mínima se obtiene con $x_0=0$. $\square$

Observación. $x^\dagger$ es una elección entre todos los ajustes óptimos. Si $x_{\mathrm{true}}$ tiene componente en $\ker A$, ningún método basado solo en $b$ puede recuperarla. En nuestro problema $\ker A=\{0\}$ en aritmética exacta; lo que se pierde en la práctica son las componentes $v_i^Tx_{\mathrm{true}}$ cuyo $\sigma_i$ queda por debajo del nivel del ruido o del redondeo.

## 4. Por qué se amplifica el ruido

Por linealidad, $x^\dagger = A^\dagger Ax_{\mathrm{true}} + A^\dagger e$, y la parte que viene del ruido es

$$
A^{\dagger}e=\sum_{i=1}^{r}\frac{u_i^Te}{\sigma_i}\,v_i ,
\qquad
\lVert A^\dagger e\rVert=\Big(\sum_{i=1}^{r} \frac{(u_i^Te)^2}{\sigma_i^2}\Big)^{1/2}\le \frac{1}{\sigma_r}\Big(\sum_{i=1}^{r}(u_i^Te)^2\Big)^{1/2}\le \frac{\lVert e\rVert}{\sigma_r},
$$

donde la última desigualdad es la de Bessel ($\{u_i\}_{i\le r}$ es ortonormal).

El coeficiente de $v_i$ es el coeficiente del ruido en la dirección $u_i$ multiplicado por $1/\sigma_i$. Nuestro ruido $e=\eta\lVert Ax_{\mathrm{true}}\rVert z/\lVert z\rVert$ tiene norma fija y dirección uniforme en la esfera (es isotrópico, aunque no exactamente gaussiano). Si completamos $\{u_i\}$ a una base ortonormal de $\mathbb R^m$, $\sum_{i=1}^{m}(u_i^Te)^2=\lVert e\rVert^2$ y, por simetría, $\mathbb E[(u_i^Te)^2]=\lVert e\rVert^2/m$: los $u_i^Te$ tienen todos el mismo tamaño típico $\lVert e\rVert/\sqrt m$ (cada uno fluctúa alrededor de ese valor; aquí $m=n=128$). En cambio, para datos limpios $u_i^T(Ax_{\mathrm{true}})=\sigma_i\,v_i^Tx_{\mathrm{true}}$ decae al menos tan rápido como $\sigma_i$ (condición discreta de Picard [3]). Hay entonces un índice a partir del cual los coeficientes de $b$ están dominados por el ruido y dividirlos entre $\sigma_i$ produce oscilaciones enormes. Esto es lo que muestra el gráfico de Picard (`resultados/figuras/03_espectro_picard.png`).

### 4.1 Ejemplo diagonal (cálculo a mano)

$$
A=\begin{pmatrix}1&0\\0&10^{-3}\end{pmatrix},\quad x_{\mathrm{true}}=\begin{pmatrix}1\\1\end{pmatrix},\quad e=\begin{pmatrix}0\\0.01\end{pmatrix},\quad b=\begin{pmatrix}1\\0.011\end{pmatrix}.
$$

Aquí $u_i=v_i=e_i$, $\sigma_1=1$, $\sigma_2=10^{-3}$.

- Pseudoinversa: $x^\dagger=(1/1,\ 0.011/10^{-3})=(1,11)$. Residuo $0$, error $\lVert(0,10)\rVert=10$. La dirección 2 multiplica el ruido por $1/\sigma_2=1000$: un error de $0.01$ en el dato se vuelve un error de $10$ en la solución.
- SVD truncada, $k=1$: $x_1=(1,0)$. Error $1$, residuo $\lVert(0,-0.011)\rVert=0.011$.
- Tikhonov, $\lambda=10^{-4}$: $x_\lambda=\Big(\tfrac{1}{1+10^{-4}},\ \tfrac{10^{-3}\cdot0.011}{10^{-6}+10^{-4}}\Big)=\Big(\tfrac{1}{1.0001},\tfrac{0.000011}{0.000101}\Big)\approx(0.99990,\ 0.10891)$. Error $\approx0.8911$.

Conclusión: un residuo pequeño no certifica una buena reconstrucción. El script `scripts/01_ejemplo_diagonal.py` y las pruebas `test_diagonal_*` verifican estos números.

## 5. SVD truncada

Para $k\le r$,

$$
x_k=\sum_{i=1}^{k}\frac{u_i^Tb}{\sigma_i}\,v_i .
$$

Se descartan las direcciones con $\sigma_i$ pequeño: mejora la estabilidad, pero se pierden las componentes verdaderas $v_i^Tx_{\mathrm{true}}$ con $i>k$. El residuo es

$$
\lVert Ax_k-b\rVert^2=\sum_{i=k+1}^{r}(u_i^Tb)^2+\lVert b_\perp\rVert^2 ,
$$

que no crece al aumentar $k$ en aritmética exacta (lo comprueba `test_residuo_tsvd_no_crece_con_k` en una matriz bien escalada). En el cálculo con nuestra $A$, el residuo calculado deja de ser monótono desde $k=94$, cuando entran direcciones con $\sigma_i$ al nivel del redondeo (tabla `07_diagnosticos`).

### 5.1 Pseudoinversa numérica ≠ inversa sin regularizar

`scipy.linalg.pinv` descarta los $\sigma_i<\texttt{rtol}\cdot\sigma_1$ con $\texttt{rtol}=\max(m,n)\,\varepsilon$ por defecto [5]. En nuestro problema ese umbral es $\approx2.84\times10^{-14}$ y conserva 89 de 128 valores singulares. Es decir, la "pseudoinversa numérica" ya es una SVD truncada con $k=89$; lo informamos explícitamente. Para comparar, también calculamos la suma sobre los 128 $\sigma_i$ positivos ("pseudoinversa completa"). Desde $i=98$ los $\sigma_i$ calculados forman una meseta en $\approx7\times10^{-17}$ ($\approx0.3\,\varepsilon_{\text{máq}}\sigma_1$): los verdaderos son todavía menores y lo calculado refleja solo el redondeo, no información del operador.

## 6. Regularización de Tikhonov

La idea, debida a Tikhonov [6] y estudiada para problemas discretos en [2], [3], es penalizar la norma de la solución:

$$
x_\lambda=\arg\min_x\{\lVert Ax-b\rVert^2+\lambda\lVert x\rVert^2\},\qquad \lambda>0 .
$$

Teorema 2 (unicidad y fórmula espectral). Para todo $\lambda>0$ el problema tiene solución única

$$
x_\lambda=(A^TA+\lambda I)^{-1}A^Tb=\sum_{i=1}^{r}\frac{\sigma_i}{\sigma_i^2+\lambda}\,(u_i^Tb)\,v_i .
$$

Demostración. La función $J(x)=\lVert Ax-b\rVert^2+\lambda\lVert x\rVert^2$ tiene gradiente $\nabla J=2A^T(Ax-b)+2\lambda x$ y hessiana $2(A^TA+\lambda I)$. Para $z\neq0$,

$$
z^T(A^TA+\lambda I)z=\lVert Az\rVert^2+\lambda\lVert z\rVert^2\ \ge\ \lambda\lVert z\rVert^2>0,
$$

así que la matriz es definida positiva: es invertible, $J$ es estrictamente convexa y el único punto estacionario $(A^TA+\lambda I)x=A^Tb$ es el mínimo global. Para la fórmula, escribimos $x=\sum_{i=1}^n z_iv_i$ (base ortonormal completa $V$). Como $A^TAv_i=\sigma_i^2v_i$ y $A^Tb=\sum_{i\le r}\sigma_i(u_i^Tb)v_i$, comparando coeficientes:

$$
(\sigma_i^2+\lambda)z_i=\sigma_i(u_i^Tb)\quad(i\le r),\qquad \lambda z_i=0\quad(i>r).
$$

Luego $z_i=\sigma_i(u_i^Tb)/(\sigma_i^2+\lambda)$ y las componentes del núcleo se anulan. $\square$

Lo que la positividad usa es solo $\lambda>0$: aunque $A$ tenga núcleo no trivial, $\lambda\lVert z\rVert^2$ es positivo para todo $z\ne0$. (Con una penalización $\lambda\lVert Lx\rVert^2$ haría falta $\ker A\cap\ker L=\{0\}$.)

### 6.1 Factores de filtro

$$
x_\lambda=\sum_{i=1}^{r}f_i(\lambda)\,\frac{u_i^Tb}{\sigma_i}\,v_i,\qquad f_i(\lambda)=\frac{\sigma_i^2}{\sigma_i^2+\lambda}\in(0,1).
$$

- Si $\sigma_i^2\gg\lambda$, $f_i\approx1$: la dirección se conserva.
- Si $\sigma_i^2\ll\lambda$, $f_i\approx\sigma_i^2/\lambda$: la dirección se atenúa y el término queda $\approx\sigma_i(u_i^Tb)/\lambda$, que no explota.

La SVD truncada corresponde a $f_i\in\{0,1\}$ (corte abrupto); Tikhonov es un filtro suave [2] (`resultados/figuras/08_factores_filtro.png`).

Calculamos $\sigma_i/(\sigma_i^2+\lambda)\cdot(u_i^Tb)$ en lugar de $f_i\cdot(u_i^Tb)/\sigma_i$ para no dividir entre $\sigma_i\approx10^{-17}$. Contrastamos con el sistema aumentado

$$
\min_x\left\lVert\begin{pmatrix}A\\ \sqrt\lambda\,I\end{pmatrix}x-\begin{pmatrix}b\\0\end{pmatrix}\right\rVert^2,
$$

cuyas ecuaciones normales son exactamente $(A^TA+\lambda I)x=A^Tb$, resuelto con QR sin formar $A^TA$. En 13 valores de la malla (uno de cada cinco, con $\eta=10^{-3}$) la diferencia relativa máxima es $\approx3\times10^{-10}$.

### 6.2 Residuo como función de $\lambda$

$$
\lVert Ax_\lambda-b\rVert^2=\sum_{i=1}^{r}\left(\frac{\lambda}{\sigma_i^2+\lambda}\right)^2(u_i^Tb)^2+\lVert b_\perp\rVert^2 .
$$

Cada factor $\lambda/(\sigma_i^2+\lambda)$ es creciente en $\lambda$, así que el residuo es no decreciente en $\lambda$ (estrictamente creciente si algún $u_i^Tb\neq0$ con $i\le r$). Esto hace que la regla de discrepancia de la sección 8 esté bien definida.

## 7. Estabilidad, sesgo y descomposición del error

### 7.1 Cota de estabilidad

Si los datos cambian de $b$ a $b+e$, por linealidad

$$
\lVert x_\lambda(b+e)-x_\lambda(b)\rVert=\Big\lVert\sum_i\frac{\sigma_i}{\sigma_i^2+\lambda}(u_i^Te)v_i\Big\rVert\le\max_{s\ge0}\frac{s}{s^2+\lambda}\,\lVert e\rVert=\frac{\lVert e\rVert}{2\sqrt\lambda},
$$

usando la ortonormalidad de los $v_i$ y otra vez la desigualdad de Bessel.

El máximo se obtiene de $(s-\sqrt\lambda)^2\ge0\iff s^2+\lambda\ge2s\sqrt\lambda\iff \dfrac{s}{s^2+\lambda}\le\dfrac{1}{2\sqrt\lambda}$, con igualdad en $s=\sqrt\lambda$. Comparado con la pseudoinversa (amplificación $1/\sigma_r$), Tikhonov amplifica a lo sumo $1/(2\sqrt\lambda)$, independiente de cuán pequeño sea $\sigma_r$.

### 7.2 Descomposición sesgo + ruido

Con $b=Ax_{\mathrm{true}}+e$ y $u_i^TAx_{\mathrm{true}}=\sigma_iv_i^Tx_{\mathrm{true}}$:

$$
x_\lambda-x_{\mathrm{true}}=\underbrace{\sum_{i=1}^{r}\big(f_i(\lambda)-1\big)(v_i^Tx_{\mathrm{true}})\,v_i\;-\;P_{\ker A}\,x_{\mathrm{true}}}_{\text{sesgo (no depende del ruido)}}\;+\;\underbrace{\sum_{i=1}^{r}f_i(\lambda)\frac{u_i^Te}{\sigma_i}\,v_i}_{\text{ruido propagado}} .
$$

- Al aumentar $\lambda$: cada $f_i$ baja, así que la norma del término de ruido, $\big(\sum_i f_i^2(u_i^Te)^2/\sigma_i^2\big)^{1/2}$, baja (y también su cota $\lVert e\rVert/(2\sqrt\lambda)$), mientras que la del sesgo, $\big(\sum_i(1-f_i)^2(v_i^Tx_{\mathrm{true}})^2+\lVert P_{\ker A}x_{\mathrm{true}}\rVert^2\big)^{1/2}$, sube.
- Al disminuir $\lambda$: lo contrario.
- $P_{\ker A}x_{\mathrm{true}}$ no depende de $\lambda$: es información perdida. En nuestro problema es cero ($A$ es invertible), pero el mismo papel lo cumplen en la práctica las direcciones con $\sigma_i$ diminuto, cuyo $f_i\approx0$ para cualquier $\lambda$ razonable.

Límites: $\lambda\to0^+$ da $x_\lambda\to A^\dagger b$; $\lambda\to\infty$ da $x_\lambda\to0$ (`test_limites_de_lambda`). El mejor $\lambda$ depende de la escala de $A$, del ruido y de la señal; por eso la malla se expresa como $\lambda/\sigma_1^2\in[10^{-10},1]$.

## 8. Elección del parámetro sin conocer la verdad

### 8.1 Oráculo

En un experimento sintético podemos calcular el error para cada $\lambda$ de la malla y tomar el mínimo. Esto usa $x_{\mathrm{true}}$, que no existe en un problema real. Lo usamos solo como referencia para evaluar la regla que sí es aplicable.

### 8.2 Principio de discrepancia (Morozov [7])

Si se conoce una cota $\delta\approx\lVert e\rVert$, no tiene sentido ajustar los datos mejor que el nivel del ruido: la propia señal verdadera tiene residuo $\lVert Ax_{\mathrm{true}}-b\rVert=\lVert e\rVert$, así que un residuo claramente menor que $\lVert e\rVert$ indica que la solución está reproduciendo parte del ruido. La regla es:

> Elegir el mayor $\lambda$ de la malla tal que $\lVert Ax_\lambda-b\rVert\le\tau\delta$, con $\tau=1.1$.
> Para la SVD truncada: el menor $k$ tal que $\lVert Ax_k-b\rVert\le\tau\delta$.

Por la monotonía de la sección 6.2, el conjunto de $\lambda$ viables es un intervalo inicial de la malla y la regla elige su extremo derecho. Para la TSVD, en aritmética exacta los $k$ viables forman un intervalo final; numéricamente eso puede fallar para $k\ge94$ (sección 5), por lo que el código no supone monotonía: busca el máximo $\lambda$ (o el mínimo $k$) entre todos los valores viables. Supuestos: (i) $\delta$ es una buena estimación del ruido; (ii) los errores de modelo son menores que el ruido; (iii) el umbral es alcanzable en la malla, lo que exige $\min_\lambda\lVert Ax_\lambda-b\rVert\le\tau\delta$. Si no se cumple (iii), la función devuelve "no viable" en lugar de inventar un valor. Ojo: como $A$ es invertible, en aritmética exacta el residuo tiende a $0$ cuando $\lambda\to0^+$; con una malla que bajara más el umbral siempre se alcanzaría, pero con un $\lambda$ diminuto que amplifica el ruido. La "no viabilidad" depende del extremo de la malla y es una señal de alarma útil, no una propiedad del problema. En el código, las funciones de discrepancia no reciben $x_{\mathrm{true}}$ (lo comprueba `test_discrepancia_no_recibe_la_verdad`).

En este experimento tomamos $\delta=\lVert e\rVert$ (conocido porque generamos el ruido) y estudiamos aparte qué ocurre si $\delta$ está mal estimado ($\delta=0.5\lVert e\rVert$ y $2\lVert e\rVert$).

### 8.3 Curva L (exploratoria)

Dibujamos $\log\lVert Ax_\lambda-b\rVert$ frente a $\log\lVert x_\lambda\rVert$. Su esquina es una heurística de compromiso; en esta entrega solo se grafica, no se usa para elegir.

### 8.4 Estimar $\delta$ a partir de los datos (entrega final)

En un problema real no se conoce $\lVert e\rVert$. Pero el gráfico de Picard sugiere cómo estimarlo: en las últimas direcciones los $\sigma_i$ son tan pequeños que $u_i^T(Ax_{\mathrm{true}})=\sigma_iv_i^Tx_{\mathrm{true}}$ es despreciable y, por tanto, $u_i^Tb\approx u_i^Te$: la cola de los coeficientes es solo ruido. Por la sección 4, $\mathbb E[(u_i^Te)^2]=\lVert e\rVert^2/m$, así que

$$
\hat\delta=\sqrt{m\cdot\operatorname{media}_{i>i_0}\,(u_i^Tb)^2},\qquad i_0=n/2 .
$$

Supuesto: las direcciones $i>i_0$ no contienen señal. Con nuestra $A$, $\sigma_{64}\approx10^{-8}$, así que vale para cualquier $\eta\ge10^{-8}$. Con $q=64$ coeficientes en la cola, la desviación relativa del estimador es $\approx1/\sqrt{2q}\approx9\,\%$: en 50 semillas, $\hat\delta/\lVert e\rVert$ queda entre $0.89$ y $1.15$ (tabla `09_estimacion_delta`). Usar $\hat\delta$ en lugar de $\lVert e\rVert$ no cambia la mediana del error de la discrepancia (sección 4 de los resultados). Con esto, toda la cadena de reconstrucción se puede ejecutar sin conocer ni $e$ ni $x_{\mathrm{true}}$; la verdad solo se usa para calificar.

**Refinamiento de $\lambda$.** La malla tiene resolución $10^{1/6}\approx1.47$. Como el residuo es no decreciente en $\lambda$ (sección 6.2), entre el último $\lambda$ viable y el siguiente hay exactamente un punto donde $\lVert Ax_\lambda-b\rVert=\tau\delta$, que se encuentra por bisección (`discrepancia_tikhonov_refinada`). El cambio en el error es menor que $10^{-4}$: la malla ya era suficiente.

## 9. Ejemplos prácticos: una imagen y un sonido

Incorporamos la retroalimentación del docente a la entrega parcial: mostrar la reconstrucción sobre datos que se puedan *ver* y *escuchar*. No cambiamos de método: usamos las mismas reconstrucciones y la misma regla de discrepancia de las secciones 3 a 8 (en la imagen, las tres: pseudoinversa, SVD truncada y Tikhonov; en el audio, la pseudoinversa completa y Tikhonov). Lo nuevo es cómo obtener la SVD cuando la matriz es demasiado grande para formarla. Usamos un suavizado con el mismo ancho del experimento principal ($0.03\times127\approx3.8$ píxeles o muestras), así que la "dificultad" del operador es comparable.

### 9.1 Imagen: desenfoque separable y producto de Kronecker

Una imagen en escala de grises es una matriz $X\in\mathbb R^{n\times n}$ ($n=256$, valores en $[0,1]$). Desenfocamos cada columna con $A_v$ (dirección vertical) y cada fila con $A_h$ (horizontal), ambas construidas con `datos.matriz_suavizado`:

$$
B=A_vXA_h^T+E .
$$

**Vectorización.** Si $\operatorname{vec}(X)\in\mathbb R^{n^2}$ apila las columnas de $X$, entonces

$$
\operatorname{vec}(A_vXA_h^T)=(A_h\otimes A_v)\operatorname{vec}(X),
$$

porque la columna $j$ de $A_vXA_h^T$ es $\sum_\ell (A_h)_{j\ell}\,A_vx_\ell$, y el bloque $(j,\ell)$ de $A_h\otimes A_v$ es justamente $(A_h)_{j\ell}A_v$. Así, la imagen es un caso del problema $b=Ax+e$ con $A=A_h\otimes A_v\in\mathbb R^{65\,536\times65\,536}$. Esa matriz ocuparía unos 34 GB en doble precisión, y su SVD directa sería inviable.

**Proposición 3 (SVD de un producto de Kronecker).** Si $A_v=U_v\Sigma_vV_v^T$ y $A_h=U_h\Sigma_hV_h^T$ son SVD, entonces

$$
A_h\otimes A_v=(U_h\otimes U_v)\,(\Sigma_h\otimes\Sigma_v)\,(V_h\otimes V_v)^T
$$

es una SVD (salvo el orden de los valores singulares), y los valores singulares de $A_h\otimes A_v$ son los productos $S_{ij}=\sigma_i(A_v)\,\sigma_j(A_h)$.

*Demostración.* Por la propiedad del producto mixto $(P\otimes Q)(R\otimes T)=PR\otimes QT$,

$$
(U_h\otimes U_v)(\Sigma_h\otimes\Sigma_v)(V_h^T\otimes V_v^T)=U_h\Sigma_hV_h^T\otimes U_v\Sigma_vV_v^T=A_h\otimes A_v .
$$

Además, $(U_h\otimes U_v)^T(U_h\otimes U_v)=U_h^TU_h\otimes U_v^TU_v=I\otimes I=I$, y lo mismo para $V_h\otimes V_v$, así que ambas son ortogonales. $\Sigma_h\otimes\Sigma_v$ es diagonal con entradas $\sigma_j(A_h)\sigma_i(A_v)\ge0$. Reordenando las columnas se obtienen los valores singulares en orden decreciente. $\square$

El par de vectores singulares asociado a $S_{ij}$ es $u^h_j\otimes u^v_i$ y $v^h_j\otimes v^v_i$, que como imágenes son las matrices de rango uno $u^v_i(u^h_j)^T$ y $v^v_i(v^h_j)^T$. El coeficiente de Picard es

$$
(u^h_j\otimes u^v_i)^T\operatorname{vec}(B)=(u^v_i)^TB\,u^h_j=C_{ij},\qquad C=U_v^TBU_h .
$$

**Consecuencia.** Las fórmulas de las secciones 3, 5 y 6 se aplican entrada por entrada sobre la matriz $C$ ($\circ$ es el producto entrada por entrada):

$$
X^\dagger=V_v\big(C\oslash S\big)V_h^T,\qquad
X_k=V_v\big(M_k\circ C\oslash S\big)V_h^T,\qquad
X_\lambda=V_v\Big(\frac{S}{S^{\circ2}+\lambda}\circ C\Big)V_h^T ,
$$

donde $M_k$ es la máscara de los $k$ valores $S_{ij}$ más grandes. Solo se calculan dos SVD de $256\times256$ y productos de matrices de $256\times256$. Como $U_v$ y $U_h$ son cuadradas, el residuo de Tikhonov es $\lVert(\lambda/(S^{\circ2}+\lambda))\circ C\rVert_F$, que es la sección 6.2 escrita en forma matricial. Por eso la discrepancia sigue bien definida.

**Por qué la imagen es más difícil.** $\kappa_2(A_h\otimes A_v)=\kappa_2(A_h)\,\kappa_2(A_v)$: el condicionamiento *se eleva al cuadrado* al pasar a dos dimensiones. Los $S_{ij}$ calculados llegan a $\sim10^{-34}$. El umbral por defecto de `pinv`, $\max(m,n)\,\varepsilon\,S_{\max}$ con $m=n=65\,536$, conserva 18 817 de las 65 536 direcciones, y aun así la reconstrucción es ruido puro (sección 7 de los resultados).

**Comprobación.** `tests/test_aplicaciones.py` forma $A_h\otimes A_v$ explícitamente para un caso $9\times7$ y verifica que los valores singulares son los productos, y que pinv, SVD truncada y Tikhonov en forma matricial coinciden con `metodos.py` aplicado a la matriz grande y con el sistema aumentado (QR).

### 9.2 Sonido: suavizado circular y transformada de Fourier

Un audio de 2.048 s a 8000 muestras por segundo es un vector $x\in\mathbb R^{N}$ con $N=16\,384$. Tampoco formamos $A$ ($N\times N$).

**Primer intento descartado.** Reutilizar la misma $A$ de $128\times128$ por bloques de 128 muestras (16 ms) es exactamente el problema 1D aplicado a cada columna. Es el caso $A_h=I$ de la sección 9.1, y lo comprueba `test_audio_por_bloques_es_el_caso_1d_columna_a_columna`. Pero como modelo de observación es malo: al suavizar cada bloque por separado aparecen saltos de hasta $0.58$ entre muestras vecinas en los bordes de los bloques (con el suavizado continuo, a lo sumo $0.06$), 62.5 veces por segundo, y el audio "observado" se oiría con un zumbido que un suavizado real no produce. El resultado se aparta hasta $0.26$ del suavizado continuo, un 14 % en norma. Por eso usamos un suavizado continuo sobre toda la señal.

**Modelo.** $Ax=h\circledast x$, convolución circular con el núcleo gaussiano $h_j=\exp(-d_j^2/2\sigma^2)/\sum_\ell\exp(-d_\ell^2/2\sigma^2)$, $d_j=\min(j,N-j)$, $\sigma=3.81$ muestras. Como $h$ suma 1, se conservan las constantes, igual que en las filas de la $A$ del experimento principal. $A$ es una matriz *circulante*: cada fila es la anterior desplazada un lugar.

**Proposición 4 (las circulantes se diagonalizan con Fourier).** Sea $F$ la matriz de la transformada discreta de Fourier, $F_{k\ell}=\omega^{k\ell}$ con $\omega=e^{-2\pi\mathrm i/N}$, y $\hat h=Fh$. Entonces $A=F^{-1}\operatorname{diag}(\hat h)\,F$. Si además $h$ es simétrico ($h_j=h_{N-j}$), $\hat h$ es real y los valores singulares de $A$ son $|\hat h_k|$.

*Demostración.* La columna $\ell$ de $F^{-1}$ es $f_\ell=\frac1N(\omega^{-k\ell})_k$, una oscilación pura. Al convolucionar $h$ con una oscilación pura se obtiene la misma oscilación multiplicada por $\hat h_\ell$, así que $Af_\ell=\hat h_\ell f_\ell$: las columnas de $F^{-1}$ son vectores propios. Con $\tilde F=F/\sqrt N$ unitaria ($\tilde F^*\tilde F=I$, producto interno complejo de la Unidad 2 del curso) se tiene $A=\tilde F^*\operatorname{diag}(\hat h)\tilde F$. Si $h$ es real y simétrico, $\hat h_k=\sum_jh_j\cos(2\pi jk/N)\in\mathbb R$. Escribiendo $\hat h_k=|\hat h_k|\operatorname{signo}(\hat h_k)$ y absorbiendo el signo en uno de los factores unitarios, se obtiene una SVD con valores singulares $|\hat h_k|$. $\square$

Con esto, la fórmula de Tikhonov (Teorema 2, válido también con $A^*$ en lugar de $A^T$) queda, frecuencia por frecuencia,

$$
\hat x_{\lambda,k}=\frac{\hat h_k}{\hat h_k^2+\lambda}\,\hat b_k,\qquad f_k(\lambda)=\frac{\hat h_k^2}{\hat h_k^2+\lambda},
$$

y el residuo también se puede calcular en frecuencia gracias a la identidad de Parseval, $\lVert r\rVert=\lVert\hat r\rVert/\sqrt N$ (el código lo calcula en el tiempo, que da lo mismo). Todo cuesta $O(N\log N)$ con la FFT.

**Interpretación en sonido.** El índice $k$ corresponde a la frecuencia $kf_s/N$ Hz. Analíticamente $\hat h(f)\approx\exp\!\big(-\tfrac12(2\pi f\sigma/f_s)^2\big)$: vale $1$ en los graves, $10^{-2}$ hacia los 1000 Hz y $\sim10^{-31}$ a 4000 Hz. Suavizar es "apagar" los agudos, como oír a través de una pared. Para Tikhonov, $f_k(\lambda)$ funciona como un ecualizador: devuelve los agudos hasta donde $\hat h_k^2\gg\lambda$, es decir, mientras la señal sobresale del ruido, y luego los corta. La pseudoinversa divide entre $\hat h_k\approx10^{-31}$ (en la práctica, entre redondeo: 901 de los $\hat h_k$ calculados son exactamente $0$ y 2757 son negativos). Es la misma meseta de redondeo de la sección 5.1.

**Comprobación.** `tests/test_aplicaciones.py` forma la circulante con `scipy.linalg.circulant` para $N=48$ y verifica que los valores singulares son $|\hat h_k|$ y que Tikhonov por Fourier coincide con la SVD y con el sistema aumentado (QR).

**Limitación del modelo.** La convolución circular "pega" el final del audio con su inicio. Afecta solo a unas $3\sigma\approx11$ muestras (1.5 ms) en cada extremo y no se nota al escuchar. Un modelo sin esta periodicidad necesitaría una matriz de Toeplitz y otra forma de calcular la SVD.

## 10. Ampliación: penalización de diferencias $\lambda\lVert Lx\rVert^2$

### 10.1 Motivación

Tikhonov clásico penaliza $\lVert x\rVert^2$: prefiere soluciones *pequeñas*. Pero lo que delata al ruido amplificado no es su tamaño sino sus *oscilaciones*: saltos grandes entre puntos vecinos. Penalizar las diferencias

$$
x_\lambda=\arg\min_x\Big\{\lVert Ax-b\rVert^2+\lambda\lVert Lx\rVert^2\Big\},\qquad
L_1=\begin{pmatrix}-1&1&&\\&-1&1&\\&&\ddots&\ddots\\&&&-1&1\end{pmatrix}\in\mathbb R^{(n-1)\times n},
$$

castiga la rugosidad y deja libres las constantes ($L_1\mathbf 1=0$). Con $L_2$ (segundas diferencias, $(L_2x)_j=x_j-2x_{j+1}+x_{j+2}$) se dejan libres las rectas. La penalización codifica qué consideramos una solución "simple".

### 10.2 Unicidad

**Proposición 5.** El problema tiene solución única para todo $\lambda>0$ si y solo si $\ker A\cap\ker L=\{0\}$.

*Demostración.* La hessiana de la función objetivo es $2(A^TA+\lambda L^TL)$ y

$$
z^T(A^TA+\lambda L^TL)z=\lVert Az\rVert^2+\lambda\lVert Lz\rVert^2\ \ge 0,
$$

con igualdad si y solo si $Az=0$ y $Lz=0$, es decir, $z\in\ker A\cap\ker L$. Si la intersección es trivial, la matriz es definida positiva y hay un único mínimo, solución de $(A^TA+\lambda L^TL)x=A^Tb$. Si hay $z\neq0$ en la intersección, $x+tz$ tiene el mismo valor objetivo para todo $t$ y el mínimo no es único. $\square$

En nuestro problema $\ker L_1=\operatorname{span}\{\mathbf 1\}$ y $A\mathbf 1=\mathbf 1\neq0$ (las filas suman 1), así que la condición se cumple aunque $A$ fuera singular. Para $L_2$, $\ker L_2=\operatorname{span}\{\mathbf 1,t\}$ y $A$ no anula ninguna recta. `test_nucleos_disjuntos` lo comprueba y muestra un contraejemplo: si $A$ anulara las constantes (una $A$ "centrada"), con $L_1$ la solución no sería única.

### 10.3 Fórmula espectral: la SVD generalizada

Con $L\neq I$ la SVD de $A$ ya no diagonaliza el problema, porque $L^TL$ no es la identidad en la base $V$. La herramienta adecuada es la **SVD generalizada** (GSVD, Paige–Saunders [2], [3]) del par $(A,L)$:

$$
A=U\,C\,X^{-1},\qquad L=V\,S\,X^{-1},\qquad C=\operatorname{diag}(c_i),\ S=\operatorname{diag}(s_i),\ c_i^2+s_i^2=1,
$$

con $U$, $V$ de columnas ortonormales y $X$ invertible (no ortogonal). Los cocientes $\gamma_i=c_i/s_i$ son los *valores singulares generalizados*; $\gamma_i=\infty$ en las direcciones de $\ker L$ ($s_i=0$).

**Cómo la calculamos (sin bibliotecas externas).** Sea $\begin{pmatrix}A\\L\end{pmatrix}=QR$ una QR reducida ($R$ invertible por la Proposición 5) y $Q=\begin{pmatrix}Q_1\\Q_2\end{pmatrix}$. Como $Q^TQ=I$, $Q_1^TQ_1+Q_2^TQ_2=I$. Tomamos la SVD $Q_1=U\,C\,W^T$; entonces $W^TQ_2^TQ_2W=I-C^2$ es diagonal, así que las columnas de $Q_2W$ son ortogonales con normas $s_i=\sqrt{1-c_i^2}$, y $Q_2W=VS$. Con $X=R^{-1}W$ queda $A=Q_1R=UC\,W^TR=UCX^{-1}$ y $L=Q_2R=VSX^{-1}$. (En el código $s_i$ se calcula como norma de las columnas de $Q_2W$, porque $\sqrt{1-c_i^2}$ pierde precisión cuando $c_i\approx1$.)

**Teorema 3 (fórmula espectral generalizada).** Si $\ker A\cap\ker L=\{0\}$,

$$
x_\lambda=\sum_{i=1}^{n}\frac{c_i}{c_i^2+\lambda s_i^2}\,(u_i^Tb)\,x_i
=\sum_{i=1}^{n}f_i(\lambda)\,\frac{u_i^Tb}{c_i}\,x_i,\qquad
f_i(\lambda)=\frac{c_i^2}{c_i^2+\lambda s_i^2}=\frac{\gamma_i^2}{\gamma_i^2+\lambda},
$$

donde $x_i$ son las columnas de $X$.

*Demostración.* Sustituyendo la GSVD en las ecuaciones normales, $A^TA+\lambda L^TL=X^{-T}(C^2+\lambda S^2)X^{-1}$ y $A^Tb=X^{-T}CU^Tb$. Como $X^{-T}$ es invertible, el sistema equivale a $(C^2+\lambda S^2)X^{-1}x=CU^Tb$, diagonal: $X^{-1}x=\operatorname{diag}\big(c_i/(c_i^2+\lambda s_i^2)\big)U^Tb$. Multiplicando por $X$ se obtiene la fórmula. La matriz diagonal es invertible porque $c_i^2+\lambda s_i^2>0$ (si $c_i=s_i=0$, la columna $i$ de $\begin{pmatrix}A\\L\end{pmatrix}X$ sería nula, contra la invertibilidad de $R$). $\square$

Con $L=I$ se tiene $c_i=\sigma_i/\sqrt{1+\sigma_i^2}$, $s_i=1/\sqrt{1+\sigma_i^2}$, $\gamma_i=\sigma_i$ y $X=VS$ (las columnas de $V$ escaladas por $s_i$; al sustituir, los factores $s_i$ se cancelan), y se recupera el Teorema 2 (`test_con_L_identidad_es_tikhonov_clasico`). Los factores de filtro siguen la misma lógica: $f_i\approx1$ donde $\gamma_i^2\gg\lambda$, $f_i\approx\gamma_i^2/\lambda$ donde $\gamma_i^2\ll\lambda$, y $f_i=1$ exactamente en $\ker L$: esas direcciones **nunca se penalizan**. El residuo es $\sum_i\big(\lambda s_i^2/(c_i^2+\lambda s_i^2)\big)^2(u_i^Tb)^2$, no decreciente en $\lambda$, así que la discrepancia se aplica sin cambios. Límites: $\lambda\to0^+$ da $A^\dagger b$; $\lambda\to\infty$ da la proyección sobre $\ker L$ que mejor ajusta los datos, que para $L_1$ es la constante $\operatorname{media}(b)$ porque $A\mathbf 1=\mathbf 1$ (`test_limite_lambda_grande_es_la_media`).

**Comprobaciones.** La GSVD calculada reconstruye $A$ y $L$ con error $10^{-14}$; la fórmula espectral coincide con el sistema aumentado $\min\lVert\begin{pmatrix}A\\\sqrt\lambda L\end{pmatrix}x-\begin{pmatrix}b\\0\end{pmatrix}\rVert$ resuelto por QR hasta $10^{-7}$ relativo en toda la malla; el gradiente se anula en $x_\lambda$ (`tests/test_general.py`, 14 pruebas).

### 10.4 Qué cambia y qué no (resultados en la sección 7 de los resultados)

- **Para la señal con pulso, nada.** Con las tres penalizaciones el oráculo da el mismo error ($0.1008$ para $\eta=10^{-3}$) y la discrepancia, errores dentro del 4 % entre sí. Los factores de filtro muestran por qué (figura 17): los tres filtros dejan pasar las mismas $\approx30$ direcciones, porque el corte lo impone la caída de los $\sigma_i$ de $A$ frente al ruido, no la penalización. Un salto no es "simple" para ninguna penalización cuadrática: ni $\lVert x\rVert^2$, ni $\lVert L_1x\rVert^2$, ni $\lVert L_2x\rVert^2$ lo favorecen. Recuperar bordes exige una penalización no cuadrática (variación total, $\lVert L_1x\rVert_1$), que sale del álgebra lineal y de este proyecto.
- **Para una señal suave, mucho.** Con $x(t)=0.5+0.4\sin(2\pi t)+0.2t$ (sin saltos y distinta de cero en los extremos), el error pasa de $1.0\,\%$ ($L=I$) a $0.44\,\%$ ($L_1$) y $0.08\,\%$ ($L_2$) con $\eta=10^{-3}$; en 50 semillas, de $0.95\,\%$ a $0.09\,\%$. La penalización $\lVert x\rVert^2$ empuja hacia cero una señal que no es pequeña; $\lVert L_2x\rVert^2$ solo castiga la curvatura, que en esta señal es pequeña de verdad.

**Conclusión de la ampliación.** La penalización es un *modelo a priori* de la señal. Acierta cuando describe lo que la señal es (suave) y no ayuda cuando no (un pulso). El límite del pulso no está en el método sino en lo que la observación dejó en los datos.

## 11. Preguntas de defensa, con nuestras respuestas

Las seis de la sección 11 del artículo [1] y otras que anticipamos.

1. **¿Cómo puede una solución de residuo cero ser mala?** El residuo mide el ajuste a $b$, que contiene ruido. Ajustar exactamente obliga a reproducir el ruido, y en direcciones con $\sigma_i$ pequeño eso exige coeficientes $u_i^Te/\sigma_i$ enormes. Ejemplo: $(1,11)$ con residuo 0 y error 10.
2. **¿Qué dirección del ejemplo diagonal amplifica más el ruido?** La segunda ($v_2=e_2$), porque $\sigma_2=10^{-3}$: el factor es $1/\sigma_2=1000$.
3. **¿Pseudoinversa vs. truncamiento numérico?** $A^\dagger$ divide entre todos los $\sigma_i>0$. `pinv` usa una tolerancia y trata como cero los $\sigma_i$ por debajo (conserva 89 de 128). Es una SVD truncada con un umbral pensado para el redondeo, no para el ruido de los datos, que es $10^{10}$ veces mayor; por eso aún da errores de $10^{8}$–$10^{10}$.
4. **¿Por qué $\lambda>0$ asegura unicidad?** Porque $z^T(A^TA+\lambda I)z=\lVert Az\rVert^2+\lambda\lVert z\rVert^2>0$ para $z\neq0$, aunque $A$ tenga núcleo. Con $\lambda\lVert Lx\rVert^2$ hace falta además $\ker A\cap\ker L=\{0\}$ (Proposición 5).
5. **¿Qué usa el oráculo que no se tiene en la realidad?** La señal verdadera $x_{\mathrm{true}}$, para calcular el error. La discrepancia usa solo $A$, $b$ y $\delta$; y $\delta$ se puede estimar de los datos (sección 8.4).
6. **¿Qué se gana y qué se pierde al aumentar $\lambda$?** Se gana estabilidad (el ruido propagado baja; cota $\lVert e\rVert/2\sqrt\lambda$); se pierde fidelidad: más sesgo, mayor residuo y la señal se suaviza (los bordes del pulso se redondean).
7. **¿Por qué no usan las ecuaciones normales?** Al formar $A^TA$ se pierde todo $\sigma_i<\sqrt{\varepsilon}\,\sigma_1\approx10^{-8}$: el rango numérico cae de 89 a 63. Con la SVD de $A$ o con QR no se pierde nada.
8. **¿Es $A$ singular?** En aritmética exacta no (núcleo gaussiano definido positivo, filas normalizadas). Numéricamente sí: sus $\sigma_i$ pequeños están bajo $\varepsilon\sigma_1$, y el $\kappa_2\approx2\times10^{16}$ que devuelve el cálculo es un artefacto del redondeo, no el valor verdadero.
9. **¿Qué significa que el residuo sea menor que $\lVert e\rVert$?** Que la solución explica parte del ruido como si fuera señal. La propia $x_{\mathrm{true}}$ tiene residuo exactamente $\lVert e\rVert$.
10. **¿Por qué el error no mejora al bajar el ruido?** Sí mejora, pero despacio (12 % → 6.5 % al pasar de $\eta=10^{-2}$ a $10^{-10}$): los $\sigma_i$ caen como una gaussiana, así que cada década de ruido suma pocas direcciones útiles, y el 96–99 % del error está en los bordes del pulso, donde hacen falta direcciones que $A$ ya borró.
11. **¿Por qué la penalización de diferencias no arregla los bordes?** Porque un salto no es "suave" bajo ninguna penalización cuadrática: los tres filtros cortan en las mismas direcciones. Sí ayuda (10–20 veces) cuando la señal es realmente suave.
12. **¿Cómo elegirían $\lambda$ con datos reales?** Estimar $\hat\delta$ con la cola de Picard, aplicar la discrepancia con $\tau=1.1$ y, si hay duda, preferir sobrestimar el ruido: cuesta poco; subestimarlo lleva a reproducir ruido.
13. **¿Cómo escala esto a una imagen?** Con el producto de Kronecker: la SVD de $A_h\otimes A_v$ sale de dos SVD de $256\times256$ (Proposición 3) y todas las fórmulas se aplican entrada por entrada sobre $C=U_v^TBU_h$.
14. **¿Y a un audio?** La convolución circular es una matriz circulante, que la DFT diagonaliza (Proposición 4): Tikhonov es un ecualizador frecuencia por frecuencia, en $O(N\log N)$.

## Referencias

[1] O. Velásquez Castañón, "Reconstrucción de señales con ruido: pseudoinversa, valores singulares y regularización de Tikhonov," material de referencia para proyectos, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, oct. 2026.

[2] P. C. Hansen, *Rank-Deficient and Discrete Ill-Posed Problems: Numerical Aspects of Linear Inversion*. Philadelphia, PA, USA: SIAM, 1998, doi: 10.1137/1.9780898719697.

[3] P. C. Hansen, *Discrete Inverse Problems: Insight and Algorithms*. Philadelphia, PA, USA: SIAM, 2010, doi: 10.1137/1.9780898718836.

[4] S. Gazzola, P. C. Hansen y J. G. Nagy, "IR Tools: a MATLAB package of iterative regularization methods and large-scale test problems," *Numerical Algorithms*, vol. 81, pp. 773–811, 2019, doi: 10.1007/s11075-018-0570-7.

[5] SciPy, "scipy.linalg.pinv" y "scipy.linalg.lstsq," documentación de SciPy. [En línea]. Disponible: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinv.html

[6] A. N. Tikhonov, "Solution of incorrectly formulated problems and the regularization method," *Soviet Mathematics Doklady*, vol. 4, pp. 1035–1038, 1963.

[7] V. A. Morozov, "On the solution of functional equations by the method of regularization," *Soviet Mathematics Doklady*, vol. 7, pp. 414–417, 1966.

[8] G. H. Golub y C. F. Van Loan, *Matrix Computations*, 4.ª ed. Baltimore, MD, USA: Johns Hopkins University Press, 2013.

[9] P. C. Hansen, J. G. Nagy y D. P. O'Leary, *Deblurring Images: Matrices, Spectra, and Filtering*. Philadelphia, PA, USA: SIAM, 2006, doi: 10.1137/1.9780898718874.

[10] C. F. Van Loan, "The ubiquitous Kronecker product," *Journal of Computational and Applied Mathematics*, vol. 123, no. 1–2, pp. 85–100, 2000, doi: 10.1016/S0377-0427(00)00393-9.

[11] S. van der Walt *et al.*, "scikit-image: image processing in Python," *PeerJ*, vol. 2, e453, 2014, doi: 10.7717/peerj.453. (Fuente de la fotografía de prueba `chelsea`, de S. van der Walt, publicada con licencia CC0.)
