# Desarrollo matemático

**Proyecto:** Reconstrucción de señales con ruido pseudoinversa, valores singulares y regularización de Tikhonov
**Curso:** Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico
**Estado:** entrega parcial

Notación: $\lVert\cdot\rVert$ es la norma euclídea; $A^{\dagger}$, la pseudoinversa de Moore–Penrose; $\kappa_2(A)=\sigma_1/\sigma_n$, el número de condición en norma 2. En Tikhonov, $\lambda$ es **directamente** el coeficiente de la penalización (no $\lambda^2$), igual que en el artículo de referencia [1].

---

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

1. Reproducir el ejemplo diagonal de dos coordenadas y comprobar que el código coincide con los cálculos a mano (tolerancia $10^{-10}$).
2. Para cada $\eta$, comparar pseudoinversa numérica, SVD truncada y Tikhonov midiendo error relativo $\lVert x-x_{\mathrm{true}}\rVert/\lVert x_{\mathrm{true}}\rVert$, residuo $\lVert Ax-b\rVert$ y norma $\lVert x\rVert$.
3. Elegir el parámetro con una regla que no usa $x_{\mathrm{true}}$ (principio de discrepancia) y cuantificar cuánto se aleja del mejor parámetro posible (oráculo), con 50 semillas.
4. Explicar con la SVD por qué la pseudoinversa falla y qué parte del error no puede eliminarse regularizando.

### 1.4 Alcance

Incluye: señal unidimensional, penalización identidad $\lambda\lVert x\rVert^2$, SVD calculada con NumPy (no se programa desde cero), principio de discrepancia y curva L solo como gráfico exploratorio.

---

## 2. Conceptos del curso que se utilizan

| Concepto | Dónde aparece en el proyecto |
|---|---|
| Mínimos cuadrados y ecuaciones normales | Planteamiento del ajuste; por qué no conviene usar $A^TA$ |
| Núcleo, imagen y rango | Solución de mínima norma; componentes que ningún método recupera |
| Bases ortonormales y proyección ortogonal | Coordenadas en $\{u_i\}$ y $\{v_i\}$; residuo fuera de la imagen |
| Formas cuadráticas definidas positivas | Unicidad de la solución de Tikhonov |
| SVD y pseudoinversa | Fórmulas de las tres reconstrucciones, filtrado espectral |
| Número de condición | Amplificación del ruido |

---

## 3. Mínimos cuadrados y pseudoinversa

El ajuste básico es

$$
\min_{x}\ \lVert Ax-b\rVert^2 .
$$

Sus ecuaciones normales son $A^TAx=A^Tb$. No las usamos para calcular: para columnas independientes $\kappa_2(A^TA)=\kappa_2(A)^2$. En nuestro problema el $\kappa_2(A)$ calculado es del orden de $2\times10^{16}$, así que $\kappa_2(A^TA)\approx5\times10^{32}$ excede por mucho $1/\varepsilon_{\text{máq}}\approx4.5\times10^{15}$ y formar $A^TA$ destruye la información.

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

Observación. $x^\dagger$ es una elección entre todos los ajustes óptimos. Si $x_{\mathrm{true}}$ tiene componente en $\ker A$, ningún método basado solo en $b$ puede recuperarla.

---

## 4. Por qué se amplifica el ruido

Por linealidad, $x^\dagger = A^\dagger Ax_{\mathrm{true}} + A^\dagger e$, y la parte que viene del ruido es

$$
A^{\dagger}e=\sum_{i=1}^{r}\frac{u_i^Te}{\sigma_i}\,v_i ,
\qquad
\lVert A^\dagger e\rVert=\Big(\sum_i \frac{(u_i^Te)^2}{\sigma_i^2}\Big)^{1/2}\le \frac{\lVert e\rVert}{\sigma_r}.
$$

El coeficiente de $v_i$ es el coeficiente del ruido en la dirección $u_i$ multiplicado por $1/\sigma_i$. Como $U$ es ortogonal y $e$ es gaussiano isotrópico, los $u_i^Te$ tienen todos aproximadamente el mismo tamaño, $\lVert e\rVert/\sqrt n$. En cambio, para datos limpios $u_i^T(Ax_{\mathrm{true}})=\sigma_i\,v_i^Tx_{\mathrm{true}}$ decae con $\sigma_i$. Hay entonces un índice a partir del cual los coeficientes de $b$ son solo ruido y dividirlos entre $\sigma_i$ produce oscilaciones enormes. Esto es lo que muestra el gráfico de Picard (`resultados/figuras/03_espectro_picard.png`).

### 4.1 Ejemplo diagonal (cálculo a mano)

$$
A=\begin{pmatrix}1&0\\0&10^{-3}\end{pmatrix},\quad x_{\mathrm{true}}=\begin{pmatrix}1\\1\end{pmatrix},\quad e=\begin{pmatrix}0\\0.01\end{pmatrix},\quad b=\begin{pmatrix}1\\0.011\end{pmatrix}.
$$

Aquí $u_i=v_i=e_i$, $\sigma_1=1$, $\sigma_2=10^{-3}$.

- Pseudoinversa: $x^\dagger=(1/1,\ 0.011/10^{-3})=(1,11)$. Residuo $0$, error $\lVert(0,10)\rVert=10$. La dirección 2 multiplica el ruido por $1/\sigma_2=1000$: un error de $0.01$ en el dato se vuelve un error de $10$ en la solución.
- SVD truncada, $k=1$: $x_1=(1,0)$. Error $1$, residuo $\lVert(0,-0.011)\rVert=0.011$.
- Tikhonov, $\lambda=10^{-4}$: $x_\lambda=\Big(\tfrac{1}{1+10^{-4}},\ \tfrac{10^{-3}\cdot0.011}{10^{-6}+10^{-4}}\Big)=\Big(\tfrac{1}{1.0001},\tfrac{0.000011}{0.000101}\Big)\approx(0.99990,\ 0.10891)$. Error $\approx0.8911$.

Conclusión: un residuo pequeño no certifica una buena reconstrucción. El script `scripts/01_ejemplo_diagonal.py` y las pruebas `test_diagonal_*` verifican estos números.

---

## 5. SVD truncada

Para $k\le r$,

$$
x_k=\sum_{i=1}^{k}\frac{u_i^Tb}{\sigma_i}\,v_i .
$$

Se descartan las direcciones con $\sigma_i$ pequeño: mejora la estabilidad, pero se pierden las componentes verdaderas $v_i^Tx_{\mathrm{true}}$ con $i>k$. El residuo es

$$
\lVert Ax_k-b\rVert^2=\sum_{i=k+1}^{r}(u_i^Tb)^2+\lVert b_\perp\rVert^2 ,
$$

que no crece al aumentar $k$ (lo comprueba `test_residuo_tsvd_no_crece_con_k`).

### 5.1 Pseudoinversa numérica ≠ inversa sin regularizar

`scipy.linalg.pinv` descarta los $\sigma_i<\texttt{rtol}\cdot\sigma_1$ con $\texttt{rtol}=\max(m,n)\,\varepsilon$ por defecto [5]. En nuestro problema ese umbral es $\approx2.84\times10^{-14}$ y conserva 89 de 128 valores singulares. Es decir, la "pseudoinversa numérica" ya es una SVD truncada con $k=89$; lo informamos explícitamente. Para comparar, también calculamos la suma sobre los 128 $\sigma_i$ positivos ("pseudoinversa completa"). Los $\sigma_i$ con $i\gtrsim95$ forman una meseta cerca de $10^{-16}$: son errores de redondeo de la SVD, no información del operador.

---

## 6. Regularización de Tikhonov

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

La SVD truncada corresponde a $f_i\in\{0,1\}$ (corte abrupto); Tikhonov es un filtro suave (`resultados/figuras/08_factores_filtro.png`).

Calculamos $\sigma_i/(\sigma_i^2+\lambda)\cdot(u_i^Tb)$ en lugar de $f_i\cdot(u_i^Tb)/\sigma_i$ para no dividir entre $\sigma_i\approx10^{-17}$. Contrastamos con el sistema aumentado

$$
\min_x\left\lVert\begin{pmatrix}A\\ \sqrt\lambda\,I\end{pmatrix}x-\begin{pmatrix}b\\0\end{pmatrix}\right\rVert^2,
$$

cuyas ecuaciones normales son exactamente $(A^TA+\lambda I)x=A^Tb$, resuelto con QR sin formar $A^TA$. La diferencia relativa máxima en la malla es $\approx3\times10^{-10}$.

### 6.2 Residuo como función de $\lambda$

$$
\lVert Ax_\lambda-b\rVert^2=\sum_{i=1}^{r}\left(\frac{\lambda}{\sigma_i^2+\lambda}\right)^2(u_i^Tb)^2+\lVert b_\perp\rVert^2 .
$$

Cada factor $\lambda/(\sigma_i^2+\lambda)$ es creciente en $\lambda$, así que el residuo crece con $\lambda$. Esto hace que la regla de discrepancia de la sección 8 esté bien definida.

---

## 7. Estabilidad, sesgo y descomposición del error

### 7.1 Cota de estabilidad

Si los datos cambian de $b$ a $b+e$, por linealidad

$$
\lVert x_\lambda(b+e)-x_\lambda(b)\rVert=\Big\lVert\sum_i\frac{\sigma_i}{\sigma_i^2+\lambda}(u_i^Te)v_i\Big\rVert\le\max_{s\ge0}\frac{s}{s^2+\lambda}\,\lVert e\rVert=\frac{\lVert e\rVert}{2\sqrt\lambda}.
$$

El máximo se obtiene de $(s-\sqrt\lambda)^2\ge0\iff s^2+\lambda\ge2s\sqrt\lambda\iff \dfrac{s}{s^2+\lambda}\le\dfrac{1}{2\sqrt\lambda}$, con igualdad en $s=\sqrt\lambda$. Comparado con la pseudoinversa (amplificación $1/\sigma_r$), Tikhonov amplifica a lo sumo $1/(2\sqrt\lambda)$, independiente de cuán pequeño sea $\sigma_r$.

### 7.2 Descomposición sesgo + ruido

Con $b=Ax_{\mathrm{true}}+e$ y $u_i^TAx_{\mathrm{true}}=\sigma_iv_i^Tx_{\mathrm{true}}$:

$$
x_\lambda-x_{\mathrm{true}}=\underbrace{\sum_{i=1}^{r}\big(f_i(\lambda)-1\big)(v_i^Tx_{\mathrm{true}})\,v_i\;-\;P_{\ker A}\,x_{\mathrm{true}}}_{\text{sesgo (no depende del ruido)}}\;+\;\underbrace{\sum_{i=1}^{r}f_i(\lambda)\frac{u_i^Te}{\sigma_i}\,v_i}_{\text{ruido propagado}} .
$$

- Al aumentar $\lambda$: los $f_i$ bajan, el término de ruido baja (cota $\lVert e\rVert/(2\sqrt\lambda)$) y el sesgo sube.
- Al disminuir $\lambda$: lo contrario.
- $P_{\ker A}x_{\mathrm{true}}$ no depende de $\lambda$: es información perdida.

Límites: $\lambda\to0^+$ da $x_\lambda\to A^\dagger b$; $\lambda\to\infty$ da $x_\lambda\to0$ (`test_limites_de_lambda`). El mejor $\lambda$ depende de la escala de $A$, del ruido y de la señal; por eso la malla se expresa como $\lambda/\sigma_1^2\in[10^{-10},1]$.

---

## 8. Elección del parámetro sin conocer la verdad

### 8.1 Oráculo

En un experimento sintético podemos calcular el error para cada $\lambda$ de la malla y tomar el mínimo. Esto usa $x_{\mathrm{true}}$, que no existe en un problema real. Lo usamos solo como referencia para evaluar la regla que sí es aplicable.

### 8.2 Principio de discrepancia (Morozov [7])

Si se conoce una cota $\delta\approx\lVert e\rVert$, no tiene sentido ajustar los datos mejor que el nivel del ruido: un residuo menor que $\lVert e\rVert$ significa que la solución está reproduciendo ruido. La regla es:

> Elegir el mayor $\lambda$ de la malla tal que $\lVert Ax_\lambda-b\rVert\le\tau\delta$, con $\tau=1.1$.
> Para la SVD truncada: el menor $k$ tal que $\lVert Ax_k-b\rVert\le\tau\delta$.

Por la monotonía de la sección 6.2 (y 5), el conjunto de $\lambda$ viables es un intervalo inicial de la malla y la regla elige su extremo derecho. Supuestos: (i) $\delta$ es una buena estimación del ruido; (ii) los errores de modelo son menores que el ruido; (iii) el umbral es alcanzable en la malla, lo que exige $\min_\lambda\lVert Ax_\lambda-b\rVert\le\tau\delta$. Si no se cumple (iii), la función devuelve "no viable" en lugar de inventar un valor. En el código, las funciones de discrepancia no reciben $x_{\mathrm{true}}$ (lo comprueba `test_discrepancia_no_recibe_la_verdad`).

En este experimento tomamos $\delta=\lVert e\rVert$ (conocido porque generamos el ruido) y estudiamos aparte qué ocurre si $\delta$ está mal estimado ($\delta=0.5\lVert e\rVert$ y $2\lVert e\rVert$).

### 8.3 Curva L (exploratoria)

Dibujamos $\log\lVert Ax_\lambda-b\rVert$ frente a $\log\lVert x_\lambda\rVert$. Su esquina es una heurística de compromiso; en esta entrega solo se grafica, no se usa para elegir.

---

## Referencias

[1] O. Velásquez Castañón, "Reconstrucción de señales con ruido: pseudoinversa, valores singulares y regularización de Tikhonov," material de referencia para proyectos, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, oct. 2026.

[2] P. C. Hansen, *Rank-Deficient and Discrete Ill-Posed Problems: Numerical Aspects of Linear Inversion*. Philadelphia, PA, USA: SIAM, 1998, doi: 10.1137/1.9780898719697.

[3] P. C. Hansen, *Discrete Inverse Problems: Insight and Algorithms*. Philadelphia, PA, USA: SIAM, 2010, doi: 10.1137/1.9780898718836.

[4] S. Gazzola, P. C. Hansen y J. G. Nagy, "IR Tools: a MATLAB package of iterative regularization methods and large-scale test problems," *Numerical Algorithms*, vol. 81, pp. 773–811, 2019, doi: 10.1007/s11075-018-0570-7.

[5] SciPy, "scipy.linalg.pinv" y "scipy.linalg.lstsq," documentación de SciPy. [En línea]. Disponible: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinv.html

[6] A. N. Tikhonov, "Solution of incorrectly formulated problems and the regularization method," *Soviet Mathematics Doklady*, vol. 4, pp. 1035–1038, 1963.

[7] V. A. Morozov, "On the solution of functional equations by the method of regularization," *Soviet Mathematics Doklady*, vol. 7, pp. 414–417, 1966.

[8] G. H. Golub y C. F. Van Loan, *Matrix Computations*, 4.ª ed. Baltimore, MD, USA: Johns Hopkins University Press, 2013.
