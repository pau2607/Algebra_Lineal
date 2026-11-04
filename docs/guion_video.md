# Guion del video final (versión preliminar)

Duración objetivo: **10:00** (tolerancia 9–11 min, incluye identificación y créditos).
Distribución sugerida por el artículo: 1 min problema · 3 min explicación espectral y Tikhonov · 4 min reconstrucciones · 2 min conclusiones.
Cada integrante explica una parte **sustantiva**; los nombres o créditos no cuentan como intervención.

> Reemplazar *Integrante A–D* por los nombres reales. La asignación coincide con la propuesta de `docs/organizacion.md`, pero el grupo puede cambiarla.

| # | Tiempo | Quién | Contenido | Recurso en pantalla |
|---|---|---|---|---|
| 0 | 0:00–0:20 | A | Tema, título, curso e integrantes | Portada |
| 1 | 0:20–1:20 | A | **Problema.** $b=Ax_{\mathrm{true}}+e$; qué es $A$ (suavizado gaussiano, $s=0.03$), la señal (zona suave + pulso), el ruido $\eta$. Objetivo y alcance. Pregunta guía: ¿cómo reconstruir sin amplificar lo que no sabemos? | `02_senal_y_datos.png` |
| 2 | 1:20–2:20 | A | **Ejemplo diagonal.** Cálculo a mano: $x^\dagger=(1,11)$ con residuo 0 y error 10; $1/\sigma_2=1000$. "Residuo pequeño no certifica buena reconstrucción". | `01_ejemplo_diagonal.png` + pizarra |
| 3 | 2:20–3:20 | B | **SVD y pseudoinversa.** Fórmula $x^\dagger=\sum (u_i^Tb/\sigma_i)v_i$, por qué es la de mínima norma, término $A^\dagger e$ y cota $\lVert e\rVert/\sigma_r$. Por qué no usamos ecuaciones normales ($\kappa^2$). | Desarrollo matemático §3–4 |
| 4 | 3:20–4:20 | B | **Tikhonov.** Función objetivo, unicidad ($z^T(A^TA+\lambda I)z>0$), fórmula espectral y factores de filtro $f_i$; comparación con el corte de la TSVD. | §6 + `08_factores_filtro.png` |
| 5 | 4:20–5:00 | C | **Estabilidad y sesgo.** Cota $\lVert e\rVert/2\sqrt\lambda$ deducida de $(s-\sqrt\lambda)^2\ge0$; descomposición sesgo + ruido; límites $\lambda\to0^+$ y $\lambda\to\infty$. | §7 |
| 6 | 5:00–6:00 | C | **Espectro de nuestro $A$.** Gráfico de Picard: dónde los coeficientes se vuelven ruido; rango numérico 89 de la `pinv` de SciPy y su umbral. Resultado de la pinv: error $10^8$–$10^{10}$ con residuo menor que el ruido. | `03_espectro_picard.png`, `06b_pseudoinversa.png` |
| 7 | 6:00–7:00 | D | **Curvas de error, residuo y norma** frente a $\lambda/\sigma_1^2$ y frente a $k$. Oráculo vs. discrepancia: qué usa cada uno. | `04_curvas_tikhonov.png`, `05_curvas_tsvd.png` |
| 8 | 7:00–8:00 | D | **Reconstrucciones** elegidas por discrepancia para los tres $\eta$; tabla resumen; robustez con 50 semillas (mediana a 0.5–14 % del oráculo); qué pasa si $\delta$ está mal estimado. | `06_reconstrucciones.png`, tablas 02 y 03 |
| 9 | 8:00–9:00 | B / C | **Interpretación.** El error no baja de ≈9 % porque domina el sesgo en los bordes del pulso (96–99 % del error); sin ruido el mínimo es 5.4 %. Variación de $s$. | `09_variacion_s.png`, tabla 03 sin ruido |
| 10 | 9:00–9:40 | A / D | **Conclusiones y limitaciones.** Truncamiento y Tikhonov controlan la amplificación a costa de sesgo; la selección del parámetro necesita una buena estimación del ruido; la regularización no recupera lo perdido; señal sintética ≠ datos reales. | Diapositiva final |
| 11 | 9:40–10:00 | Todos | Referencias principales y enlace al repositorio | Créditos |

## Pendiente para la versión final del guion

- [ ] Ajustar tiempos después del primer ensayo cronometrado.
- [ ] Decidir si se incluye la ampliación (penalización de diferencias o curva L) y en qué bloque.
- [ ] Preparar diapositivas con fórmulas legibles (tamaño mínimo de letra) y las figuras exportadas en alta resolución.
- [ ] Verificar que cada integrante tenga al menos ~2 minutos de contenido técnico.
