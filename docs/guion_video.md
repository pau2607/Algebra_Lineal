# Guion del video final (versión 3, entrega final)

Duración objetivo: **9:40** (tolerancia 9–11 min, incluye identificación y créditos). Cada integrante habla unos 3:10 de contenido técnico. Las diapositivas están en `docs/diapositivas/` (una por bloque, numeradas igual que aquí).

**Cambios respecto de la versión 2:** se agregan la ampliación (bloque 10, penalización de diferencias) y la estimación del nivel de ruido (dentro del bloque 6); se corrige la interpretación del bloque 8 ("el error baja despacio", no "no baja de 9 %") y la del bloque 3 (ecuaciones normales); se comprimen los bloques 1D para hacer espacio.

## Reparto y tiempos

| # | Tiempo | Quién | Bloque | Diapositiva / recurso |
|---|---|---|---|---|
| 0 | 0:00–0:15 | Paula | Portada | 1 |
| 1 | 0:15–1:00 | Paula | El problema | 2–3 (`texto_observada_eta1e-3.png`, `02_senal_y_datos.png`) |
| 2 | 1:00–1:45 | Paula | Ejemplo diagonal | 4 (`01_ejemplo_diagonal.png`) |
| 3 | 1:45–2:30 | Matías | SVD y pseudoinversa | 5–6 |
| 4 | 2:30–3:20 | Matías | Tikhonov: unicidad, filtros, estabilidad | 7–8 (`08_factores_filtro.png`) |
| 5 | 3:20–4:00 | Carlos | El espectro de nuestra A y la pinv numérica | 9 (`03_espectro_picard.png`, `06b_pseudoinversa.png`) |
| 6 | 4:00–4:45 | Carlos | Elegir λ sin conocer la verdad; estimar δ | 10–11 (`04_curvas_tikhonov.png`, `18_estimacion_delta.png`) |
| 7 | 4:45–5:20 | Carlos | Reconstrucciones 1D y robustez | 12 (`06_reconstrucciones.png`) |
| 8 | 5:20–5:50 | Paula | Interpretación: dónde está el error que queda | 13 |
| 9a | 5:50–6:35 | Paula | Demostración: imagen | 14 (`11_imagen_texto.png`, `13_imagen_ampliacion.png`) |
| 9b | 6:35–7:20 | Matías | Demostración: audio | 15 (`14_audio.png`, archivos `.wav`) |
| 10 | 7:20–8:20 | Carlos | Ampliación: penalizar diferencias | 16–17 (`15_penalizaciones_reconstrucciones.png`, `17_penalizaciones_filtros.png`) |
| 11 | 8:20–9:15 | Matías | Conclusiones y limitaciones | 18 |
| 12 | 9:15–9:40 | Carlos | Referencias y repositorio | 19 |

## Texto de cada bloque

Está escrito para decirse en voz alta (unas 120 palabras por minuto). Cambien lo que no suene a ustedes; lo importante son las cifras y el orden de las ideas.

### 0. Portada (Paula, 0:15)

Hola. Somos Paula Romero, Carlos Ledo y Matías Desmonteix, y este es nuestro proyecto de Álgebra Lineal Aplicada: reconstrucción de señales con ruido usando la pseudoinversa, los valores singulares y la regularización de Tikhonov.

### 1. El problema (Paula, 0:45)

¿Alguien puede leer esto? *(diapositiva 2: el texto desenfocado)*. Es una imagen que pasó por un desenfoque y tiene un poco de ruido. Ese es nuestro problema: medir algo que no podemos observar directamente. Lo escribimos como b = A·x + e: x es la señal verdadera, A es el proceso que la emborrona, e es el ruido, y b es lo único que tenemos.

Para estudiarlo de forma controlada inventamos una señal de 128 puntos *(diapositiva 3)*: una montaña suave y un escalón con bordes rectos. La pasamos por un suavizado gaussiano y le sumamos ruido en tres niveles: una diezmilésima, una milésima y una centésima. Como conocemos la señal original, podemos medir cuánto se equivoca cada método. Pero nos pusimos una regla: la verdad se usa solo para calificar, nunca para decidir.

### 2. Ejemplo diagonal (Paula, 0:45)

La dificultad se ve con dos coordenadas. A = diag(1, 0.001), señal (1, 1), ruido de una centésima en la segunda coordenada. El dato medido ahí es 0.011: la señal aportó 0.001 y el ruido 0.01, diez veces más.

Si invertimos, dividimos 0.011 entre 0.001 y obtenemos 11, cuando la verdad era 1. Esa solución explica los datos a la perfección, residuo cero, y se equivoca en la señal por 10. Un ruido de una centésima se multiplicó por mil. Conclusión central del proyecto: ajustar bien los datos no significa recuperar bien la señal. Las otras dos barras son los remedios que explicará Matías.

### 3. SVD y pseudoinversa (Matías, 0:45)

La herramienta para entender esto es la descomposición en valores singulares: A = UΣVᵀ. A toma cada dirección vᵢ, la convierte en uᵢ y la multiplica por σᵢ. En nuestra A, las ondas lentas tienen σ cerca de 1 y los detalles finos, σ de 10⁻¹⁴ o menos.

La pseudoinversa deshace eso dirección por dirección: x = Σ (uᵢᵀb / σᵢ) vᵢ. Demostramos que es la solución de mínimos cuadrados de norma mínima. Y por linealidad, la parte que viene del ruido es Σ (uᵢᵀe / σᵢ) vᵢ, acotada por ‖e‖/σ_min. Los coeficientes del ruido son todos del mismo tamaño, ‖e‖/√n, pero los de la señal decaen con σᵢ. Así que a partir de cierto índice los datos son ruido, y dividir ruido entre 10⁻¹⁴ es multiplicarlo por 10¹⁴.

Un aparte: no usamos las ecuaciones normales AᵀAx = Aᵀb, porque al formar AᵀA se pierde todo σ menor que la raíz de épsilon, unos 10⁻⁸. Lo medimos: el rango numérico cae de 89 a 63.

### 4. Tikhonov (Matías, 0:50)

Tikhonov cambia el objetivo: minimizar ‖Ax − b‖² más λ‖x‖². El segundo término castiga soluciones grandes. Demostramos que para todo λ > 0 la solución es única, porque zᵀ(AᵀA + λI)z = ‖Az‖² + λ‖z‖² es positivo para todo z distinto de cero, aunque A sea singular.

En la base de la SVD, la solución es la misma suma de antes, pero cada término multiplicado por un factor de filtro fᵢ = σᵢ²/(σᵢ² + λ) *(diapositiva 8)*: vale 1 donde σᵢ² es mucho mayor que λ y se apaga suavemente donde es menor. La SVD truncada hace lo mismo con un corte abrupto: fᵢ = 1 hasta k y 0 después.

Lo que se gana es estabilidad: demostramos que el ruido se amplifica a lo sumo por 1/(2√λ), sin importar cuán pequeño sea el σ más chico. Lo que se pierde es fidelidad: el error se descompone en sesgo más ruido propagado, y al subir λ el ruido baja pero el sesgo sube.

### 5. El espectro de nuestra A (Carlos, 0:40)

Este es el gráfico de Picard de nuestro problema *(diapositiva 9)*. La curva negra son los σᵢ, que caen como una gaussiana. Los puntos son los coeficientes de los datos, uᵢᵀb. Siguen a la señal hasta que chocan con el nivel del ruido, y de ahí en adelante quedan planos: son solo ruido. Para η = 10⁻³ eso ocurre cerca de la dirección 37.

A es numéricamente singular: desde la dirección 98 los σ calculados son redondeo. La pseudoinversa de SciPy ya descarta direcciones, conserva 89 de 128, pero su umbral está pensado para el redondeo, no para el ruido, que es 10¹⁰ veces mayor. Resultado: error relativo de 10⁸ a 10¹⁰, con un residuo *menor* que el ruido. Reproduce el ruido perfectamente y la señal, nada.

### 6. Elegir λ sin conocer la verdad (Carlos, 0:45)

Estas son las curvas de error, residuo y norma frente a λ *(diapositiva 10)*. En nuestro experimento podríamos elegir el λ de menor error, el oráculo, pero eso usa la señal verdadera. Usamos el principio de discrepancia: la propia señal verdadera tiene residuo igual a ‖e‖, así que no tiene sentido ajustar mejor. Tomamos el mayor λ cuyo residuo no supere 1.1 veces δ. Solo usa A, b y δ; demostramos que el residuo crece con λ, así que la regla está bien definida.

¿Y δ? *(diapositiva 11)*. En datos reales no se conoce. Pero la cola del gráfico de Picard es solo ruido, así que su media cuadrática lo estima: en 50 semillas, δ estimado queda entre 0.89 y 1.15 veces el real, y el error final es el mismo. Toda la cadena corre sin conocer ni el ruido ni la verdad.

### 7. Reconstrucciones y robustez (Carlos, 0:35)

Estas son las reconstrucciones con λ elegido por discrepancia *(diapositiva 12)*: de un error de 10⁸ pasamos a entre 9 y 12 %. La montaña sale casi perfecta; el escalón, con los bordes redondeados y pequeñas ondulaciones.

Repetimos con 50 ruidos distintos: la discrepancia queda, en mediana, entre 0.5 y 14 % por encima del oráculo. Y probamos equivocarnos en δ: sobrestimar el ruido al doble cuesta entre 6 y 27 % más de error; subestimarlo a la mitad es grave, la regla intenta ajustar más de lo que el ruido permite y termina reproduciéndolo.

### 8. Dónde está el error que queda (Paula, 0:30)

¿Por qué queda un 10 %? Lo medimos: entre el 96 y el 99 % del error está en los bordes del escalón, 32 puntos de 128. Y es sesgo, no ruido: separando los dos términos, el sesgo aporta 0.108 y el ruido 0.007. Los bordes necesitan detalles finos que el suavizado ya destruyó.

Bajar el ruido ayuda, pero despacio: de η = 10⁻² a 10⁻¹⁰ el error pasa de 12 a 6.5 %, porque los σ caen como una gaussiana y cada década de ruido suma pocas direcciones. La regularización no recupera lo que la observación borró.

### 9a. Demostración: imagen (Paula, 0:45)

Volvamos a la imagen del inicio *(diapositiva 14)*. Una imagen de 256 × 256 desenfocada en las dos direcciones: B = A·X·Aᵀ + E. La matriz completa tendría 65 536 × 65 536, 34 gigabytes. No hace falta formarla: demostramos que la SVD del producto de Kronecker sale de las dos SVD de A, y todas las fórmulas se aplican entrada por entrada.

Miren en orden: la observada, la pseudoinversa, que es ruido puro, y Tikhonov con λ elegido por discrepancia, sin mirar la original: el texto vuelve a leerse, de 16.6 a 20.1 decibeles. Lo que no vuelve son los bordes más finos de las letras, igual que el escalón en una dimensión.

### 9b. Demostración: audio (Matías, 0:45)

Y ahora un sonido. *(Reproducir ~5 s de cada archivo: original, observado, Tikhonov.)* El suavizado apaga los agudos: se oye el tono de cada nota pero sin brillo. Aquí la matriz es circulante, y la transformada de Fourier la diagonaliza: los valores singulares son los coeficientes de Fourier del núcleo. Tikhonov queda como un ecualizador frecuencia por frecuencia, en O(N log N).

En el espectro *(diapositiva 15)* se ve que devuelve los armónicos hasta cerca de mil hertz, donde el núcleo todavía sobresale del ruido, y corta más arriba. La pseudoinversa, en cambio, es ruido blanco de amplitud 10¹³. El archivo está en el repositorio con el volumen reducido, para no dañar parlantes.

### 10. Ampliación: penalizar diferencias (Carlos, 1:00)

Para la entrega final probamos cambiar la penalización: en vez de castigar el tamaño de x, castigar sus saltos, λ‖Lx‖², con L de primeras o segundas diferencias. La unicidad ahora exige que los núcleos de A y de L no se crucen; lo demostramos y se cumple, porque L anula las constantes y A no. Y la fórmula espectral necesita la SVD generalizada del par (A, L), que calculamos con una QR y una SVD; los factores de filtro quedan γᵢ²/(γᵢ² + λ).

Resultado honesto *(diapositiva 16)*: para el escalón no cambia nada. Las tres penalizaciones dan el mismo error, porque los tres filtros dejan pasar las mismas 30 direcciones *(diapositiva 17)*. Un salto no es "suave" para ninguna penalización cuadrática. Pero para una señal realmente suave, el error baja diez veces: de 0.95 a 0.09 %. La penalización es un modelo a priori de la señal: ayuda cuando describe lo que la señal es. Para recuperar bordes haría falta un prior no cuadrático, como la variación total, que ya no es un filtro lineal.

### 11. Conclusiones y limitaciones (Matías, 0:55)

Tres conclusiones. Primera: una solución que ajusta los datos a la perfección puede ser inservible, porque ajustar es reproducir el ruido en las direcciones que A aplasta. Segunda: truncar o penalizar controla esa amplificación a costa de sesgo, y el parámetro se puede elegir sin conocer la verdad, estimando el ruido de los propios datos. Tercera: hay un límite que ningún método cruza. Lo que la observación destruyó, los bordes del escalón, los trazos finos de las letras, los agudos, no está en los datos.

Limitaciones: el ruido es sintético e isotrópico; la imagen y el audio usan una sola semilla; los desenfoques reales no siempre son separables; y probamos solo penalizaciones cuadráticas, que son filtros lineales. Y un buen resultado con una señal inventada no garantiza el mismo desempeño con datos reales.

### 12. Cierre (Carlos, 0:25)

Todo el código, las 61 pruebas automáticas, las figuras y los documentos están en el repositorio, con la etiqueta entrega-final; con `ejecutar_todo.py` se regeneran todos los resultados en menos de un minuto. Referencias principales: el material del curso, los libros de Hansen sobre problemas inversos discretos, y Golub y Van Loan. Gracias.

## Notas de producción

- Los audios están a 8000 Hz: reproducirlos desde el archivo, no grabarlos con micrófono. Avisar antes de la pseudoinversa o solo mostrarla.
- Decir explícitamente que la imagen y el audio se reconstruyen sin mirar la original: es lo que haría alguien con datos reales.
- Ensayar con cronómetro. Si sobra tiempo, el bloque 10 puede crecer; si falta, recortar el bloque 3 (dejar solo la fórmula y la amplificación).
- En pantalla, las fórmulas deben leerse a tamaño de video: usar las diapositivas, no el PDF del desarrollo.
