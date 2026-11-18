# Guion del video

**Duración total:** 10 minutos

**Documentos de apoyo:** las secciones citadas como "desarrollo, sección N" son de [`desarrollo_matematico.md`](desarrollo_matematico.md); las tablas y secciones de resultados son de [`Interpretacion_resultados.md`](Interpretacion_resultados.md). Los experimentos A, B y C corresponden a los Exp. 2, 3 y 4 del código y de las figuras.

## Resumen de tiempos e intervenciones


| Integrante | Bloques | Contenido | Tiempo | Duración |
|---|---|---|---|---|
| Yanira Sanchez | 1 y 2 | Problema, formulación y reducción a traza | 0:00 – 2:30 | 2 min 30 s |
| Alyssa Trujillo | 3 y 9 | Solución por SVD, unicidad y rotaciones propias; limitaciones y cierre | 2:30 – 4:00 y 9:00 – 10:00 | 2 min 30 s |
| Daniel Seas | 4, 5 y 6a | Implementación, experimento A (ruido) y diseño del experimento B (reflexión) | 4:00 – 6:30 | 2 min 30 s |
| Clara Palma | 6b, 7 y 8 | Resultado del experimento B, experimento C (puntos casi colineales) e interpretación de los resultados | 6:30 – 9:00 | 2 min 30 s |

---


## Bloque 1. Presentación del problema

- **Tiempo aproximado:** 0:00 – 1:00 (1 min)
- **Presenta:** Yanira Sanchez
- **Apoyo visual:** título del proyecto y figura de los puntos antes del alineamiento.

**Secuencia de contenidos**
1. Presentación del tema.
2. Ejemplo de dos nubes de puntos que representan la misma forma.
3. Pregunta del proyecto.

**Intervención prevista**

> Hola, soy Yanira Sanchez. Junto con Clara Palma, Daniel Seas y Alyssa Trujillo presentamos nuestro proyecto de Álgebra Lineal Aplicada: "Recuperación de transformaciones rígidas bidimensionales bajo ruido mediante Procrustes ortogonal".
>
> Imaginen que tenemos dos nubes de puntos, como las de la pantalla. A simple vista parecen distintas, pero en realidad representan la misma forma: una de ellas está rotada y trasladada respecto de la otra. Esto pasa mucho en la práctica, por ejemplo cuando dos cámaras observan el mismo objeto desde posiciones diferentes, o cuando comparamos formas medidas en distintos sistemas de referencia.
>
> Nuestra pregunta es doble. Primero: ¿cómo recuperar esa rotación y esa traslación a partir de los puntos correspondientes? Y segundo: ¿qué tan confiable es esa respuesta cuando los datos tienen ruido o cuando los puntos casi no nos dan información en alguna dirección?
>
> Nuestro objetivo es demostrar la solución con la descomposición en valores singulares, implementarla en Python y contrastar el error de ajuste con el error real de recuperación de la transformación.

---

## Bloque 2. Formulación y reducción a traza

- **Tiempo aproximado:** 1:00 – 2:30 (1 min 30 s)
- **Presenta:** Yanira Sanchez
- **Apoyo visual:** notación y desarrollo de las secciones 1 y 2 del desarrollo matemático (centrado y reducción a traza), línea por línea.

**Secuencia de contenidos**
1. Representación de los puntos por filas en X e Y.
2. Matrices ortogonales, O(2) y SO(2).
3. Planteamiento del problema de minimización.
4. Expansión de la norma y uso de la ortogonalidad.
5. Equivalencia con maximizar una traza.

**Intervención prevista**

> Empecemos con el problema matemático. Guardamos los puntos en dos matrices, X e Y: cada fila es un punto del plano, y la fila i de X corresponde a la fila i de Y.
>
> Suponemos que los datos cumplen Y igual a X por Q cero, más una traslación b cero, más un ruido E. No conocemos ni Q cero ni b cero: queremos estimarlos.
>
> Q debe ser ortogonal, es decir, Q transpuesta por Q igual a la identidad, para que conserve las distancias y no deforme la figura. Su determinante solo puede ser más uno, una rotación, o menos uno, una reflexión. A todas ellas las llamamos O(2), y a las que solo rotan, SO(2).
>
> El problema es encontrar Q y b que hagan a XQ más b lo más parecida posible a Y, minimizando la norma de Frobenius de la diferencia al cuadrado. La traslación se elimina primero: si centramos ambas nubes restando su promedio, la b óptima sale de los centroides y solo queda encontrar la Q que acerque Xc Q a Yc.
>
> Para resolverlo, escribimos esa norma como una traza y la expandimos. Aparecen tres partes: la norma de X, la norma de Y y un término cruzado. El paso clave está en la norma de XQ: como Q por Q transpuesta es la identidad, la traza se simplifica y queda solo la norma de X. Es decir, rotar o reflejar no cambia el tamaño de la figura.
>
> Por eso, lo único que depende de Q es el término menos dos por la traza de Q transpuesta por M, donde M es Xc transpuesta por Yc, es decir, se calcula con las nubes centradas. Así que minimizar la distancia es lo mismo que maximizar esa traza.

---

## Bloque 3. Solución por SVD, unicidad y rotaciones propias

- **Tiempo aproximado:** 2:30 – 4:00 (1 min 30 s)
- **Presenta:** Alyssa Trujillo
- **Apoyo visual:** secciones 3, 4 y 5 del desarrollo matemático (solución por SVD, unicidad y rotaciones propias).

**Secuencia de contenidos**
1. Descomposición en valores singulares de M.
2. Matriz auxiliar Z y cota de la traza.
3. Solución Q = U·Vᵀ.
4. Condición de unicidad: en el plano, dos soluciones cuando σ₂ = 0.
5. Restricción a rotaciones con la matriz D y costo de 4σ₂.

**Intervención prevista**

> Ahora tenemos que maximizar la traza de Q transpuesta por M, y para eso usamos la descomposición en valores singulares: M es igual a U por Sigma por V transpuesta, donde Sigma tiene los valores singulares sigma uno y sigma dos.
>
> Definimos una matriz auxiliar Z, igual a U transpuesta por Q por V, que también es ortogonal. Con la propiedad cíclica de la traza, todo se reduce a sigma uno por z uno uno, más sigma dos por z dos dos. Como Z es ortogonal, sus elementos diagonales son como máximo uno, así que la traza nunca supera sigma uno más sigma dos. Ese máximo se alcanza cuando Z es la identidad, y de ahí despejamos la solución: Q igual a U por V transpuesta.
>
> Sobre la unicidad: si los dos valores singulares son positivos, la solución es única. Pero si sigma dos es cero, hay una dirección que los datos no ven, y en el plano aparecen exactamente dos soluciones con el mismo costo.
>
> Por último, a veces la Q óptima resulta ser una reflexión y nosotros queremos solo rotaciones. En ese caso cambiamos el signo de la última dirección con una matriz D. La traza baja de sigma uno más sigma dos a sigma uno menos sigma dos, y por eso el costo de prohibir la reflexión es exactamente cuatro veces sigma dos.

---

## Bloque 4. Implementación

- **Tiempo aproximado:** 4:00 – 5:00 (1 min)
- **Presenta:** Daniel Seas
- **Apoyo visual:** tabla de correspondencia entre la consigna y el código (README, "Implementación y experimentos desarrollados", punto 1), tabla de la sección 7 del desarrollo matemático y la función `procrustes()`.

**Secuencia de contenidos**
1. Organización del código por secciones.
2. Los cinco pasos de la función `procrustes()`.
3. Semilla fija, ruido exacto y repeticiones.

**Intervención prevista**

> Ahora pasamos de la teoría al código. Implementamos todo en Python con NumPy, en un archivo organizado por secciones, donde cada una cumple un paso del experimento obligatorio, como muestra la tabla.
>
> El núcleo es la función procrustes, que sigue cinco pasos. Primero, centra los puntos restando el promedio, para eliminar la traslación. Segundo, forma la matriz M como Xc transpuesta por Yc. Tercero, calcula su descomposición en valores singulares. Cuarto, si queremos solo rotaciones, corrige el signo con la matriz D. Y quinto, recupera la traslación a partir de los centroides.
>
> Cuidamos dos detalles para que los resultados sean confiables. Usamos una semilla fija, así que cualquiera puede repetir exactamente nuestros experimentos. Y generamos el ruido con un tamaño exacto, proporcional al tamaño de la nube de puntos, para comparar niveles de ruido de forma justa. Además, cada experimento se repite cuatrocientas veces.

---

## Bloque 5. Experimento A: efecto del ruido

- **Tiempo aproximado:** 5:00 – 6:00 (1 min)
- **Presenta:** Daniel Seas
- **Apoyo visual:** tabla 3.2 de la interpretación y `fig2_errores_vs_ruido.png` (panel izquierdo).

**Secuencia de contenidos**
1. Diseño del experimento.
2. Resultado sin ruido.
3. Crecimiento del error con el ruido.
4. Diferencia entre residuo de ajuste y error de recuperación.

**Intervención prevista**

> En el primer experimento tomamos treinta puntos, les aplicamos una rotación de cincuenta grados y una traslación conocidas, y agregamos niveles crecientes de ruido. El algoritmo debía resolver solo con rotaciones, en SO(2).
>
> Sin ruido, el algoritmo recupera la rotación exacta: el error es del orden de diez a la menos quince, es decir, precisión de máquina. Con ruido, el error crece de forma gradual y casi proporcional al ruido: con uno por ciento de ruido, el ángulo se equivoca en cerca de siete centésimas de grado; con cinco por ciento, en casi tres décimas de grado; y con veinte por ciento de ruido, en poco más de un grado.
>
> Algo importante es que medimos dos cosas por separado: el residuo de ajuste, que compara contra los datos con ruido, y el error de recuperación, que compara contra la transformación verdadera. Con una nube bien repartida, el método es estable frente a un ruido moderado; con cincuenta por ciento de ruido ya no es tan fiable, porque la dispersión entre corridas casi iguala al propio error.

---

## Bloque 6a. Experimento B: diseño del experimento de reflexión

- **Tiempo aproximado:** 6:00 – 6:30 (30 s)
- **Presenta:** Daniel Seas
- **Apoyo visual:** fila inferior de `fig1_puntos_antes_despues.png` (reflexión: antes, O(2) y SO(2)).

**Secuencia de contenidos**
1. Diseño del experimento con una reflexión verdadera.
2. Comparación entre O(2) y SO(2).

**Intervención prevista**

> En el segundo experimento la transformación verdadera es una reflexión, con determinante menos uno, y pusimos a competir dos versiones del algoritmo. Si permitimos cualquier matriz ortogonal, es decir, en O(2), el algoritmo siempre detecta correctamente la reflexión: devuelve determinante menos uno y alinea los puntos casi perfectamente. Pero si lo obligamos a usar solo rotaciones, en SO(2), no puede reproducir el espejo: devuelve una rotación válida, pero el ajuste empeora claramente y el error de la matriz Q queda fijo en dos.

---

## Bloque 6b. Experimento B: verificación del costo 4σ₂

- **Tiempo aproximado:** 6:30 – 7:00 (30 s)
- **Presenta:** Clara Palma
- **Apoyo visual:** tabla 4.2 de la interpretación.

**Secuencia de contenidos**
1. Verificación numérica del costo 4σ₂.
2. Conexión con la sección 5 del desarrollo matemático.

**Intervención prevista**

> Cuando la verdad es una reflexión y obligamos al algoritmo a usar solo rotaciones, el ajuste empeora, y la teoría nos dice exactamente cuánto: el costo sube cuatro sigma dos, donde sigma dos es el segundo valor singular de la matriz M. Al comparar la diferencia de costos entre SO(2) y O(2) con el valor de cuatro sigma dos, vemos que coinciden en todos los niveles de ruido. De hecho, en las ciento cuarenta corridas la diferencia máxima fue de apenas diez a la menos trece, que es solo redondeo de la máquina, así que la fórmula se cumple.

---

## Bloque 7. Experimento C: puntos casi colineales

- **Tiempo aproximado:** 7:00 – 8:00 (1 min)
- **Presenta:** Clara Palma
- **Apoyo visual:** tabla 5.2 de la interpretación, `fig3_casi_colineal.png` y `fig4_colineal_puntos.png`.

**Secuencia de contenidos**
1. Diseño del experimento con nubes cada vez más delgadas.
2. Comportamiento del residuo de ajuste.
3. Error de la matriz Q y reflexiones equivocadas en O(2).

**Intervención prevista**

> En el tercer experimento mantuvimos el ruido fijo y fuimos aplanando la nube de puntos. Su espesor, que llamamos delta, baja de uno, que es una nube normal, hasta una milésima, cuando ya es casi una línea recta. Lo curioso es que, aunque la nube cambia tanto, el residuo de ajuste se mantiene en tres centésimas de principio a fin, así que el ajuste parece igual de bueno en todos los casos.
> Estos tres gráficos muestran lo que pasa en realidad, y en todos la nube se aplana a medida que avanzamos hacia la derecha. En el primero vemos el error de la transformación recuperada: mientras SO(2), en verde, se mantiene en una milésima, O(2), en naranja, se dispara hasta 1,4. Sin embargo, en el segundo el residuo de ajuste sigue plano, es decir, no nos avisa de que algo anda mal. El tercero explica ese salto: con la nube más plana, O(2) empieza a elegir una reflexión en vez de la rotación, en catorce de veinte corridas, y cerca del 40% de las veces cuando repetimos 400 veces. Por eso decimos que un buen ajuste no garantiza haber encontrado la transformación correcta.
---

## Bloque 8. Interpretación de los resultados

- **Tiempo aproximado:** 8:00 – 9:00 (1 min)
- **Presenta:** Clara Palma
- **Apoyo visual:** sección 5.3 de la interpretación (puntos 3 a 7), tabla 5.2 y tabla de la sección 5.4 (fracción de reflexiones erróneas según r).

**Secuencia de contenidos**
1. Explicación del fenómeno con el valor singular σ₂.
2. Conexión con la sección 4 del desarrollo: dos soluciones posibles en el plano.
3. Por qué SO(2) sigue estable y O(2) no.
4. Criterio relativo r = σ₂(Xc)/s: espesor de la nube comparado con el ruido (simulación auxiliar, sección 5.4).
5. Conclusión principal: ajuste no es identificación; r como diagnóstico.

**Intervención prevista**


> Esto ocurre porque, cuando los puntos están casi en una recta, una rotación y una reflexión los dejan prácticamente en el mismo lugar, y lo único que las distingue es el grosor de la nube. Si ese grosor es tan pequeño como el ruido, termina siendo el ruido el que decide cuál elige o de dos. Ese o de dos evita este problema porque descarta la reflexión, aunque eso solo sirve si sabemos de antemano que la verdad es una rotación.
> Para medirlo definimos r, que es el grosor de la nube, es decir, su segundo valor singular, dividido entre el tamaño del ruido. Con una simulación de 400 repeticiones comprobamos que es justamente r lo que decide. Cuando r es mayor que tres, o de dos es confiable; cuando está cerca de uno, falla entre el 10% y el 20% de las veces; y cuando baja de 0.3, elige casi al azar. En conclusión, para confiar en el resultado no basta con mirar el residuo: también hay que mirar r.

---

## Bloque 9. Limitaciones y cierre

- **Tiempo aproximado:** 9:00 – 10:00 (1 min)
- **Presenta:** Alyssa Trujillo
- **Apoyo visual:** diapositiva de cierre con las conclusiones.

**Secuencia de contenidos**
1. Limitaciones del trabajo.
2. Resumen de los resultados.
3. Siguiente paso y despedida.

**Intervención prevista**

> Para cerrar, queremos ser claros sobre los límites de lo que mostramos. Trabajamos con datos sintéticos en el plano, sin estimar cambios de escala, con ruido igual en todas las direcciones y sabiendo de antemano qué punto corresponde a cuál. En problemas reales esto no siempre se cumple: puede haber varias configuraciones a la vez, correspondencias desconocidas o puntos atípicos, que son extensiones que dejamos fuera de esta etapa.
>
> En resumen: demostramos que el problema de Procrustes se reduce a maximizar una traza y que la SVD lo resuelve con Q igual a U por V transpuesta; implementamos la solución y verificamos numéricamente los resultados teóricos principales, como la fórmula del costo mínimo y el sobrecosto de cuatro sigma dos; y encontramos que, con reflexiones permitidas, la confiabilidad depende de r, el espesor de la nube comparado con el ruido.
>
> Las referencias principales están en pantalla, y todo el código, los datos y los resultados están en nuestro repositorio de GitHub, cuyo enlace aparece aquí y en el README. Muchas gracias.
