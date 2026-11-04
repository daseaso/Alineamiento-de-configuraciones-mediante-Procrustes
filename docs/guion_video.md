# Guion del video

**Duración total:** 10 minutos

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

> Hola, nuestro proyecto es sobre el alineamiento de configuraciones mediante Procrustes.
>
> Imaginen que tenemos dos nubes de puntos, como las de la pantalla. A simple vista parecen distintas, pero en realidad representan la misma forma: una de ellas está rotada y trasladada respecto de la otra. Esto pasa mucho en la práctica, por ejemplo cuando dos cámaras observan el mismo objeto desde posiciones diferentes, o cuando comparamos formas medidas en distintos sistemas de referencia.
>
> Nuestra pregunta es doble. Primero: ¿cómo recuperar exactamente esa rotación y esa traslación a partir de los puntos correspondientes? Y segundo: ¿qué tan confiable es esa respuesta cuando los datos tienen ruido o cuando los puntos casi no nos dan información en alguna dirección?
>
> Vamos a responder con álgebra lineal: matrices ortogonales, la traza y la descomposición en valores singulares.

---

## Bloque 2. Formulación y reducción a traza

- **Tiempo aproximado:** 1:00 – 2:30 (1 min 30 s)
- **Presenta:** Yanira Sanchez
- **Apoyo visual:** definiciones y desarrollo del inciso 2.a, línea por línea.

**Secuencia de contenidos**
1. Representación de los puntos por filas en X e Y.
2. Matrices ortogonales, O(2) y SO(2).
3. Planteamiento del problema de minimización.
4. Expansión de la norma y uso de la ortogonalidad.
5. Equivalencia con maximizar una traza.

**Intervención prevista**

> Empecemos formulando el problema. Guardamos los puntos en dos matrices, X e Y, donde cada fila es un punto del plano, y la fila i de X corresponde a la fila i de Y.
>
> Buscamos una matriz Q ortogonal, es decir, que cumpla Q transpuesta por Q igual a la identidad. Estas matrices conservan las distancias, así que no deforman la figura. Además, su determinante solo puede ser más uno, que corresponde a una rotación, o menos uno, que corresponde a una reflexión. Al conjunto de todas ellas lo llamamos O(2), y a las que solo rotan, SO(2).
>
> El problema es encontrar la Q que haga a XQ lo más parecida posible a Y, minimizando la norma de Frobenius de XQ menos Y al cuadrado.
>
> Para resolverlo, escribimos esa norma como una traza y la expandimos. Aparecen tres partes: la norma de X, la norma de Y y un término cruzado. El paso clave está en la norma de XQ: como Q por Q transpuesta es la identidad, la traza se simplifica y queda solo la norma de X. Es decir, rotar o reflejar no cambia el tamaño de la figura.
>
> Por eso, lo único que depende de Q es el término menos dos por la traza de Q transpuesta por M, donde M es X transpuesta por Y. Así que minimizar la distancia es lo mismo que maximizar esa traza.

---

## Bloque 3. Solución por SVD, unicidad y rotaciones propias

- **Tiempo aproximado:** 2:30 – 4:00 (1 min 30 s)
- **Presenta:** Alyssa Trujillo
- **Apoyo visual:** incisos 2.b, 2.c y 2.d.

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
- **Apoyo visual:** tabla de la sección 3.1 y la función `procrustes()`.

**Secuencia de contenidos**
1. Organización del código por secciones.
2. Los cinco pasos de la función `procrustes()`.
3. Semilla fija, ruido exacto y repeticiones.

**Intervención prevista**

> Ahora pasamos de la teoría al código. Implementamos todo en Python con NumPy, en un archivo organizado por secciones, donde cada una cumple un paso del experimento obligatorio, como muestra la tabla.
>
> El núcleo es la función procrustes, que sigue cinco pasos. Primero, centra los puntos restando el promedio, para eliminar la traslación. Segundo, forma la matriz M como Xc transpuesta por Yc. Tercero, calcula su descomposición en valores singulares. Cuarto, si queremos solo rotaciones, corrige el signo con la matriz D. Y quinto, recupera la traslación a partir de los centroides.
>
> Cuidamos dos detalles para que los resultados sean confiables. Usamos una semilla fija, así que cualquiera puede repetir exactamente nuestros experimentos. Y generamos el ruido con un tamaño exacto, proporcional al tamaño de la nube de puntos, para comparar niveles de ruido de forma justa. Además, cada nivel se repite veinte veces.

---

## Bloque 5. Experimento A: efecto del ruido

- **Tiempo aproximado:** 5:00 – 6:00 (1 min)
- **Presenta:** Daniel Seas
- **Apoyo visual:** tabla 4.1.

**Secuencia de contenidos**
1. Diseño del experimento.
2. Resultado sin ruido.
3. Crecimiento del error con el ruido.
4. Diferencia entre residuo de ajuste y error de recuperación.

**Intervención prevista**

> En el primer experimento tomamos treinta puntos, les aplicamos una rotación de cincuenta grados y una traslación conocidas, y agregamos niveles crecientes de ruido. El algoritmo debía resolver solo con rotaciones, en SO(2).
>
> Sin ruido, el algoritmo recupera la rotación exacta: el error es cero. Con ruido, el error crece de forma gradual: con uno por ciento de ruido, el ángulo se equivoca en cerca de siete centésimas de grado; con cinco por ciento, en casi un tercio de grado; y aun con veinte por ciento de ruido, el error apenas supera un grado.
>
> Algo importante es que medimos dos cosas por separado: el residuo de ajuste, que compara contra los datos con ruido, y el error de recuperación, que compara contra la transformación verdadera. El método resulta muy estable frente al ruido.

---

## Bloque 6a. Experimento B: diseño del experimento de reflexión

- **Tiempo aproximado:** 6:00 – 6:30 (30 s)
- **Presenta:** Daniel Seas
- **Apoyo visual:** resultados del punto 4.2.

**Secuencia de contenidos**
1. Diseño del experimento con una reflexión verdadera.
2. Comparación entre O(2) y SO(2).

**Intervención prevista**

> En el segundo experimento la transformación verdadera es una reflexión, con determinante menos uno, y pusimos a competir dos versiones del algoritmo. Si permitimos cualquier matriz ortogonal, es decir, en O(2), el algoritmo siempre detecta correctamente la reflexión: devuelve determinante menos uno y alinea los puntos casi perfectamente. Pero si lo obligamos a usar solo rotaciones, en SO(2), no puede reproducir el espejo y el ajuste empeora.

---

## Bloque 6b. Experimento B: verificación del costo 4σ₂

- **Tiempo aproximado:** 6:30 – 7:00 (30 s)
- **Presenta:** Clara Palma
- **Apoyo visual:** resultados del punto 4.2.

**Secuencia de contenidos**
1. Verificación numérica del costo 4σ₂.
2. Conexión con el inciso 2.d.

**Intervención prevista**

> Lo interesante es cuánto empeora. La teoría del inciso 2.d dice que el costo debe subir exactamente cuatro veces sigma dos. Al calcular la diferencia entre el costo de SO(2) y el de O(2), coincidió con cuatro sigma dos con un margen de error del orden de diez a la menos trece, prácticamente cero. Así comprobamos numéricamente nuestra demostración.

---

## Bloque 7. Experimento C: puntos casi colineales

- **Tiempo aproximado:** 7:00 – 8:00 (1 min)
- **Presenta:** Clara Palma
- **Apoyo visual:** tabla 4.3.

**Secuencia de contenidos**
1. Diseño del experimento con nubes cada vez más delgadas.
2. Comportamiento del residuo de ajuste.
3. Error de la matriz Q y reflexiones equivocadas en O(2).

**Intervención prevista**

> El tercer experimento es el más revelador. Mantuvimos el ruido fijo y fuimos reduciendo el ancho de la nube de puntos, que llamamos delta, hasta que casi formaron una línea recta.
>
> Lo primero que observamos es que el residuo de ajuste se mantiene prácticamente constante. El algoritmo parece funcionar igual de bien. Pero cuando miramos la transformación recuperada, la historia es otra. Con la nube normal, el error de la matriz Q es de apenas una milésima y nunca se equivoca. Pero cuando delta baja a tres milésimas, el error sube a cero coma veinte y el algoritmo elige la reflexión equivocada en el diez por ciento de las corridas. Y cuando los puntos casi forman una línea, con delta igual a una milésima, el error se dispara a uno coma cuarenta, y se equivoca en el setenta por ciento de las veces.

---

## Bloque 8. Interpretación de los resultados

- **Tiempo aproximado:** 8:00 – 9:00 (1 min)
- **Presenta:** Clara Palma
- **Apoyo visual:** interpretación del punto 4.3.

**Secuencia de contenidos**
1. Explicación del fenómeno con el valor singular σ₂.
2. Conexión con el inciso 2.c: dos soluciones posibles en el plano.
3. Conclusión principal: ajuste no es identificación.
4. El cociente σ₂/σ₁ como diagnóstico.

**Intervención prevista**

> ¿Por qué pasa esto? Cuando los puntos casi forman una recta, el cociente entre sigma dos y sigma uno tiende a cero, y la matriz M se acerca a ser singular. Eso significa que los datos casi no dan información en la dirección perpendicular a la recta.
>
> Como vimos en el inciso 2.c, cuando sigma dos es cero, en el plano existen exactamente dos soluciones con el mismo costo: la rotación correcta y una reflexión. Si sigma dos es casi cero, las dos opciones cuestan casi lo mismo, y el algoritmo ya no puede distinguirlas. Eso explica las reflexiones equivocadas.
>
> La conclusión principal de nuestro proyecto es que un buen residuo de ajuste no garantiza haber identificado correctamente la transformación. Por eso, el cociente entre sigma dos y sigma uno funciona como un diagnóstico indispensable: si se acerca a cero, el resultado deja de ser confiable, aunque el ajuste se vea bien.

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
> En resumen: demostramos que el problema de Procrustes se reduce a maximizar una traza y que la SVD lo resuelve con Q igual a U por V transpuesta; implementamos la solución y verificamos numéricamente cada resultado teórico, incluido el costo de cuatro sigma dos; y encontramos que la estabilidad depende del segundo valor singular.
>
> Como siguiente paso, vamos a estudiar con más detalle cómo cambia Q ante perturbaciones de los datos. Muchas gracias.
