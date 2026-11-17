# Desarrollo matemático: problema de Procrustes ortogonal en el plano

## Notación y convención

- $`X, Y \in \mathbb{R}^{n\times 2}`$: cada **fila** es un punto del plano; la fila $`i`$ de $`X`$ corresponde a la fila $`i`$ de $`Y`$.
- $`\mathbf{1}\in\mathbb{R}^n`$ es el vector de unos; $`\bar x = \tfrac1n X^\top\mathbf{1}`$ y $`\bar y = \tfrac1n Y^\top\mathbf{1}`$ son los centroides (como vectores columna).
- $`X_c = X-\mathbf{1}\bar x^\top`$ e $`Y_c = Y-\mathbf{1}\bar y^\top`$ son las matrices centradas.
- $`\lVert A\rVert_F^2=\mathrm{tr}(A^\top A)`$ (norma de Frobenius).
- $`O(2)=\lbrace Q: Q^\top Q=I\rbrace`$ y $`SO(2)=\lbrace Q\in O(2):\det Q=+1\rbrace`$. Una matriz de $`O(2)`$ tiene determinante $`+1`$ (rotación) o $`-1`$ (reflexión).
- Convención por filas: el modelo es $`Y\approx XQ+\mathbf{1}b^\top`$; cada punto $`x_i^\top`$ pasa a $`x_i^\top Q`$, es decir $`x\mapsto Q^\top x`$.
- **Sentido de giro.** Por la convención por filas, con $`Q_0=R(50^\circ)`$ cada punto pasa a $`Q_0^\top x=R(-50^\circ)\,x`$: **la nube gira 50° en sentido horario** (se ve en las figuras 1 y 4). El error angular se mide con $`Q_0^\top Q`$ y no depende de esta convención. La convención por filas sigue a [[8, secc. 3]](fuentes.md#ref8).

**Problema.** Hallar $`Q\in O(2)`$ (o $`SO(2)`$) y $`b\in\mathbb{R}^2`$ que minimicen

```math
f(Q,b)=\lVert XQ+\mathbf{1}b^\top-Y\rVert_F^2 .
```

(La escala se fija en $`s=1`$; ver la [sección 6](#sec6).) Planteamiento del tema según [[8]](fuentes.md#ref8).

<a id="sec0"></a>

## 0. Supuestos

| Supuesto | Dónde se usa | Qué pasa si no se cumple |
|---|---|---|
| **S1. Correspondencia conocida:** la fila $`i`$ de $`X`$ corresponde a la fila $`i`$ de $`Y`$. | Toda la formulación: $`f`$ compara fila con fila. | El problema pasa a ser Procrustes–Wasserstein [[4]](fuentes.md#ref4), fuera del alcance. |
| **S2. Modelo rígido sin escala:** $`Y=XQ_0+\mathbf{1}b_0^\top+E`$ con $`Q_0\in O(2)`$ y $`s=1`$. | Secciones 1 a 5. | Con escala se usa la [sección 6](#sec6). |
| **S3. Ruido aditivo, independiente e isotrópico** (misma varianza en las dos direcciones). | Justifica la norma de Frobenius: con errores gaussianos isotrópicos, minimizarla equivale a máxima verosimilitud [[8, secc. 11.3]](fuentes.md#ref8). También se usa en las secciones 8 y 9. | Con ruido anisotrópico o atípicos habría que usar pesos o pérdidas robustas (fuera del alcance). |
| **S4. Nube no colineal** ($`\mathrm{rango}\,X_c=2`$). | Unicidad en $`O(2)`$ ([sección 4](#sec4)). | Si $`X_c`$ es colineal, $`\sigma_2(M)=0`$ y hay dos soluciones en $`O(2)`$. |
| **S5. Dimensión $`p=2`$.** | Las secciones 4, 5, 8 y 9 usan que hay exactamente dos valores singulares. | Las ideas se generalizan, pero las fórmulas cambian. |

---

<a id="sec1"></a>

## 1. Eliminación de la traslación (centrado)

Sea $`A=XQ-Y`$. Para $`Q`$ fija, $`f(b)=\lVert A+\mathbf{1}b^\top\rVert_F^2`$ es una función cuadrática y convexa de $`b`$ (su hessiano es $`2nI`$). Su gradiente es

```math
\frac{\partial f}{\partial b}=2\,(A+\mathbf{1}b^\top)^\top\mathbf{1}=2\,(Q^\top X^\top\mathbf{1}-Y^\top\mathbf{1}+n\,b).
```

Igualando a cero: $`b=\bar y-Q^\top\bar x`$, es decir,

```math
b^\top=\bar y^\top-\bar x^\top Q .
```

Al sustituir,

```math
XQ+\mathbf{1}b^\top-Y=(X-\mathbf{1}\bar x^\top)Q-(Y-\mathbf{1}\bar y^\top)=X_cQ-Y_c ,
```

de modo que $`Q`$ se obtiene resolviendo $`\min_Q\lVert X_cQ-Y_c\rVert_F^2`$, sin depender de $`b`$. Una vez hallada $`Q`$, la traslación es $`b=\bar y-Q^\top\bar x`$. Este centrado se trata en [[3]](fuentes.md#ref3) y [[8]](fuentes.md#ref8).

<a id="sec2"></a>

## 2. Reducción a una maximización de traza

Para $`Q`$ ortogonal,

```math
\lVert X_cQ-Y_c\rVert_F^2=\mathrm{tr}(Q^\top X_c^\top X_cQ)-2\,\mathrm{tr}(Q^\top X_c^\top Y_c)+\mathrm{tr}(Y_c^\top Y_c).
```

Los dos términos cruzados son iguales porque $`\mathrm{tr}(A)=\mathrm{tr}(A^\top)`$. **La ortogonalidad se usa en el primer término:** por la propiedad cíclica de la traza,

```math
\mathrm{tr}(Q^\top X_c^\top X_cQ)=\mathrm{tr}(X_c^\top X_c\,QQ^\top)=\mathrm{tr}(X_c^\top X_c)=\lVert X_c\rVert_F^2 .
```

Por lo tanto

```math
\lVert X_cQ-Y_c\rVert_F^2=\lVert X_c\rVert_F^2+\lVert Y_c\rVert_F^2-2\,\mathrm{tr}(Q^\top M),\qquad M=X_c^\top Y_c .
```

Como $`\lVert X_c\rVert_F^2`$ y $`\lVert Y_c\rVert_F^2`$ no dependen de $`Q`$, **minimizar la distancia equivale a maximizar $`\mathrm{tr}(Q^\top M)`$.** Es la Proposición 4.1 de [[8]](fuentes.md#ref8), aplicada a los datos centrados.

<a id="sec3"></a>

## 3. Solución en $`O(2)`$ mediante la SVD

Esta es la solución general del problema ortogonal de [[1]](fuentes.md#ref1). Sea $`M=U\Sigma V^\top`$ con $`U,V`$ ortogonales y $`\Sigma=\mathrm{diag}(\sigma_1,\sigma_2)`$, $`\sigma_1\ge\sigma_2\ge0`$. Definimos $`Z=U^\top QV`$, que es ortogonal por ser producto de ortogonales. Por la propiedad cíclica,

```math
\mathrm{tr}(Q^\top M)=\mathrm{tr}(Q^\top U\Sigma V^\top)=\mathrm{tr}(V^\top Q^\top U\,\Sigma)=\mathrm{tr}(Z^\top\Sigma)=\sigma_1z_{11}+\sigma_2z_{22}.
```

Las columnas de $`Z`$ tienen norma 1, así que $`|z_{ii}|\le1`$; como $`\sigma_i\ge0`$,

```math
\mathrm{tr}(Q^\top M)\le\sigma_1+\sigma_2 .
```

La cota se alcanza con $`Z=I`$, es decir $`U^\top QV=I`$, de donde

**Solución:**

```math
Q^*=UV^\top
```

y el costo mínimo es

```math
\lVert X_c\rVert_F^2+\lVert Y_c\rVert_F^2-2(\sigma_1+\sigma_2).
```

Es el Teorema 5.1 de [[8]](fuentes.md#ref8) con $`p=2`$ (la matriz auxiliar $`Z`$ y la cota de la traza siguen su demostración); el resultado original, válido aunque los datos no tengan rango completo, es de Schönemann [[1]](fuentes.md#ref1). Además, si $`M`$ es no singular, $`Q^*=UV^\top=M(M^\top M)^{-1/2}`$ es el factor ortogonal de la descomposición polar de $`M`$ [[8, secc. 6]](fuentes.md#ref8).

<a id="sec4"></a>

## 4. Unicidad y colinealidad

La caracterización de la unicidad es el Corolario 5.2 de [[8]](fuentes.md#ref8) con $`p=2`$: allí las soluciones son $`Q=U_rV_r^\top+U_0WV_0^\top`$ con $`W\in O(p-r)`$; con $`r=1`$, $`W\in O(1)=\lbrace\pm1\rbrace`$, y por eso aparecen exactamente dos soluciones. La pérdida respecto del máximo es

```math
(\sigma_1+\sigma_2)-\mathrm{tr}(Q^\top M)=\sigma_1(1-z_{11})+\sigma_2(1-z_{22}),
```

con ambos sumandos $`\ge0`$. Para alcanzar el máximo se necesita $`z_{ii}=1`$ en toda dirección con $`\sigma_i>0`$. En una matriz ortogonal, $`z_{ii}=1`$ obliga a que la fila $`i`$ y la columna $`i`$ sean el vector canónico $`e_i`$ (la fila tiene norma 1, así que sus demás entradas son nulas; lo mismo vale para la columna). Entonces:

- **Si $`\sigma_1,\sigma_2>0`$:** $`Z=I`$ y $`Q^*=UV^\top`$ es la **única** solución.
- **Si $`\sigma_2=0<\sigma_1`$:** solo $`z_{11}=1`$ es obligatorio, y $`Z=\mathrm{diag}(1,\pm1)`$. Hay **exactamente dos** soluciones, $`Q=UV^\top`$ y $`Q=U\,\mathrm{diag}(1,-1)V^\top`$, con el mismo costo y **determinantes opuestos** (una rotación y una reflexión). No son infinitas porque el núcleo de $`M`$ tiene dimensión 1 y $`O(1)=\lbrace\pm1\rbrace`$. (Si $`M=0`$, todo $`Q`$ ortogonal es solución.)

Por tanto, **en $`O(2)`$ la solución es única si y solo si $`M`$ es no singular.** Si $`X_c`$ (o $`Y_c`$) es colineal, $`\mathrm{rango}(M)\le1`$ y $`\sigma_2=0`$.

**Importante (casi colinealidad).** La no unicidad ocurre solo con $`\sigma_2=0`$ exacto. Si $`\sigma_2`$ es pequeño pero positivo, la solución sigue siendo única, pero las dos candidatas ($`UV^\top`$ y $`U\,\mathrm{diag}(1,-1)V^\top`$) difieren en el costo solo por $`4\sigma_2`$ (ver [sección 5](#sec5)), así que una perturbación pequeña de los datos basta para cambiar el **signo de $`\det(UV^\top)`$**: el algoritmo en $`O(2)`$ siempre devuelve $`UV^\top`$, pero esa matriz pasa de ser una rotación a ser una reflexión. La [sección 9](#sec9) cuantifica cuándo ocurre. Es lo que muestran los experimentos con nubes casi colineales (ver [Experimento C](Interpretacion_resultados.md#s5)).

<a id="sec5"></a>

## 5. Restricción a rotaciones propias ($`SO(2)`$)

Es el Teorema 7.1 de [[8]](fuentes.md#ref8), que lo demuestra para $`p`$ general con un argumento de valores propios; aquí, para $`p=2`$, se usa la forma explícita de las reflexiones. La matriz de corrección $`D`$ es la de Umeyama [[3]](fuentes.md#ref3), con justificación puramente algebraica en [[7]](fuentes.md#ref7); el problema de la "mejor rotación" proviene de Kabsch [[2]](fuentes.md#ref2), y por eso se habla del algoritmo de Kabsch–Umeyama. Sea $`Q\in SO(2)`$ y $`Z=U^\top QV`$. Entonces

```math
\det Z=\det U\cdot\det Q\cdot\det V=\det(U)\det(V)=\det(UV^\top).
```

- **Si $`\det(UV^\top)=+1`$:** $`Z=I`$ es alcanzable y la solución es $`Q^*=UV^\top`$ (ya es una rotación).
- **Si $`\det(UV^\top)=-1`$:** $`Z`$ debe tener determinante $`-1`$ y, por lo tanto, no puede ser $`I`$. Una matriz ortogonal $`2\times2`$ con determinante $`-1`$ es una reflexión, de la forma

```math
Z=\begin{pmatrix}c&s\\ s&-c\end{pmatrix},\qquad c^2+s^2=1,
```

  de modo que $`z_{22}=-z_{11}=-c`$. Entonces

```math
\mathrm{tr}(Z^\top\Sigma)=\sigma_1c-\sigma_2c=(\sigma_1-\sigma_2)\,c\le\sigma_1-\sigma_2,
```

  con igualdad para $`c=1`$, es decir $`Z=\mathrm{diag}(1,-1)`$.

En ambos casos la solución se escribe

```math
Q_{SO}=U\,D\,V^\top,\qquad D=\mathrm{diag}\big(1,\det(UV^\top)\big)
```

(y $`\det Q_{SO}=\det U\,\det D\,\det V=\det(UV^\top)^2=+1`$).

**Costo de prohibir la reflexión.** Cuando $`\det(UV^\top)=-1`$, la traza máxima baja de $`\sigma_1+\sigma_2`$ a $`\sigma_1-\sigma_2`$ (pierde $`2\sigma_2`$). Como en la [sección 2](#sec2) la traza aparece multiplicada por $`-2`$,

```math
\text{costo}_{SO(2)}-\text{costo}_{O(2)}=2\big[(\sigma_1+\sigma_2)-(\sigma_1-\sigma_2)\big]=4\sigma_2 .
```

Si $`\det(UV^\top)=+1`$ la diferencia es $`0`$ (caso de control del experimento; ver [Experimento B](Interpretacion_resultados.md#s4)).

**Unicidad en $`SO(2)`$.** Si $`\sigma_2=0<\sigma_1`$, las dos soluciones de $`O(2)`$ tienen determinantes opuestos, así que exactamente una es una rotación: **el problema restringido a $`SO(2)`$ sí tiene solución única** aunque $`M`$ sea singular. (Si $`\det(UV^\top)=-1`$ y $`\sigma_1=\sigma_2`$, en cambio, todas las rotaciones empatan, pues $`c`$ queda libre.) La condición de unicidad de $`O(2)`$ no se traslada a $`SO(2)`$ sin esta justificación.

<a id="sec6"></a>

## 6. Escala (no estimada en este proyecto)

El centrado y la escala óptima son la Proposición 8.1 de [[8]](fuentes.md#ref8) y se tratan también en [[3]](fuentes.md#ref3). Con escala $`s\ge0`$, el mismo desarrollo da $`b^\top=\bar y^\top-s\,\bar x^\top Q`$ y, para $`Q`$ fija,

```math
\lVert sX_cQ-Y_c\rVert_F^2=s^2\lVert X_c\rVert_F^2-2s\,\mathrm{tr}(Q^\top X_c^\top Y_c)+\lVert Y_c\rVert_F^2,
```

una parábola en $`s`$ cuyo mínimo está en

```math
s^*=\frac{\mathrm{tr}(Q^\top X_c^\top Y_c)}{\lVert X_c\rVert_F^2}.
```

En los experimentos se fija $`s=1`$ (solo rotación/reflexión y traslación).

<a id="sec7"></a>

## 7. Algoritmo y su correspondencia con el código

La implementación usa NumPy [[5]](fuentes.md#ref5). La función `procrustes()` está en [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py); la simulación auxiliar [`src/simulacion_criterio_relativo.py`](../src/simulacion_criterio_relativo.py) la reutiliza.

| Paso | Fórmula | En `procrustes()` |
|---|---|---|
| 1. Centrar | $`X_c=X-\mathbf{1}\bar x^\top`$, $`Y_c=Y-\mathbf{1}\bar y^\top`$ | `Xc, Yc = X - mx, Y - my` |
| 2. Matriz de relación | $`M=X_c^\top Y_c`$ | `M = Xc.T @ Yc` |
| 3. SVD | $`M=U\Sigma V^\top`$ | `U, S, Vt = np.linalg.svd(M)` |
| 4. Ortogonal / rotación | $`Q=UV^\top`$ ó $`Q=UDV^\top`$ | `Q = U @ D @ Vt` (con `D = I` en $`O(2)`$) |
| 5. Traslación | $`b^\top=\bar y^\top-\bar x^\top Q`$ | `b = my - mx @ Q` |

**Verificaciones teóricas que se comprueban numéricamente** (función `verificaciones()` de [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py); resultados en [`results/resultados.txt`](../results/resultados.txt), [`exp3_reflexion_resumen.csv`](../results/exp3_reflexion_resumen.csv) y la [interpretación](Interpretacion_resultados.md#s2)): costo mínimo $`=\lVert X_c\rVert_F^2+\lVert Y_c\rVert_F^2-2(\sigma_1+\sigma_2)`$ ([sección 3](#sec3)); diferencia de costos $`=4\sigma_2(M)`$ ([sección 5](#sec5)); coincidencia con `scipy.linalg.orthogonal_procrustes` [[6]](fuentes.md#ref6) aplicada a datos centrados (SciPy no centra ni estima traslación); ortogonalidad $`\lVert Q^\top Q-I\rVert_F=4.2\times10^{-16}`$; y, con datos exactos, recuperación de $`Q_0`$ y no de su transpuesta ($`\lVert Q-Q_0\rVert_F=7.5\times10^{-16}`$, $`\lVert Q-Q_0^\top\rVert_F=2.167`$).

<a id="sec8"></a>

## 8. Sensibilidad: por qué SO(2) es estable y O(2) no con puntos casi colineales

Esta sección y la siguiente desarrollan la observación de [[8, secc. 11.2]](fuentes.md#ref8) ("si $`\sigma_p`$ es pequeño, la solución es única pero sensible") y abordan la pregunta P1 de su agenda [[8, secc. 16]](fuentes.md#ref8): cómo se relaciona el menor valor singular de $`X^\top Y`$ con la sensibilidad de $`Q^*`$. Son desarrollo propio del grupo.

**Costo al alejarse del óptimo con una rotación.** Supongamos $`\det(UV^\top)=+1`$ (la verdad es una rotación) y tomemos una rotación cualquiera $`Q=U\,R(\theta)\,V^\top`$, es decir $`Z=R(\theta)`$. Entonces $`z_{11}=z_{22}=\cos\theta`$ y

```math
\mathrm{tr}(Z^\top\Sigma)=(\sigma_1+\sigma_2)\cos\theta
\;\Longrightarrow\;
\text{costo}(Q)-\text{costo}_{\min}=2(\sigma_1+\sigma_2)(1-\cos\theta)\approx(\sigma_1+\sigma_2)\,\theta^2 .
```

**Lectura.** Dentro de $`SO(2)`$, lo que fija el ángulo es la **suma** $`\sigma_1+\sigma_2`$, que sigue siendo grande aunque $`\sigma_2\approx0`$: girar la recta cuesta mucho, porque los puntos se separan de los datos a lo largo de toda la recta. En cambio, en $`O(2)`$ existe una competidora, la reflexión, a una distancia de costo de solo $`4\sigma_2`$ ([sección 5](#sec5)), y ese costo sí tiende a cero.

| Configuración | Curvatura en ángulo, $`SO(2)`$ | Brecha rotación–reflexión, $`O(2)`$ |
|---|---|---|
| Nube del Experimento C, $`\delta=0.001`$ (sin ruido) | $`\sigma_1+\sigma_2\approx324`$ | $`4\sigma_2\approx1.7\times10^{-4}`$ |

Por eso, **la casi colinealidad no degrada el ángulo; degrada la elección del signo del determinante.** Exigir $`\det Q=+1`$ elimina esa elección, pero solo es legítimo si se sabe de antemano que la verdad es una rotación; si fuera una reflexión, $`SO(2)`$ daría el error constante ≈ 2 del Experimento B.

<a id="sec9"></a>

## 9. Criterio relativo r: cuándo el ruido decide el signo del determinante

La [sección 4](#sec4) dice que el peligro aparece cuando $`\sigma_2(M)`$ es "pequeño". La pregunta es: ¿pequeño comparado con qué? Con el ruido.

**Señal.** Sin ruido, $`Y_c=X_cQ_0`$ y $`M=X_c^\top X_cQ_0`$. Como $`Q_0`$ es ortogonal, $`\sigma_2(M_{\text{limpia}})=\sigma_2(X_c)^2`$.

**Perturbación.** Con ruido, $`M=X_c^\top X_cQ_0+X_c^\top E_c`$. En la dirección débil de $`X_c`$ (vectores singulares $`a_2`$, $`b_2`$), la componente del ruido es $`b_2^\top X_c^\top E_c=\sigma_2(X_c)\,a_2^\top E_c`$, de tamaño típico $`\sigma_2(X_c)\,s`$, donde $`s=\varepsilon\lVert X_c\rVert_F/\sqrt{2N}`$ es la desviación del ruido por coordenada (supuesto S3).

**Cociente.**

```math
r=\frac{\sigma_2(M_{\text{limpia}})}{s\,\sigma_2(X_c)}=\frac{\sigma_2(X_c)}{s}
```

es el espesor de la nube en la dirección débil medido en unidades de ruido. Si $`r\gg1`$, el ruido no alcanza a cambiar el signo de $`\det(UV^\top)`$; si $`r\lesssim1`$, el signo queda a merced del ruido y, cuando $`r\to0`$, $`O(2)`$ elige rotación o reflexión casi al azar (≈ 50 %).

Esta deducción es de orden de magnitud. La confirmación es numérica: la simulación [`src/simulacion_criterio_relativo.py`](../src/simulacion_criterio_relativo.py) (4 niveles de ruido × 9 espesores × 400 repeticiones) muestra que la frecuencia de reflexiones erróneas depende solo de $`r`$: ≈ 0 % para $`r\ge3`$ y entre 40 y 50 % para $`r<0.3`$ (ver la [interpretación, §5.4](Interpretacion_resultados.md#s5-4)). Así, la observación de [[8, secc. 11.2]](fuentes.md#ref8) de que el espectro singular sirve como diagnóstico se precisa: hay que compararlo con el ruido.

## Referencias

La lista completa, con la numeración única del repositorio, está en [`fuentes.md`](fuentes.md): [[1]](fuentes.md#ref1), [[2]](fuentes.md#ref2), [[3]](fuentes.md#ref3), [[4]](fuentes.md#ref4), [[5]](fuentes.md#ref5), [[6]](fuentes.md#ref6), [[7]](fuentes.md#ref7), [[8]](fuentes.md#ref8).
