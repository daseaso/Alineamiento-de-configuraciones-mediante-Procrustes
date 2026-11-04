# Desarrollo Matemático: Problema de Procrustes Ortogonal

**Objetivo:** Encontrar la matriz ortogonal $Q$ y el vector de traslación $b$ que minimizan la discrepancia cuadrática $f(Q, b) = ||XQ + \mathbf{1}b^T - Y||_F^2$.

**1. Eliminación de la traslación (Centrado)**
Para encontrar el óptimo respecto a $b$, derivamos $f(Q, b)$ e igualamos a cero. La solución indica que la traslación óptima compensa la diferencia entre los centroides de ambos conjuntos:

$$ 
b^T = \bar{y}^T - \bar{x}^T Q 
$$

Al sustituir $b$ en la función original, la dependencia de la ubicación espacial desaparece. El problema se reduce a alinear las matrices centradas $X_c$ e $Y_c$:

$$ 
\min_{Q} ||X_c Q - Y_c||_F^2 
$$

**2. Expansión y reducción a traza**
Desarrollamos la norma de Frobenius al cuadrado. Usando la propiedad de ortogonalidad $Q^T Q = I$, la norma de $X_c Q$ se preserva:

$$ 
||X_c Q - Y_c||_F^2 = ||X_c||_F^2 + ||Y_c||_F^2 - 2 \text{tr}(Q^T X_c^T Y_c) 
$$

Dado que las normas de las nubes de puntos son constantes, minimizar la distancia total es matemáticamente equivalente a **maximizar la traza** del producto cruzado. Definimos la matriz de covarianza cruzada $M = X_c^T Y_c$, quedando el problema como:

$$ 
\max_{Q} \text{tr}(Q^T M) 
$$

**3. Solución general en O(2) mediante SVD**
Aplicamos la Descomposición en Valores Singulares a la matriz de relación: $M = U \Sigma V^T$. Sustituyendo en la traza:

$$ 
\text{tr}(Q^T U \Sigma V^T) = \text{tr}((V^T Q^T U) \Sigma) 
$$

Definimos $Z = V^T Q^T U$, que también es una matriz ortogonal. Como los elementos diagonales de $Z$ no pueden superar 1, la traza máxima posible es la suma de los valores singulares $\sigma_1 + \sigma_2$. Este valor máximo se alcanza únicamente cuando $Z = I$, lo que permite despejar $Q$:

$$ 
Q^* = U V^T 
$$

**4. Condición de unicidad y colinealidad**
Si $\sigma_2 > 0$, la solución $Q^*$ es única. Sin embargo, si los puntos son colineales, $\sigma_2 \to 0$. Matemáticamente, esto significa que el sistema pierde información direccional, generando dos soluciones en el plano (una rotación y una reflexión) que minimizan el costo por igual.

**5. Restricción a rotaciones propias (SO(2))**
El grupo $O(2)$ admite reflexiones ($\det(Q) = -1$). Si se requiere estrictamente una rotación pura ($SO(2)$ con $\det(Q) = +1$) y la matriz $Q^*$ resultó ser una reflexión, ajustamos el signo del valor singular más pequeño mediante la matriz $D = \text{diag}(1, -1)$:

$$ 
Q_{SO} = U D V^T 
$$

La traza máxima se reduce de $\sigma_1 + \sigma_2$ a $\sigma_1 - \sigma_2$. Al reemplazar en la función de costo, se demuestra analíticamente que forzar una rotación cuando los datos provienen de una reflexión eleva el error cuadrático agregando una penalización exacta de $4\sigma_2$.
