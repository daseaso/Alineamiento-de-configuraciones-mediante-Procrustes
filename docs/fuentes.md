# Referencias

La numeración es única para todo el repositorio: los números entre corchetes que aparecen en [`desarrollo_matematico.md`](desarrollo_matematico.md) y en [`Interpretacion_resultados.md`](Interpretacion_resultados.md) remiten a esta lista. Debajo de cada fuente se indica qué parte concreta se usó y dónde.

<a id="ref1"></a>[1] P. H. Schönemann, "A generalized solution of the orthogonal Procrustes problem," *Psychometrika*, vol. 31, no. 1, pp. 1–10, Mar. 1966, doi: 10.1007/BF02289451.

- **Uso:** planteamiento del problema ortogonal y forma de su solución general, $`Q^*=UV^\top`$, válida aunque los datos no tengan rango completo ([desarrollo, sección 3](desarrollo_matematico.md#sec3)).

<a id="ref2"></a>[2] W. Kabsch, "A solution for the best rotation to relate two sets of vectors," *Acta Crystallogr. A, Cryst. Phys. Diffr. Theor. Gen. Crystallogr.*, vol. 32, no. 5, pp. 922–923, Sep. 1976, doi: 10.1107/S0567739476001873.

- **Uso:** origen del problema de la "mejor rotación" entre dos conjuntos de vectores y nombre del algoritmo de Kabsch–Umeyama ([desarrollo, sección 5](desarrollo_matematico.md#sec5)). No usamos sus pesos ni su método con vectores propios.

<a id="ref3"></a>[3] S. Umeyama, "Least-squares estimation of transformation parameters between two point patterns," *IEEE Trans. Pattern Anal. Mach. Intell.*, vol. 13, no. 4, pp. 376–380, Apr. 1991, doi: 10.1109/34.88573.

- **Uso:** matriz de corrección $`D=\mathrm{diag}(1,\det(UV^\top))`$ para obtener rotaciones propias ([desarrollo, sección 5](desarrollo_matematico.md#sec5)), traslación por centrado $`b=\bar y-Q^\top\bar x`$ ([sección 1](desarrollo_matematico.md#sec1)) y fórmula de la escala óptima, documentada pero no estimada ([sección 6](desarrollo_matematico.md#sec6)). Su advertencia de que, sin esa corrección, se puede obtener una reflexión con datos corrompidos es el fenómeno del Experimento C.

<a id="ref4"></a>[4] M. Even, L. Ganassali, J. Maier, and L. Massoulié, "Aligning embeddings and geometric random graphs: Informational results and computational approaches for the Procrustes–Wasserstein problem," 2024, *arXiv:2405.14532*. [Online]. Available: https://arxiv.org/abs/2405.14532

- **Uso:** delimitación del alcance: el caso de correspondencias desconocidas (Procrustes–Wasserstein) queda excluido por el supuesto S1 ([desarrollo, sección 0](desarrollo_matematico.md#sec0); [interpretación, §7](Interpretacion_resultados.md#s7)).

<a id="ref5"></a>[5] C. R. Harris *et al.*, "Array programming with NumPy," *Nature*, vol. 585, no. 7825, pp. 357–362, Sep. 2020, doi: 10.1038/s41586-020-2649-2.

- **Uso:** cálculo de la SVD (`numpy.linalg.svd`), determinantes, normas y generadores aleatorios con semilla en [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py) y [`src/simulacion_criterio_relativo.py`](../src/simulacion_criterio_relativo.py).

<a id="ref6"></a>[6] SciPy Developers, "scipy.linalg.orthogonal_procrustes," *SciPy Documentation*. Accessed: Nov. 1, 2026. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html

- **Uso:** implementación independiente para verificar la nuestra en `verificaciones()` ([interpretación, §2](Interpretacion_resultados.md#s2)). De su documentación tomamos que no centra los datos (por eso le pasamos $`X_c`$ e $`Y_c`$) y que su `scale` es la suma de los valores singulares, no la escala geométrica.

<a id="ref7"></a>[7] J. Lawrence, J. Bernal, and C. Witzgall, "A purely algebraic justification of the Kabsch–Umeyama algorithm," *J. Res. Natl. Inst. Stand. Technol.*, vol. 124, Art. no. 124028, 2019, doi: 10.6028/jres.124.028.

- **Uso:** justificación puramente algebraica (sin multiplicadores de Lagrange) de que $`UDV^\top`$ es la rotación óptima, y de que el problema con traslación se reduce a centrar los datos ([desarrollo, secciones 1 y 5](desarrollo_matematico.md#sec5)).

<a id="ref8"></a>[8] O. Velásquez Castañón, "Alineamiento de configuraciones mediante Procrustes: Matrices ortogonales, SVD y estabilidad," material de referencia, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, Lima, Perú, Oct. 2026.

- **Uso:** artículo de referencia del tema y fuente principal del desarrollo:
  - Prop. 4.1, reducción a traza ([desarrollo, sección 2](desarrollo_matematico.md#sec2));
  - Teo. 5.1, solución $`Q^*=UV^\top`$ con la matriz auxiliar $`Z`$, y secc. 6, descomposición polar ([sección 3](desarrollo_matematico.md#sec3));
  - Cor. 5.2, unicidad ([sección 4](desarrollo_matematico.md#sec4));
  - Teo. 7.1, rotaciones propias y costo $`4\sigma_2`$ ([sección 5](desarrollo_matematico.md#sec5));
  - Prop. 8.1, traslación y escala ([secciones 1 y 6](desarrollo_matematico.md#sec1));
  - secc. 11.2 y 11.3, sensibilidad y justificación de la norma de Frobenius ([secciones 0 y 8](desarrollo_matematico.md#sec8));
  - secc. 13 y 14, algoritmo, verificaciones y experimento obligatorio (código);
  - secc. 16, pregunta P1, que motiva las [secciones 8 y 9](desarrollo_matematico.md#sec9).

## Aporte propio del grupo

No provienen de ninguna fuente: la comparación entre la curvatura $`(\sigma_1+\sigma_2)\theta^2`$ de SO(2) y la brecha $`4\sigma_2`$ de O(2) ([desarrollo, sección 8](desarrollo_matematico.md#sec8)); el criterio $`r=\sigma_2(X_c)/s`$ con su simulación de 400 repeticiones ([sección 9](desarrollo_matematico.md#sec9)); y la predicción $`\sigma_\theta\approx\varepsilon/\sqrt{2N}`$ del error angular ([interpretación, §3.3](Interpretacion_resultados.md#s3)).
