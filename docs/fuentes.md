# Referencias

La numeración es única para todo el repositorio: los números entre corchetes que aparecen en [`desarrollo_matematico.md`](desarrollo_matematico.md) y en [`Interpretacion_resultados.md`](Interpretacion_resultados.md) remiten a esta lista. Debajo de cada fuente se indica en qué parte del trabajo se usa.

<a id="ref4"></a>[1] M. Even, L. Ganassali, J. Maier, and L. Massoulié, "Aligning embeddings and geometric random graphs: Informational results and computational approaches for the Procrustes–Wasserstein problem," 2024, *arXiv:2405.14532*. [Online]. Available: https://arxiv.org/abs/2405.14532

- **Uso:** delimitación del alcance; el caso de correspondencias desconocidas (Procrustes–Wasserstein) queda fuera de esta etapa (README, "Fuera del alcance"; interpretación, [§7 Limitaciones](Interpretacion_resultados.md#s7)).

<a id="ref5"></a>[2] C. R. Harris *et al.*, "Array programming with NumPy," *Nature*, vol. 585, no. 7825, pp. 357–362, Sep. 2020, doi: 10.1038/s41586-020-2649-2.

- **Uso:** implementación: SVD (`numpy.linalg.svd`), determinantes, normas y generación de números aleatorios en [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py) (desarrollo matemático, [sección 7](desarrollo_matematico.md#sec7); interpretación, [§1.1](Interpretacion_resultados.md#s1-1)).

<a id="ref6"></a>[3] SciPy Developers, "scipy.linalg.orthogonal_procrustes," *SciPy Documentation*. Accessed: Nov. 1, 2026. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html

- **Uso:** implementación independiente con la que se compara la solución propia (función `verificaciones()` de [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py); desarrollo matemático, [sección 7](desarrollo_matematico.md#sec7); interpretación, [§2](Interpretacion_resultados.md#s2)).

<a id="ref7"></a>[4] J. Lawrence, J. Bernal, and C. Witzgall, "A purely algebraic justification of the Kabsch–Umeyama algorithm," *J. Res. Natl. Inst. Stand. Technol.*, vol. 124, Art. no. 124028, 2019, doi: 10.6028/jres.124.028.

- **Uso:** justificación puramente algebraica (SVD y traza) de la corrección con la matriz $`D`$ para rotaciones propias (desarrollo matemático, [sección 5](desarrollo_matematico.md#sec5)).

<a id="ref8"></a>[8] O. Velásquez Castañón, "Alineamiento de configuraciones mediante Procrustes: Matrices ortogonales, SVD y estabilidad," material de referencia, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, Lima, Perú, Oct. 2026.

- **Uso:** planteamiento del tema; matriz auxiliar, caracterización de la unicidad y proposiciones sobre traslación y escala (desarrollo matemático, [sección 1](desarrollo_matematico.md#sec1), [sección 3](desarrollo_matematico.md#sec3), [sección 4](desarrollo_matematico.md#sec4) y [sección 6](desarrollo_matematico.md#sec6)).
