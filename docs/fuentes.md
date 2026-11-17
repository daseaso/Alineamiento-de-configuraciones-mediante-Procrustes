# Fuentes citadas en el proyecto (formato IEEE)

La numeración es única para todo el repositorio: los números entre corchetes que aparecen en desarrollo_matematico.md y en `Interpretacion_resultados.md` remiten a esta lista. Debajo de cada fuente se indica en qué parte del trabajo se usa.

[1] P. H. Schönemann, "A generalized solution of the orthogonal Procrustes problem," *Psychometrika*, vol. 31, no. 1, pp. 1–10, Mar. 1966, doi: 10.1007/BF02289451.

- **Uso:** formulación general del problema ortogonal y su solución por SVD, $`Q^*=UV^\top`$, y condición de unicidad (desarrollo matemático, secciones 3 y 4). En la interpretación: algoritmo (§1.1) y unicidad teórica frente a ambigüedad práctica (§5.3, punto 6).

[2] W. Kabsch, "A solution for the best rotation to relate two sets of vectors," *Acta Crystallogr. A, Cryst. Phys. Diffr. Theor. Gen. Crystallogr.*, vol. 32, no. 5, pp. 922–923, Sep. 1976, doi: 10.1107/S0567739476001873.

- **Uso:** corrección del signo del determinante (matriz $`D`$) para obtener rotaciones propias (desarrollo matemático, sección 5). En la interpretación: algoritmo en SO(2) (§1.1) y verificación del sobrecosto $`4\sigma_2(M)`$ (§4.3, punto 2).

[3] S. Umeyama, "Least-squares estimation of transformation parameters between two point patterns," *IEEE Trans. Pattern Anal. Mach. Intell.*, vol. 13, no. 4, pp. 376–380, Apr. 1991, doi: 10.1109/34.88573.

- **Uso:** igual que [2] para la matriz $`D`$; además, centrado (eliminación de la traslación) y escala óptima (desarrollo matemático, secciones 1, 5 y 6). En la interpretación: §1.1 y §4.3.

[4] M. Even, L. Ganassali, J. Maier, and L. Massoulié, "Aligning embeddings and geometric random graphs: Informational results and computational approaches for the Procrustes–Wasserstein problem," 2024, *arXiv:2405.14532*. [Online]. Available: https://arxiv.org/abs/2405.14532

- **Uso:** delimitación del alcance; el caso de correspondencias desconocidas (Procrustes–Wasserstein) queda fuera de esta etapa (README, "Fuera del alcance"; interpretación, §7 Limitaciones).

[5] C. R. Harris *et al.*, "Array programming with NumPy," *Nature*, vol. 585, no. 7825, pp. 357–362, Sep. 2020, doi: 10.1038/s41586-020-2649-2.

- **Uso:** implementación: SVD (`numpy.linalg.svd`), determinantes, normas y generación de números aleatorios en `src/procrustes_experimentos.py` (interpretación, §1.1).

[6] SciPy Developers, "scipy.linalg.orthogonal_procrustes," *SciPy Documentation*. Accessed: Nov. 1, 2026. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html

- **Uso:** implementación independiente con la que se compara la solución propia (función `verificaciones()`; interpretación, §2).

[7] J. Lawrence, J. Bernal, and C. Witzgall, "A purely algebraic justification of the Kabsch–Umeyama algorithm," *J. Res. Natl. Inst. Stand. Technol.*, vol. 124, Art. no. 124028, 2019, doi: 10.6028/jres.124.028.

- **Uso:** justificación puramente algebraica (SVD y traza) de la corrección con la matriz $`D`$ para rotaciones propias (desarrollo matemático, sección 5).

[8] O. Velásquez Castañón, "Alineamiento de configuraciones mediante Procrustes: Matrices ortogonales, SVD y estabilidad," material de referencia, Álgebra Lineal Aplicada 2026-2, Universidad del Pacífico, Lima, Perú, Oct. 2026.

- **Uso:** planteamiento del tema; matriz auxiliar, caracterización de la unicidad y proposiciones sobre traslación y escala (desarrollo matemático, secciones 1, 3, 4 y 6).
