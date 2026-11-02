# Fuentes citadas en el desarrollo

1. **Schönemann, P. H. (1966).** A generalized solution of the orthogonal Procrustes problem. *Psychometrika, 31*(1), 1–10.
   - Uso en el proyecto: base de la solución por SVD, Q\* = U·Vᵀ (inciso 2.b), y de la reducción a traza (inciso 2.a).

2. **Kabsch, W. (1976).** A solution for the best rotation to relate two sets of vectors. *Acta Crystallographica Section A, 32*, 922–923.
   - Uso en el proyecto: restricción a rotaciones propias con det Q = 1, mediante la matriz D (inciso 2.d).

3. **Umeyama, S. (1991).** Least-squares estimation of transformation parameters between two point patterns. *IEEE Transactions on Pattern Analysis and Machine Intelligence, 13*(4), 376–380.
   - Uso en el proyecto: tratamiento de la traslación mediante centrado y de la escala óptima (inciso 2.e).

4. **Even, M., Ganassali, L., Maier, J. y Massoulié, L. (2024).** Aligning embeddings and geometric random graphs: informational results and computational approaches for the Procrustes–Wasserstein problem. *arXiv:2405.14532*. https://arxiv.org/abs/2405.14532
   - Uso en el proyecto: delimitación del alcance; el caso de correspondencias desconocidas queda excluido de esta etapa.

5. **Harris, C. R. et al. (2020).** Array programming with NumPy. *Nature, 585*, 357–362.
   - Uso en el proyecto: cálculo de la SVD (`numpy.linalg.svd`), determinantes y normas en la implementación.
