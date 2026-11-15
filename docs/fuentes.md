# Fuentes citadas en el desarrollo
 
1. **Even, M., Ganassali, L., Maier, J. y Massoulié, L. (2024).** Aligning embeddings and geometric random graphs: informational results and computational    approaches for the Procrustes–Wasserstein problem. arXiv:2405.14532. https://arxiv.org/abs/2405.14532

    -Uso en el proyecto: delimitación del alcance; el caso de correspondencias desconocidas queda excluido de esta etapa.

3. **Harris, C. R. et al. (2020).** Array programming with NumPy. Nature, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2

   -Uso en el proyecto: cálculo de la SVD (numpy.linalg.svd), determinantes y normas en la implementación.
   
5. **Gower, J. C. (1975).** Generalized Procrustes analysis. Psychometrika, 40(1), 33–51. https://doi.org/10.1007/BF02291478

   -Uso en el proyecto: delimitación del alcance; el análisis de Procrustes generalizado (más de dos configuraciones, con una configuración consenso) queda como extensión fuera de esta etapa.
   
7. **ten Berge, J. M. F. (1977).** Orthogonal Procrustes rotation for two or more matrices. Psychometrika, 42(2), 267–276. https://three-mode.leidenuniv.nl/pdf/t/tenberge1977pmet.pdf

   -Uso en el proyecto: condición de optimalidad de la solución (inciso 2.b), caso singular (inciso 2.c) y delimitación del alcance.
   
9. **Lawrence, J., Bernal, J. y Witzgall, C. (2019).** A purely algebraic justification of the Kabsch–Umeyama algorithm. Journal of Research of the National Institute of Standards and Technology, 124, art. 124028. https://doi.org/10.6028/jres.124.028

   -Uso en el proyecto: justificación puramente algebraica (SVD y traza) de la corrección con la matriz D para rotaciones propias (inciso 2.d).
   
11. **Arias-Castro, E., Javanmard, A. y Pelletier, B. (2019).** Perturbation bounds for Procrustes, classical scaling, and trilateration, with applications to manifold learning. arXiv:1810.09569. https://arxiv.org/abs/1810.09569
  
  -Uso en el proyecto: estabilidad de la solución frente al ruido. Respalda los Experimentos A y C y la interpretación de resultados.
