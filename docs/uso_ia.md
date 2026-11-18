# Declaración de uso de inteligencia artificial

**Proyecto:** Recuperación de transformaciones rígidas bidimensionales bajo ruido mediante Procrustes ortogonal

**Curso:** Álgebra Lineal Aplicada 2026-2

El grupo usó herramientas de IA generativa como apoyo en las partes indicadas abajo. Todo el contenido fue revisado por los integrantes, que asumen la responsabilidad por él. El planteamiento del problema, las decisiones de diseño de los experimentos y las conclusiones son del grupo.

## 1. Resumen

| # | Herramienta y versión | Finalidad | Parte del trabajo |
|---|---|---|---|
| 1 | Claude (Anthropic), modelo Claude Sonnet 5.5 | Apoyo en la redacción de la interpretación de resultados a partir del avance y del código del grupo | [`docs/Interpretacion_resultados.md`](Interpretacion_resultados.md) |
| 2 | Claude (Anthropic), modelo Claude Opus 5.5, 16/11/2026 | Revisión de coherencia del repositorio antes de la entrega y apoyo en las correcciones | [`docs/desarrollo_matematico.md`](desarrollo_matematico.md), [`README.md`](../README.md), [`docs/fuentes.md`](fuentes.md), [`docs/guion_video.md`](guion_video.md), [`docs/contribuciones.md`](contribuciones.md), [`docs/Interpretacion_resultados.md`](Interpretacion_resultados.md), función `verificaciones()` de [`src/experimentos_procrustes.py`](../src/experimentos_procrustes.py) y carpeta `results/` |
| 3 | Claude (Anthropic), Claude Opus 5.5, 16/11/2026 | Apoyo en las diapositivas y en la redacción del guion | Diapositivas del video y [`docs/guion_video.md`](guion_video.md) |
| 4 | Claude (Anthropic),  Claude Opus 5.5, 16/11/2026 | Generación del script de la simulación auxiliar de 400 repeticiones | [`src/simulacion_criterio_relativo.py`](../src/simulacion_criterio_relativo.py), [`results/simulacion_criterio_relativo.csv`](../results/simulacion_criterio_relativo.csv) |

## 2. Instrucciones (prompts) relevantes

**Uso 1 — interpretación de resultados.** Se compartieron el avance, el código y `results/resultados.txt`, y se pidió redactar la interpretación de los tres experimentos separando el residuo de ajuste del error de recuperación y relacionando las cifras con la teoría.

**Uso 2 — revisión y corrección del repositorio.** Prompts principales, en síntesis:

1. Se compartieron los comentarios de revisión del grupo (nombres de archivo en las instrucciones, verificaciones faltantes en los resultados, valores del README, fuentes no citadas y contradicciones entre el guion y la interpretación) y el repositorio completo, con la instrucción de revisar qué faltaba corregir para la entrega final.
2. Revisión del guion del video: duración y correspondencia con las diapositivas.
   
**Uso 3 - presentación y guion.** Se compartieron la presentación y el guion con la instrucción de señalar problemas de legibilidad y duración, y de proponer una versión más clara de las diapositivas de los experimentos y de las intervenciones.


## 3. Qué hizo la IA en cada parte

- **Interpretación de resultados (uso 1):** propuesta de redacción de las secciones de interpretación a partir de las cifras del script; el grupo seleccionó y editó el contenido.
- **Revisión (uso 2):**
  - Detectó el nombre incorrecto del script en las instrucciones, enlaces a archivos inexistentes (`docs/avance_parcial.pdf`, `docs/uso_ia.md`) y que `results/resultados.txt` venía de una versión anterior del código, sin la sección de verificaciones.
  - Añadió a `verificaciones()` la impresión de det Q, del error con datos exactos y del contraste ‖Q − Q₀ᵀ‖, que ya figuraban en la interpretación pero el código no imprimía.
  - Unificó las tres listas de referencias en una sola lista IEEE (`docs/fuentes.md`) y retiró las fuentes que no se citaban en ningún documento.
  - Alineó sus afirmaciones con la interpretación (estabilidad de SO(2) en el Exp. C, diagnóstico con σ₂(M) frente al ruido, comportamiento con ruido alto).
  - Propuso la redacción de las secciones nuevas del desarrollo matemático: supuestos (sección 0), sentido de giro, sensibilidad (σ₁+σ₂)θ² frente a 4σ₂ (sección 8) y la justificación teórica del criterio r = σ₂(Xc)/s (sección 9), así como las citas a los resultados concretos de [8].
    - Propuso la redacción de la sección 5.4 de la interpretación y verificó que el script de la simulación reproduce el CSV.
- **Simulación auxiliar (uso 3):** la IA generó el script [`src/simulacion_criterio_relativo.py`](../src/simulacion_criterio_relativo.py) (4 niveles de ruido × 9 espesores × 400 repeticiones y el cálculo del cociente r) y, al ejecutarlo, el archivo [`results/simulacion_criterio_relativo.csv`](../results/simulacion_criterio_relativo.csv). El grupo planteó la pregunta (si importa más el espesor de la nube o el ruido), revisó el código, lo ejecutó para comprobar que reproduce el CSV e interpretó los resultados.
- **Presentación y guion (uso 4):** la IA creó un modelo/sugerencia de organización de las diapositivas. Se tuvo que adaptar con imágenes, tablas, resultados e interpretación del proyecto.


## 4. Cómo se verificó

- Se reejecutó `src/experimentos_procrustes.py` desde `results/` y se comprobó que todas las cifras del README, de la interpretación y del guion coinciden con el nuevo `resultados.txt` y con los CSV (semilla fija `SEED = 20261101`).
- Los máximos de |costo_SO2 − costo_O2 − 4σ₂| por nivel de ruido se recalcularon directamente desde `exp3_reflexion_corridas.csv`.
- Se reejecutó `src/simulacion_criterio_relativo.py` y reproduce exactamente `results/simulacion_criterio_relativo.csv`.
- Los datos bibliográficos de `docs/fuentes.md` se contrastaron con los DOI y enlaces de cada fuente.
- Los integrantes leyeron los documentos corregidos antes de subirlos. 
