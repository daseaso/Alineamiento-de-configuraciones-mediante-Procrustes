# Declaración de uso de inteligencia artificial

**Proyecto:** Recuperación de transformaciones rígidas bidimensionales bajo ruido mediante Procrustes ortogonal

**Curso:** Álgebra Lineal Aplicada 2026-2

El grupo usó herramientas de IA generativa como apoyo en las partes indicadas abajo. Todo el contenido fue revisado por los integrantes, que asumen la responsabilidad por él. El planteamiento del problema, las decisiones de diseño de los experimentos y las conclusiones son del grupo.

## 1. Resumen

| # | Herramienta y versión | Finalidad | Parte del trabajo |
|---|---|---|---|
| 1 | Claude (Anthropic), modelo Claude Sonnet 5.5 | Apoyo en la redacción de la interpretación de resultados a partir del avance y del código del grupo | [`docs/Interpretacion_resultados.md`](Interpretacion_resultados.md) |
| 2 | Claude (Anthropic), modelo Claude Opus 5.5, 16/11/2026 | Revisión de coherencia del repositorio antes de la entrega y apoyo en las correcciones | `README.md`, [`docs/fuentes.md`](fuentes.md), [`docs/guion_video.md`](guion_video.md), [`docs/contribuciones.md`](contribuciones.md), [`docs/Interpretacion_resultados.md`](Interpretacion_resultados.md), función `verificaciones()` de [`src/procrustes_experimentos.py`](../src/procrustes_experimentos.py) y carpeta `results/` |

## 2. Instrucciones (prompts) relevantes

**Uso 1 — interpretación de resultados.** 

**Uso 2 — revisión y corrección del repositorio.** Prompts principales, en síntesis:

1. Se compartieron los comentarios de revisión del grupo (nombres de archivo en las instrucciones, verificaciones faltantes en los resultados, valores del README, fuentes no citadas y contradicciones entre el guion y la interpretación) y el repositorio completo, con la instrucción.
2. *"¿Hay contradicciones entre el guion y las tablas que se mencionan?"*

## 3. Qué hizo la IA en cada parte

- **Interpretación de resultados (uso 1):** propuesta de redacción de las secciones de interpretación a partir de las cifras del script; el grupo seleccionó y editó el contenido.
- **Revisión (uso 2):**
  - Detectó el nombre incorrecto del script en las instrucciones, enlaces a archivos inexistentes (`docs/avance_parcial.pdf`, `docs/uso_ia.md`) y que `results/resultados.txt` venía de una versión anterior del código, sin la sección de verificaciones.
  - Añadió a `verificaciones()` la impresión de det Q, del error con datos exactos y del contraste ‖Q − Q₀ᵀ‖, que ya figuraban en la interpretación pero el código no imprimía.
  - Unificó las tres listas de referencias en una sola lista IEEE (`docs/fuentes.md`) y retiró las fuentes que no se citaban en ningún documento.
  - Alineó sus afirmaciones con la interpretación (estabilidad de SO(2) en el Exp. C, diagnóstico con σ₂(M) frente al ruido, comportamiento con ruido alto).


## 4. Cómo se verificó

- Se reejecutó `src/procrustes_experimentos.py` desde `results/` y se comprobó que todas las cifras del README, de la interpretación y del guion coinciden con el nuevo `resultados.txt` y con los CSV (semilla fija `SEED = 20261101`).
- Los máximos de |costo_SO2 − costo_O2 − 4σ₂| por nivel de ruido se recalcularon directamente desde `exp3_reflexion_corridas.csv`.
- La afirmación de la sección 5.2 de la interpretación (≈ 21 % y ≈ 38 % de reflexiones con 4000 repeticiones) se recalculó por separado y dio ≈ 20 % y ≈ 39 %.
- Los datos bibliográficos de `docs/fuentes.md` se contrastaron con los DOI y enlaces de cada fuente.
- Los integrantes leyeron los documentos corregidos antes de subirlos. 
