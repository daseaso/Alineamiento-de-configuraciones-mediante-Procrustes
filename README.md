# Recuperación de transformaciones rígidas bidimensionales bajo ruido mediante Procrustes ortogonal

**Tema asignado:** Alineamiento de configuraciones mediante Procrustes
**Estado del trabajo:** Entrega Parcial (etiqueta: `avance-parcial`)

## Integrantes
* Clara Jimena Palma Dolmos
* Yanira Thais Sanchez Cárdenas
* Daniel Alejandro Seas Ordinola
* Alyssa Antuanette Trujillo Cruzado

## Problema
El problema consiste en recuperar analíticamente una transformación rígida (rotación y traslación) que alinee dos conjuntos de puntos correspondientes, y determinar la confiabilidad de este ajuste matemático cuando los datos originales presentan ruido o direcciones poco observables.

## Objetivo
Formular y demostrar matemáticamente la solución del problema de Procrustes mediante la descomposición en valores singulares (SVD) para minimizar la discrepancia cuadrática medida en norma de Frobenius. Adicionalmente, implementar esta solución en Python para contrastar el error de ajuste geométrico frente al error real de recuperación de parámetros.

## Organización de archivos
Los archivos correspondientes a esta entrega se estructuran de la siguiente manera:
* `src/`: Contiene el script principal de Python con el algoritmo y los experimentos.
* `results/`: Archivos generados por el código (tablas `.csv`, gráficos `.png` y el registro `resultados.txt`).
* `requirements.txt`: Dependencias y versiones exactas necesarias para la ejecución.
* `docs/`: Documento PDF con el desarrollo matemático completo, interpretación de resultados y guion del video.