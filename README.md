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

```
├── README.md
├── requirements.txt
├── src/
│   └── experimentos_procrustes.py
├── results/
│   ├── resultados.txt
│   ├── puntos_X.csv
│   ├── exp2_rotacion_resumen.csv      (y exp2_rotacion_corridas.csv)
│   ├── exp3_reflexion_resumen.csv     (y exp3_reflexion_corridas.csv)
│   ├── exp4_colineal_resumen.csv      (y exp4_colineal_corridas.csv)
│   └── fig1 … fig4 (.png)
└── docs/
    └── (documentos PDF y uso_ia.md; ver tabla siguiente)
```

- `src/`: contiene el script principal de Python con el algoritmo y los experimentos.
- `results/`: archivos generados por el código (tablas `.csv`, gráficos `.png` y el registro `resultados.txt`). Se incluyen ya generados para poder revisarlos sin ejecutar nada; ejecutar el script los regenera con el mismo contenido (semilla fija).
- `requirements.txt`: dependencias y versiones exactas necesarias para la ejecución.
- `docs/`: Documento PDF con el desarrollo matemático completo, interpretación de resultados y guion del video.
