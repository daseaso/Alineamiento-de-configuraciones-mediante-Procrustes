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
    └── (documentos PDF y uso_ia.md)
```

- `src/`: contiene el script principal de Python con el algoritmo y los experimentos.
- `results/`: archivos generados por el código (tablas `.csv`, gráficos `.png` y el registro `resultados.txt`). Se incluyen ya generados para poder revisarlos sin ejecutar nada; ejecutar el script los regenera con el mismo contenido (semilla fija).
- `requirements.txt`: dependencias y versiones exactas necesarias para la ejecución.
- `docs/`: Documento PDF con el desarrollo matemático completo, interpretación de resultados y guion del video.

## Instrucciones de ejecución

**Requisitos:** Python 3.14.5 (versión con la que se obtuvieron los resultados) y las bibliotecas de `requirements.txt` (NumPy, Matplotlib y SciPy, con sus versiones).

1. Descargar o clonar el repositorio y abrir una terminal en la carpeta raíz.
2. Instalar las dependencias:

```
pip install -r requirements.txt
```

3. Entrar a la carpeta de resultados y ejecutar el script desde ahí (los archivos se generan en la carpeta desde la que se ejecuta):

```
cd results
python ../src/experimentos_procrustes.py
```

(En Windows PowerShell, escribir cada comando en su propia línea.)

El programa tarda unos segundos, imprime los resultados en pantalla y regenera en `results/`: `resultados.txt`, los archivos `.csv` y las figuras `fig1` a `fig4`. Como la semilla es fija (`SEED = 20261101`), el contenido debe ser idéntico al incluido en el repositorio, salvo diferencias de ~1e-16 en valores que valen cero.

**Cómo comprobar que la ejecución es correcta** (valores esperados en `resultados.txt`):

| Comprobación | Valor esperado |
|---|---|
| Rotación con ruido ε = 0.01: error angular | 0.066° ± 0.038° |
| Rotación con ruido ε = 0.50: error angular | 2.414° ± 2.270° |
| Reflexión: det Q en O(2) | −1 en las 140 corridas |
| Reflexión: máx. \|(costo_SO2 − costo_O2) − 4σ₂\| | 2.842e-13 |
| Casi colineal, δ = 0.001: reflexiones erróneas en O(2) | 70 % |
| Comparación con SciPy: ‖Q propio − Q SciPy‖ | 0.00e+00 |

## Resultados principales

Los experimentos usan 30 puntos sintéticos en el plano, una rotación verdadera de 50° y 20 repeticiones por nivel. Se separa el **residuo de ajuste** (contra los datos con ruido) del **error de recuperación** (contra la transformación verdadera).

**1. Rotación con ruido (SO(2)).** Sin ruido, el error es de precisión de máquina (~1e-15). Con ruido, el error angular crece casi linealmente con el nivel de ruido ε:

| ε | Error angular (°) | RMSE de ajuste | RMSE vs señal limpia |
|---|---|---|---|
| 0.01 | 0.066 ± 0.038 | 0.0472 | 0.0094 |
| 0.10 | 0.592 ± 0.451 | 0.4692 | 0.1019 |
| 0.20 | 1.079 ± 0.801 | 0.9449 | 0.1825 |
| 0.50 | 2.414 ± 2.270 | 2.341 | 0.519 |

**2. Reflexión: O(2) frente a SO(2).** Cuando la transformación verdadera es una reflexión, O(2) la detecta siempre (det = −1 en las 140 corridas). Forzar SO(2) aumenta el costo cuadrático en exactamente 4σ₂(M); el error máximo de esa igualdad fue 2.8e-13. Con verdad = rotación, la diferencia de costos es 0.

**3. Puntos casi colineales.** Al reducir el espesor δ de la nube, σ₂(M) baja de 43 a ≈1e-4 y el residuo de ajuste se mantiene constante (≈0.032). Sin embargo, en O(2) la fracción de reflexiones erróneas sube a 10 % (δ = 0.003) y 70 % (δ = 0.001), con error de recuperación de 0.20 y 1.40; en SO(2) el error se mantiene en ≈1.5e-3. **Un buen ajuste no garantiza identificar la transformación.**

**4. Verificaciones.** La solución propia coincide con `scipy.linalg.orthogonal_procrustes` (diferencia 0.00e+00), cumple ‖QᵀQ − I‖ ≈ 4e-16 y su costo coincide con la fórmula ‖Xc‖² + ‖Yc‖² − 2Σσᵢ (diferencia ≈ 2e-13).

Figuras (en `results/`):

![Puntos antes y después](results/fig1_puntos_antes_despues.png)
![Errores vs ruido](results/fig2_errores_vs_ruido.png)
![Casi colinealidad](results/fig3_casi_colineal.png)
![Puntos casi colineales](results/fig4_colineal_puntos.png)

El desarrollo matemático completo, las tablas y la interpretación están en [`docs/avance_parcial.pdf`](docs/avance_parcial.pdf).


## Enlace al video

Video final: _(pendiente; se agregará en la entrega final)_