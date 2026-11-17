"""Simulacion auxiliar: frecuencia de reflexiones erroneas en O(2) segun el criterio relativo
   r = sigma2(M) / ( s * sigma2(X) ),  s = eps * ||Xc||_F / sqrt(2N)

Se usa la misma nube casi colineal del Experimento C. Colocar este archivo junto a
experimentos_procrustes.py (carpeta src/) y ejecutar:  python src/simulacion_criterio_relativo.py
Salida: tabla en pantalla y results/simulacion_criterio_relativo.csv
"""
import csv
import os
import numpy as np
import experimentos_procrustes as ep

EPSILONS = [0.005, 0.01, 0.02, 0.05]
DELTAS = [0.0003, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.1]
REPS = 400
SEMILLA = 1

t, w = ep.make_collinear_base()
N = ep.N
rng = np.random.default_rng(SEMILLA)
filas = []
for eps in EPSILONS:
    for d in DELTAS:
        Xd = np.column_stack([t, d * w])
        Xdc = Xd - Xd.mean(0)
        sx = np.linalg.svd(Xdc, compute_uv=False)
        Yclean = Xd @ ep.Q0 + ep.B0
        _, _, S0 = ep.procrustes(Xd, Yclean, False)      # valores singulares de M sin ruido
        s = eps * np.linalg.norm(Xdc) / np.sqrt(2 * N)
        r = S0[1] / (s * sx[1])
        flips = 0
        for _ in range(REPS):
            Z = rng.normal(size=Xd.shape)
            E = eps * np.linalg.norm(Xdc) * Z / np.linalg.norm(Z)
            QO, _, _ = ep.procrustes(Xd, Yclean + E, False)
            flips += np.linalg.det(QO) < 0
        filas.append(dict(r=r, frac_reflexiones=flips / REPS, eps=eps, delta=d))

filas.sort(key=lambda f: f["r"])
print("      r   frac. reflexiones   eps    delta")
for f in filas:
    print(f"{f['r']:8.3f}   {f['frac_reflexiones']:8.2f}        {f['eps']:<6} {f['delta']}")

carpeta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
os.makedirs(carpeta, exist_ok=True)
with open(os.path.join(carpeta, "simulacion_criterio_relativo.csv"), "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=["r", "frac_reflexiones", "eps", "delta"])
    wr.writeheader(); wr.writerows(filas)
