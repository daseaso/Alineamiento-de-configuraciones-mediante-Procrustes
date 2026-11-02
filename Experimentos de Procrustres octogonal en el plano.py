"""
Experimentos de Procrustes ortogonal en el plano (sin escala).

    Modelo:  Y = X Q0 + 1 b0^T + E
    Solucion: Q = U D V^T  (SVD de M = Xc^T Yc),   b = ybar - Q^T xbar

Estructura (cada seccion cumple un punto de la consigna):
    SECCION 0 : Parametros y semilla                     -> Punto 1 (semilla fija y publicada)
    SECCION A : Nucleo (procrustes) y utilidades comunes
    SECCION B : Construccion de los puntos X             -> Punto 1
    SECCION C : Experimento rotacion + niveles de ruido  -> Punto 2
    SECCION D : Experimento reflexion, O(2) vs SO(2)     -> Punto 3 (verifica 4*sigma2)
    SECCION E : Casi colinealidad                        -> Punto 4
    SECCION F : Repeticiones, media y dispersion         -> Punto 5 (usado en C, D, E)
    SECCION G : Figuras antes/despues, ejes iguales      -> Punto 5
    SECCION H : Programa principal
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import orthogonal_procrustes

# =====================================================================
# SECCION 0. PARAMETROS Y SEMILLA                         (Punto 1)
# =====================================================================
SEED = 20261101                                   # semilla publicada: todo es reproducible
N = 30                                            # numero de puntos (>= 20)
THETA0 = np.deg2rad(50.0)                         # rotacion verdadera
B0 = np.array([2.5, -1.5])                        # traslacion verdadera

NOISE = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]    # niveles ||E||_F / ||Xc||_F   (Punto 2)
REPS = 20                                         # repeticiones por nivel (>= 10) (Punto 5)

DELTAS = [1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]   # espesor relativo           (Punto 4)
EPS_COL = 0.01                                    # ruido relativo en Exp. 4
EPS_PLOT = 0.1                                    # ruido usado en figuras de puntos


# =====================================================================
# SECCION A. NUCLEO Y UTILIDADES
# =====================================================================
def procrustes(X, Y, solo_rotacion=False):
    """Alinea X con Y: devuelve Q (ortogonal), b (traslacion) y S (valores singulares de M).
    solo_rotacion=False -> resuelve sobre O(2);  True -> sobre SO(2) (det Q = +1)."""
    mx, my = X.mean(0), Y.mean(0)
    Xc, Yc = X - mx, Y - my                       # 1. centrar
    M = Xc.T @ Yc                                 # 2. matriz de relacion
    U, S, Vt = np.linalg.svd(M)                   # 3. SVD
    D = np.eye(M.shape[0])
    if solo_rotacion:                             # 4. evitar reflexion
        D[-1, -1] = np.sign(np.linalg.det(U @ Vt))
    Q = U @ D @ Vt
    b = my - mx @ Q                               # 5. traslacion
    return Q, b, S


def rot(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


Q0 = rot(THETA0)                                  # rotacion verdadera (det = +1)
Q0_REF = Q0 @ np.diag([1.0, -1.0])                # reflexion verdadera (det = -1)


def rng_for(*key):
    """Generador independiente y reproducible para cada (experimento, nivel, repeticion)."""
    return np.random.default_rng([SEED, *key])


def make_noise(r, Xc, eps):
    """Ruido E con ||E||_F = eps * ||Xc||_F exacto."""
    Z = r.normal(size=Xc.shape)
    return eps * np.linalg.norm(Xc) * Z / np.linalg.norm(Z)


def cost(Xc, Yc, Q):
    """Costo cuadrado ||Yc - Xc Q||_F^2 sobre datos centrados."""
    return np.linalg.norm(Yc - Xc @ Q) ** 2


def angle_deg(Q, Qt):
    """Error angular en grados (nan si los determinantes difieren)."""
    if np.linalg.det(Q) * np.linalg.det(Qt) < 0:
        return np.nan
    A = Qt.T @ Q
    return abs(np.degrees(np.arctan2(A[1, 0], A[0, 0])))


def metrics(X, Y, Q, b, Qt, bt):
    """Separa RESIDUO DE AJUSTE (contra Y ruidoso) de ERROR DE RECUPERACION
    (contra la transformacion verdadera)."""
    fit = Y - (X @ Q + b)                         # residuo de ajuste
    true = (X @ Qt + bt) - (X @ Q + b)            # error contra la senal limpia
    return dict(
        q_err=np.linalg.norm(Q - Qt),             # recuperacion de Q (Frobenius)
        ang_err=angle_deg(Q, Qt),                 # recuperacion angular (grados)
        b_err=np.linalg.norm(b - bt),             # recuperacion de b
        rmse_fit=np.linalg.norm(fit) / np.sqrt(len(X)),     # ajuste
        rmse_true=np.linalg.norm(true) / np.sqrt(len(X)),   # vs senal limpia
    )


# =====================================================================
# SECCION F. ESTADISTICA Y SALIDA                         (Punto 5)
# =====================================================================
def agg(rows):
    """Promedio y desviacion estandar (ddof=1) de cada metrica sobre las repeticiones."""
    out = {}
    for k in rows[0]:
        v = np.array([r[k] for r in rows], float)
        out[k + "_mean"] = np.nanmean(v) if not np.all(np.isnan(v)) else np.nan
        out[k + "_std"] = np.nanstd(v, ddof=1) if np.sum(~np.isnan(v)) > 1 else np.nan
    return out


def save_csv(name, rows):
    with open(name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


log = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.append(s)


# =====================================================================
# SECCION B. CONSTRUCCION DE LOS PUNTOS                   (Punto 1)
# =====================================================================
def make_points():
    """N >= 20 puntos asimetricos (sesgo + curvatura) y no colineales."""
    r = rng_for(0)
    X = r.normal(size=(N, 2)) * np.array([4.0, 1.5])
    X[:, 0] += r.exponential(2.0, N)              # sesgo -> asimetria
    X[:, 1] += 0.15 * X[:, 0] ** 2 / 4            # curvatura -> sin simetria de espejo
    return X


def report_points(X, Xc):
    sv = np.linalg.svd(Xc, compute_uv=False)
    skew = np.mean(Xc[:, 0] ** 3) / np.mean(Xc[:, 0] ** 2) ** 1.5
    P("=== PARAMETROS ===")
    P(f"SEED={SEED}  N={N}  THETA0={np.degrees(THETA0):.1f} deg  B0={B0.tolist()}  REPS={REPS}")
    P(f"NOISE (||E||_F/||Xc||_F)={NOISE}")
    P(f"DELTAS={DELTAS}  EPS_COL={EPS_COL}")
    P(f"Valores singulares de Xc: {sv.round(4).tolist()}  (cociente s2/s1={sv[1]/sv[0]:.3f}) -> no colineal")
    P(f"Asimetria (skewness) de x: {skew:.3f}")
    P("Q0 =", Q0.round(6).tolist())
    P("Q0_ref (reflexion) =", Q0_REF.round(6).tolist(), " det=", round(np.linalg.det(Q0_REF), 6))
    np.savetxt("puntos_X.csv", X, delimiter=",", header="x,y", comments="")


# =====================================================================
# SECCION C. ROTACION + NIVELES DE RUIDO                  (Punto 2)
#   Y = X Q0 + 1 b0^T + E ; recuperacion de Q, b y alineamiento
#   Cada nivel se repite REPS veces                       (Punto 5)
# =====================================================================
def exp2_rotacion(X, Xc):
    P("\n=== EXP 1-2: Y = X Q0 + 1 b0^T + E (rotacion), SO(2) ===")
    rows, summ = [], []
    for li, eps in enumerate(NOISE):
        per = []
        for rep in range(REPS):
            r = rng_for(2, li, rep)
            E = make_noise(r, Xc, eps)
            Y = X @ Q0 + B0 + E
            Q, b, _ = procrustes(X, Y, solo_rotacion=True)
            m = metrics(X, Y, Q, b, Q0, B0)
            per.append(m)
            rows.append(dict(eps=eps, rep=rep, **m))
        a = agg(per)
        summ.append(dict(eps=eps, **a))
        P(f"eps={eps:5.2f} | ||Q-Q0||={a['q_err_mean']:.2e}±{a['q_err_std']:.1e}"
          f" | ang={a['ang_err_mean']:.3f}±{a['ang_err_std']:.3f} deg"
          f" | ||b-b0||={a['b_err_mean']:.2e}±{a['b_err_std']:.1e}"
          f" | RMSE ajuste={a['rmse_fit_mean']:.3e}±{a['rmse_fit_std']:.1e}"
          f" | RMSE vs limpia={a['rmse_true_mean']:.3e}±{a['rmse_true_std']:.1e}")
    save_csv("exp2_rotacion_corridas.csv", rows)
    save_csv("exp2_rotacion_resumen.csv", summ)
    return summ


# =====================================================================
# SECCION D. REFLEXION: O(2) vs SO(2)                     (Punto 3)
#   Verifica: costo_SO2 - costo_O2 = 4 * sigma2(M)  (datos centrados, sin escala)
# =====================================================================
def exp3_reflexion(X, Xc):
    P("\n=== EXP 3: Y = X Q0_ref + 1 b0^T + E (reflexion): O(2) vs SO(2) ===")
    rows, summ = [], []
    for li, eps in enumerate(NOISE):
        for rep in range(REPS):
            r = rng_for(3, li, rep)
            E = make_noise(r, Xc, eps)
            Y = X @ Q0_REF + B0 + E
            Yc = Y - Y.mean(0)
            QO, bO, S = procrustes(X, Y, False)
            QS, bS, _ = procrustes(X, Y, True)
            cO, cS = cost(Xc, Yc, QO), cost(Xc, Yc, QS)
            diff = cS - cO
            mO = metrics(X, Y, QO, bO, Q0_REF, B0)
            mS = metrics(X, Y, QS, bS, Q0_REF, B0)
            rows.append(dict(
                eps=eps, rep=rep, sigma2=S[1], cost_O2=cO, cost_SO2=cS,
                diff=diff, cuatro_sigma2=4 * S[1], err_abs=abs(diff - 4 * S[1]),
                detQ_O2=np.linalg.det(QO), q_err_O2=mO["q_err"], q_err_SO2=mS["q_err"],
                rmse_fit_O2=mO["rmse_fit"], rmse_fit_SO2=mS["rmse_fit"],
                b_err_O2=mO["b_err"], b_err_SO2=mS["b_err"]))

        sub = [x for x in rows if x["eps"] == eps]
        mean = lambda k: np.mean([s[k] for s in sub])
        std = lambda k: np.std([s[k] for s in sub], ddof=1)
        summ.append(dict(
            eps=eps,
            **{k + "_mean": mean(k) for k in
               ["diff", "cuatro_sigma2", "err_abs", "q_err_O2", "q_err_SO2",
                "rmse_fit_O2", "rmse_fit_SO2", "b_err_O2", "b_err_SO2"]},
            **{k + "_std": std(k) for k in
               ["q_err_O2", "q_err_SO2", "rmse_fit_O2", "rmse_fit_SO2"]}))
        P(f"eps={eps:5.2f} | det(Q_O2) medio={mean('detQ_O2'):+.0f}"
          f" | costo_SO2-costo_O2={mean('diff'):.5e} | 4*sigma2={mean('cuatro_sigma2'):.5e}"
          f" | |dif-4s2| max={max(s['err_abs'] for s in sub):.1e}"
          f" | ||Q-Q0|| O2={mean('q_err_O2'):.2e}, SO2={mean('q_err_SO2'):.2e}"
          f" | RMSE ajuste O2={mean('rmse_fit_O2'):.3e}, SO2={mean('rmse_fit_SO2'):.3e}")
    save_csv("exp3_reflexion_corridas.csv", rows)
    save_csv("exp3_reflexion_resumen.csv", summ)
    P("Maximo |(costo_SO2-costo_O2) - 4*sigma2| sobre todas las corridas:",
      f"{max(x['err_abs'] for x in rows):.3e}")

    # Control: si la verdad es una rotacion, ambos costos coinciden (diferencia 0)
    r = rng_for(30)
    Y = X @ Q0 + B0 + make_noise(r, Xc, 0.1)
    Yc = Y - Y.mean(0)
    QO, _, _ = procrustes(X, Y, False)
    QS, _, _ = procrustes(X, Y, True)
    P(f"Control (verdad=rotacion, eps=0.1): costo_SO2-costo_O2 = "
      f"{cost(Xc, Yc, QS) - cost(Xc, Yc, QO):.2e}")
    return summ


# =====================================================================
# SECCION E. CONFIGURACIONES CASI COLINEALES              (Punto 4)
#   X = [t, delta*w] ; delta -> 0 aplana la nube sobre una recta.
#   Se separan:
#     - residuo de ajuste       (rmse_fit,  contra Y ruidoso)
#     - error de recuperacion   (q_err, rmse_true, contra la verdad conocida)
#     - sensibilidad de Q       (perturb = ||Q(Y) - Q(Y limpia)||)
# =====================================================================
def make_collinear_base():
    rr = rng_for(4)
    t = rr.uniform(-5, 5, N) + rr.exponential(1.5, N)
    w = rr.normal(size=N)
    return t, w


def exp4_colineal(t, w):
    P("\n=== EXP 4: casi colinealidad (verdad = rotacion), ruido relativo %.3f ===" % EPS_COL)
    rows, summ = [], []
    for li, d in enumerate(DELTAS):
        Xd = np.column_stack([t, d * w])
        Xdc = Xd - Xd.mean(0)
        s_x = np.linalg.svd(Xdc, compute_uv=False)
        Yclean = Xd @ Q0 + B0
        QS_clean, _, _ = procrustes(Xd, Yclean, True)
        per = []
        for rep in range(REPS):
            r = rng_for(5, li, rep)
            Y = Yclean + make_noise(r, Xdc, EPS_COL)
            QS, bS, S = procrustes(Xd, Y, True)
            QO, bO, _ = procrustes(Xd, Y, False)
            mS = metrics(Xd, Y, QS, bS, Q0, B0)
            mO = metrics(Xd, Y, QO, bO, Q0, B0)
            per.append(dict(
                q_err_SO2=mS["q_err"], q_err_O2=mO["q_err"],
                flip_O2=float(np.linalg.det(QO) < 0),
                perturb_SO2=np.linalg.norm(QS - QS_clean),
                rmse_fit_SO2=mS["rmse_fit"], rmse_true_SO2=mS["rmse_true"],
                sigma2_M=S[1], cuatro_sigma2_M=4 * S[1]))
            rows.append(dict(delta=d, rep=rep, s2_s1_X=s_x[1] / s_x[0], **per[-1]))
        a = agg(per)
        summ.append(dict(delta=d, s2_s1_X=s_x[1] / s_x[0], **a))
        P(f"delta={d:6.3f} | s2/s1(X)={s_x[1]/s_x[0]:.4f}"
          f" | ||Q-Q0|| SO2={a['q_err_SO2_mean']:.2e}±{a['q_err_SO2_std']:.1e}"
          f" | O2={a['q_err_O2_mean']:.2e}±{a['q_err_O2_std']:.1e}"
          f" | frac. reflexiones O2={a['flip_O2_mean']:.2f}"
          f" | ||Q(Y)-Q(Y_limpia)||={a['perturb_SO2_mean']:.2e}"
          f" | RMSE ajuste={a['rmse_fit_SO2_mean']:.3e}"
          f" | RMSE vs limpia={a['rmse_true_SO2_mean']:.3e}"
          f" | 4*sigma2(M)={a['cuatro_sigma2_M_mean']:.3e}")
    save_csv("exp4_colineal_corridas.csv", rows)
    save_csv("exp4_colineal_resumen.csv", summ)
    return summ


# =====================================================================
# SECCION G. FIGURAS                                      (Punto 5)
#   Puntos antes/despues con ejes de IGUAL escala (set_aspect("equal")).
# =====================================================================
def eq(ax, title):
    ax.set_aspect("equal", adjustable="datalim")      # igual escala en x e y
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)


def fig1_puntos(X, Xc):
    """Rotacion (fila 1) y reflexion (fila 2): antes y despues."""
    li = NOISE.index(EPS_PLOT)
    r = rng_for(2, li, 0)
    Yr = X @ Q0 + B0 + make_noise(r, Xc, EPS_PLOT)
    Qa, ba, _ = procrustes(X, Yr, True)
    r = rng_for(3, li, 0)
    Yf = X @ Q0_REF + B0 + make_noise(r, Xc, EPS_PLOT)
    QOf, bOf, _ = procrustes(X, Yf, False)
    QSf, bSf, _ = procrustes(X, Yf, True)

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    ax[0, 0].scatter(*X.T, c="C0", label="X"); ax[0, 0].scatter(*Yr.T, c="C3", label="Y"); ax[0, 0].legend()
    eq(ax[0, 0], f"Rotacion, antes (eps={EPS_PLOT})")
    ax[0, 1].scatter(*Yr.T, c="C3", label="Y"); ax[0, 1].scatter(*(X @ Qa + ba).T, c="C2", marker="x", label="XQ+b")
    ax[0, 1].legend(); eq(ax[0, 1], "Rotacion, despues (SO(2))")
    ax[0, 2].scatter(*Yr.T, c="C3", label="Y")
    ax[0, 2].scatter(*(X @ Q0 + B0).T, c="k", marker="+", label="X Q0+b0 (verdad)")
    ax[0, 2].legend(); eq(ax[0, 2], "Y ruidoso vs senal limpia")
    ax[1, 0].scatter(*X.T, c="C0", label="X"); ax[1, 0].scatter(*Yf.T, c="C3", label="Y (reflejado)"); ax[1, 0].legend()
    eq(ax[1, 0], f"Reflexion, antes (eps={EPS_PLOT})")
    ax[1, 1].scatter(*Yf.T, c="C3", label="Y"); ax[1, 1].scatter(*(X @ QOf + bOf).T, c="C2", marker="x", label="O(2)")
    ax[1, 1].legend(); eq(ax[1, 1], "Reflexion, despues O(2) (det=-1)")
    ax[1, 2].scatter(*Yf.T, c="C3", label="Y"); ax[1, 2].scatter(*(X @ QSf + bSf).T, c="C1", marker="x", label="SO(2)")
    ax[1, 2].legend(); eq(ax[1, 2], "Reflexion, despues SO(2) (peor ajuste)")
    plt.tight_layout(); plt.savefig("fig1_puntos_antes_despues.png", dpi=130); plt.close()


def fig2_errores_vs_ruido(summ2, summ3):
    """Errores (media ± sd sobre REPS repeticiones) vs nivel de ruido."""
    pos = [s for s in summ2 if s["eps"] > 0]
    e = [s["eps"] for s in pos]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for k, lab in [("q_err", "||Q-Q0||_F"), ("b_err", "||b-b0||"),
                   ("rmse_fit", "RMSE ajuste (vs Y)"), ("rmse_true", "RMSE vs senal limpia")]:
        ax[0].errorbar(e, [s[k + "_mean"] for s in pos], [s[k + "_std"] for s in pos],
                       marker="o", capsize=3, label=lab)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("ruido relativo eps"); ax[0].legend()
    ax[0].set_title(f"Exp. 2: errores vs ruido (media±sd, {REPS} reps)"); ax[0].grid(alpha=.3, which="both")

    p3 = [s for s in summ3 if s["eps"] > 0]
    e3 = [s["eps"] for s in p3]
    ax[1].errorbar(e3, [s["q_err_O2_mean"] for s in p3], [s["q_err_O2_std"] for s in p3],
                   marker="o", capsize=3, label="O(2)")
    ax[1].errorbar(e3, [s["q_err_SO2_mean"] for s in p3], [s["q_err_SO2_std"] for s in p3],
                   marker="s", capsize=3, label="SO(2)")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("ruido relativo eps"); ax[1].legend()
    ax[1].set_title("Exp. 3 (verdad = reflexion): ||Q-Q0_ref||_F"); ax[1].grid(alpha=.3, which="both")
    plt.tight_layout(); plt.savefig("fig2_errores_vs_ruido.png", dpi=130); plt.close()


def fig3_casi_colineal(summ4):
    dd = [s["delta"] for s in summ4]
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    for k, lab in [("q_err_SO2", "||Q-Q0|| SO(2)"), ("q_err_O2", "||Q-Q0|| O(2)"),
                   ("perturb_SO2", "||Q(Y)-Q(Y limpia)|| SO(2)")]:
        ax[0].errorbar(dd, [s[k + "_mean"] for s in summ4], [s[k + "_std"] for s in summ4],
                       marker="o", capsize=3, label=lab)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].invert_xaxis(); ax[0].legend()
    ax[0].set_xlabel("delta (espesor)"); ax[0].set_title("Exp. 4: error de recuperacion"); ax[0].grid(alpha=.3, which="both")
    for k, lab in [("rmse_fit_SO2", "RMSE ajuste (vs Y)"), ("rmse_true_SO2", "RMSE vs senal limpia")]:
        ax[1].errorbar(dd, [s[k + "_mean"] for s in summ4], [s[k + "_std"] for s in summ4],
                       marker="o", capsize=3, label=lab)
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].invert_xaxis(); ax[1].legend()
    ax[1].set_xlabel("delta (espesor)"); ax[1].set_title("Residuo de ajuste vs error de transformacion")
    ax[1].grid(alpha=.3, which="both")
    ax[2].plot(dd, [s["flip_O2_mean"] for s in summ4], "o-")
    ax[2].set_xscale("log"); ax[2].invert_xaxis(); ax[2].set_ylim(-.05, 1.05)
    ax[2].set_xlabel("delta (espesor)"); ax[2].set_title("Fraccion de corridas con reflexion en O(2)")
    ax[2].grid(alpha=.3)
    plt.tight_layout(); plt.savefig("fig3_casi_colineal.png", dpi=130); plt.close()


def fig4_colineal_puntos(t, w, d=0.01):
    Xd = np.column_stack([t, d * w])
    r = rng_for(5, DELTAS.index(d), 0)
    Yd = Xd @ Q0 + B0 + make_noise(r, Xd - Xd.mean(0), EPS_COL)
    Qd, bd, _ = procrustes(Xd, Yd, True)
    QOd, bOd, _ = procrustes(Xd, Yd, False)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].scatter(*Xd.T, c="C0", s=15, label="X"); ax[0].scatter(*Yd.T, c="C3", s=15, label="Y"); ax[0].legend()
    eq(ax[0], f"Casi colineal (delta={d}), antes")
    ax[1].scatter(*Yd.T, c="C3", s=15, label="Y")
    ax[1].scatter(*(Xd @ Qd + bd).T, c="C2", marker="x", s=25, label="SO(2)")
    ax[1].scatter(*(Xd @ QOd + bOd).T, c="C1", marker="+", s=25, label=f"O(2) (det={np.linalg.det(QOd):+.0f})")
    ax[1].legend(); eq(ax[1], "Despues del alineamiento")
    plt.tight_layout(); plt.savefig("fig4_colineal_puntos.png", dpi=130); plt.close()


# =====================================================================
# SECCION H. PROGRAMA PRINCIPAL
# =====================================================================
def verificaciones(X, Xc):
    P("\n=== VERIFICACIONES (SciPy, ortogonalidad, formula del minimo) ===")
    r = rng_for(2, NOISE.index(0.1), 0)
    Y = X @ Q0 + B0 + make_noise(r, Xc, 0.1)
    Yc = Y - Y.mean(0)
    Q, b, S = procrustes(X, Y, solo_rotacion=False)
    R, escala = orthogonal_procrustes(Xc, Yc)  # SciPy no centra: se le pasan datos centrados
    P(f"||Q_propio - R_scipy||_F = {np.linalg.norm(Q - R):.2e}")
    P(f"scale SciPy = {escala:.10f} ; suma de valores singulares = {S.sum():.10f}")
    P(f"||Q^T Q - I||_F = {np.linalg.norm(Q.T @ Q - np.eye(2)):.2e}")
    c_directo = cost(Xc, Yc, Q)
    c_formula = np.linalg.norm(Xc) ** 2 + np.linalg.norm(Yc) ** 2 - 2 * S.sum()
    P(f"costo directo = {c_directo:.12f} ; formula = {c_formula:.12f}")
    P(f"|directo - formula| = {abs(c_directo - c_formula):.2e}")

def main():
    # Punto 1
    X = make_points()
    Xc = X - X.mean(0)
    report_points(X, Xc)

    # Puntos 2, 3, 4 (cada uno con REPS repeticiones y media ± sd: Punto 5)
    summ2 = exp2_rotacion(X, Xc)
    summ3 = exp3_reflexion(X, Xc)
    t, w = make_collinear_base()
    summ4 = exp4_colineal(t, w)
    verificaciones(X, Xc)

    # Punto 5: figuras
    fig1_puntos(X, Xc)
    fig2_errores_vs_ruido(summ2, summ3)
    fig3_casi_colineal(summ4)
    fig4_colineal_puntos(t, w)

    with open("resultados.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(log))
    print("\nListo. Archivos generados.")


if __name__ == "__main__":
    main()
