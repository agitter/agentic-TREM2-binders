"""
Render each binder Cα trace as a 2D 'tube' image, following the pipeline
documented in py2Dmol's technical_readme.md (Rendering System):

  rotate -> orthographic project -> painter's-algorithm z-sort of segments
  -> grid shadow map (darken segments below the per-cell max z) -> outline.

py2Dmol itself draws to an HTML5 canvas via JavaScript, so it cannot rasterise
headlessly here; its Python-side best_view() is used directly for the canonical
orientation, and the draw stage is reimplemented on the same spec.
"""
import numpy as np, json, re, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from py2Dmol.viewer import best_view

CA = {k: np.array(v) for k, v in json.load(open('ca.json')).items()}
aln = pd.read_json('aln_sup.json')
aln['key'] = aln.q + '|' + aln.t
A = aln.set_index('key')


def _pairs(qs, ts, cig):
    qi, ti = qs - 1, ts - 1
    P = []
    for n, op in re.findall(r'(\d+)([MID])', cig):
        n = int(n)
        if op == 'M':
            P += [(qi + k, ti + k) for k in range(n)]
            qi += n; ti += n
        elif op == 'I':
            qi += n
        else:
            ti += n
    return P


def superpose(ref, mob):
    """Rigid-superpose mob onto ref using the Foldseek TM-align correspondence."""
    if ref == mob:
        return CA[ref].copy()
    k = ref + '|' + mob
    if k not in A.index:
        return None
    r = A.loc[k]
    r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
    P = _pairs(int(r.qs), int(r.ts), r.cigar)
    if len(P) < 4:
        return None
    X = CA[ref][[i for i, _ in P]]
    Y = CA[mob][[j for _, j in P]]
    H = (Y - Y.mean(0)).T @ (X - X.mean(0))
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (CA[mob] - Y.mean(0)) @ R.T + X.mean(0)


def _spline(P, k=6):
    """Catmull-Rom subdivision of the Cα trace for a smooth tube."""
    n = len(P)
    if n < 4:
        return P
    out = []
    for i in range(n - 1):
        p0 = P[max(i - 1, 0)]; p1 = P[i]; p2 = P[i + 1]; p3 = P[min(i + 2, n - 1)]
        for t in np.linspace(0, 1, k, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t +
                              (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
                              (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(P[-1])
    return np.array(out)


def render(coords, hue, px=340, width=13.0, pad=0.10, shadow=True, frame=None):
    """Draw one structure. Returns an HxWx4 RGBA array with transparent background.

    frame: optional (R, centre) from py2Dmol.best_view computed once on a shared
    reference. Passing one frame for every structure keeps already-superposed
    folds in a common orientation, so similar folds look similar on the page.
    """
    C = _spline(coords)
    R, ctr = best_view(C) if frame is None else frame
    V = (C - ctr) @ R.T

    # normalise into a unit box, keeping aspect
    span = np.ptp(V[:, :2], axis=0).max()
    V[:, :2] = (V[:, :2] - V[:, :2].mean(0)) / (span * (1 + pad))
    z = V[:, 2]
    z = (z - z.min()) / max(np.ptp(z), 1e-9)

    segs = np.stack([V[:-1, :2], V[1:, :2]], axis=1)
    zs = 0.5 * (z[:-1] + z[1:])

    # --- grid shadow map: darken segments that sit below the per-cell max depth
    depth_cue = zs.copy()
    if shadow:
        G = 26
        gx = np.clip(((segs[:, :, 0].mean(1) + .5) * G).astype(int), 0, G - 1)
        gy = np.clip(((segs[:, :, 1].mean(1) + .5) * G).astype(int), 0, G - 1)
        top = {}
        for i, (a, b) in enumerate(zip(gx, gy)):
            top[(a, b)] = max(top.get((a, b), -1), zs[i])
        occl = np.array([top[(a, b)] - zv for a, b, zv in zip(gx, gy, zs)])
        depth_cue = np.clip(zs - 0.55 * occl, 0, 1)

    # --- colour: one hue family per workflow, N->C light to dark,
    #     further modulated by depth so the fold reads in 3D
    base = np.array(matplotlib.colors.to_rgb(hue))
    t = np.linspace(0, 1, len(segs))
    light = 1.0 - 0.45 * t                       # N-term light -> C-term dark
    shade = 0.45 + 0.55 * depth_cue              # far dark -> near bright
    cols = np.clip(base[None, :] * (light * shade)[:, None] +
                   (1 - (light * shade))[:, None] * base[None, :] * 0.30, 0, 1)
    cols = np.clip(cols * 1.25, 0, 1)
    rgba = np.concatenate([cols, np.ones((len(cols), 1))], axis=1)

    order = np.argsort(zs)                       # painter's algorithm

    fig = plt.figure(figsize=(px / 100, px / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-.53, .53); ax.set_ylim(-.53, .53)
    ax.set_axis_off()
    fig.patch.set_alpha(0)

    # black outline underneath, then the coloured tube on top
    ax.add_collection(LineCollection(segs[order], colors='black',
                                     linewidths=width * 1.5, capstyle='round',
                                     joinstyle='round', zorder=1))
    ax.add_collection(LineCollection(segs[order], colors=rgba[order],
                                     linewidths=width, capstyle='round',
                                     joinstyle='round', zorder=2))

    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)

    # circular crop
    h, w = img.shape[:2]
    yy, xx = np.mgrid[0:h, 0:w]
    r = np.sqrt((xx - w / 2) ** 2 + (yy - h / 2) ** 2)
    img[..., 3] = np.where(r <= w / 2 - 1, img[..., 3], 0)
    return img
