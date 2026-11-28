import numpy as np, json
from scipy.optimize import brentq
out = {}

# ---------- Corollary 1c: level vs log cross-partials under CES ----------
def mu_ces(R, G, kappa=1.0, xi=0.5, r=0.0):
    if abs(r) < 1e-9:
        return kappa * R**xi * G**(1 - xi)
    return kappa * (xi * R**r + (1 - xi) * G**r) ** (1 / r)

def cross(f, R, G, h=1e-4):
    return (f(R+h, G+h) - f(R+h, G-h) - f(R-h, G+h) + f(R-h, G-h)) / (4*h*h)

ok_lvl = ok_log = tot = 0
for sigma in np.linspace(0.05, 1.5, 60):
    r = 1 - 1/sigma
    for rho in np.linspace(0.02, 0.98, 60):
        R = G = 1.0
        mu = mu_ces(R, G, r=r)
        lam = rho * mu
        E = lambda R_, G_: 1.0 / (mu_ces(R_, G_, r=r) - lam)
        lE = lambda R_, G_: -np.log(mu_ces(R_, G_, r=r) - lam)
        c1, c2 = cross(E, R, G), cross(lE, R, G)
        b1, b2 = 1 - 2*sigma, 1 - sigma
        if abs(rho - b1) > 0.02 and abs(rho - b2) > 0.02:
            tot += 1
            ok_lvl += ((c1 < 0) == (rho < b1))
            ok_log += ((c2 < 0) == (rho < b2))
out['cor1c'] = dict(points=tot, level_match=ok_lvl, log_match=ok_log)

# ---------- Relaxed A5: R lowers execution time x0(R) ----------
lam_e, kappa, xi = 8.0, 3.0, 0.5
d0, d1, m, w = 1.0, 0.25, 0.5, 14.0
x0 = lambda R: 0.5 + 5.0 / R
def parts(R, G):
    mu = kappa * np.sqrt(R * G)
    tri = 1/(mu - lam_e) if mu > lam_e else np.inf
    c = 14 / G; zeta = 0.15 * (G - 1)
    dec = (1 - zeta)*(c/2 + d0) + zeta*d1
    enf = x0(R) + m*w/2
    return tri, dec, enf
def floor(G):
    c = 14 / G; zeta = 0.15 * (G - 1)
    return (1 - zeta)*(c/2 + d0) + zeta*d1 + 0.5 + m*w/2
out['floor'] = {g: round(floor(g), 2) for g in (1.5, 4.5)}
out['E_at'] = {f'G{g}_R{R}': round(sum(parts(R, g)), 2) for g in (1.5, 4.5) for R in (5, 10, 25)}
out['cadence_share_floor'] = {g: round(((1-0.15*(g-1))*(14/g)/2 + m*w/2)/floor(g), 3) for g in (1.5, 4.5)}

# ---------- Calibration consistent with Result 2 ----------
def crit_mean(a, rho, mu, s):
    W0 = rho / mu; rho1 = a*s*rho
    return 1/mu + a*W0/(1-rho1) + (1-a)*W0/((1-rho1)*(1-rho))
s = 0.15
bound = {}
for sm in (1, 2, 5):
    mu = 1/sm
    bound[sm] = round(1/mu + (1/mu)/(1 - s), 2)  # a=1, rho->1
out['prio_bound_days'] = bound
calib = []
for target_name, total in (('All KEV', 174), ('Critical KEV', 137)):
    resid = total - 26.0  # c=14,w=30,zeta=0,d0=x0=2,m=1, tau_pre=0
    for sm in (1, 2, 5):
        mu = 1/sm
        for rho in (0.90, 0.95, 0.99):
            f = lambda a: crit_mean(a, rho, mu, s) - resid
            if f(0) < 0:
                calib.append([target_name, sm, rho, None]); continue
            a = brentq(f, 0, 1)
            calib.append([target_name, sm, rho, round(1-a, 2)])
out['calib_misclass'] = calib

# ---------- Corollary 5a: speed and stability co-move across G ----------
phi0, thetaLx, Lc = 0.5, 1.0, 60.0
c5 = []
for G in (1.5, 3.0, 4.5):
    T = 7.5 / G
    tau = T*np.log(phi0*Lc/(T*thetaLx))
    c5.append([G, round(T, 2), round(tau, 2), round(phi0*np.exp(-tau/T), 3), round(thetaLx*tau + Lc*phi0*np.exp(-tau/T), 2)])
out['cor5a'] = c5
json.dump(out, open('/home/claude/analysis.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))

# ---------- Figure 4 regenerated ----------
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), dpi=200)
sg = np.linspace(0.02, 1.2, 400)
a0 = ax[0]
a0.fill_between(sg, 0, np.clip(1-2*sg, 0, 1), color='#2a6fb0', alpha=.35, lw=0, label='Complements in level and log')
a0.fill_between(sg, np.clip(1-2*sg, 0, 1), np.clip(1-sg, 0, 1), color='#2a6fb0', alpha=.12, lw=0, label='Complements in log only')
a0.fill_between(sg, np.clip(1-sg, 0, 1), 1, color='#e8743b', alpha=.18, lw=0, label='Substitutes in both')
a0.plot(sg, 1-2*sg, 'k-', lw=1.8); a0.plot(sg, 1-sg, 'k--', lw=1.5)
a0.text(0.20, 0.50, r'$\rho = 1-2\sigma$ (level, Eq. 6b)', rotation=-55, fontsize=9, ha='center', va='center', bbox=dict(fc='white', ec='none', alpha=.7))
a0.text(0.62, 0.38, r'$\rho = 1-\sigma$ (log, Eq. 6c)', rotation=-36, fontsize=9, ha='center', va='center', bbox=dict(fc='white', ec='none', alpha=.7))
a0.axvline(1.0, color='grey', ls=':', lw=1); a0.text(1.01, 0.05, 'Cobb–\nDouglas', fontsize=8, color='grey')
a0.set_xlim(0.02, 1.2); a0.set_ylim(0, 1)
a0.set_xlabel(r'Elasticity of substitution $\sigma$ between R and G'); a0.set_ylabel(r'Utilization $\rho$')
a0.set_title('(a) Sign of the R–G cross-partial of triage latency', loc='left')
a0.legend(loc='upper right', fontsize=8, frameon=False)
a1 = ax[1]
Rs = np.linspace(2, 25, 400)
for g, col, ls in ((1.5, '#e8743b', '--'), (4.5, '#2a6fb0', '-')):
    y = [sum(parts(R, g)) for R in Rs]
    a1.plot(Rs, y, ls, color=col, lw=2.2, label=f'G = {g}')
    a1.axhline(floor(g), color=col, ls=':', lw=1)
    a1.text(25, floor(g)-0.9, f'floor F(G={g}) = {floor(g):.1f} d', ha='right', fontsize=8, color='#444')
a1.set_ylim(0, 20); a1.set_xlabel('Security resources R (analyst FTE-equivalent)'); a1.set_ylabel(r'E[$\tau$] (days)')
a1.set_title('(b) Resources lower latency only towards F(G)', loc='left')
a1.legend(frameon=False); a1.grid(alpha=.3)
for a in ax:
    a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
plt.tight_layout(); plt.savefig('/home/claude/fig4_new.png')
