"""
Sabra shell model: explicit transfer matrices and pole structure of the
spectral determinant Lambda_u(s,t).

Model:  du_n/dt = B_n(u) - nu k_n^2 u_n + f_n,   n = 0..N-1
        B_n(u) = i( a k_{n+1} u_{n+2} conj(u_{n+1})
                  + b k_n     u_{n+1} conj(u_{n-1})
                  - c k_{n-1} u_{n-1} u_{n-2} )
with k_n = k0 lam^n, a+b+c = 0 (energy conservation).
"""
import json
import numpy as np

rng = np.random.default_rng(7)

N    = 22
lam  = 2.0
k0   = 1.0
k    = k0 * lam ** np.arange(N)
a, b, c = 1.0, -0.5, -0.5           # a+b+c = 0
nu   = 1e-7
f    = np.zeros(N, complex)
f[1] = 0.1 * (1.0 + 1.0j)            # constant large-scale forcing
dt   = 1e-4

def nonlinear(u):
    U = np.zeros(N + 4, complex)
    U[2:-2] = u
    up1, up2 = U[3:-1], U[4:]
    um1, um2 = U[1:-3], U[0:-4]
    return 1j * (a * lam * k * up2 * np.conj(up1)
                 + b * k * up1 * np.conj(um1)
                 - (c / lam) * k * um1 * um2)

def rhs_nl(u):                       # nonlinearity + forcing
    return nonlinear(u) + f

E  = np.exp(-nu * k**2 * dt)
E2 = np.exp(-nu * k**2 * dt / 2)

def step(u):
    k1 = dt * rhs_nl(u)
    k2 = dt * rhs_nl(E2 * (u + k1 / 2))
    k3 = dt * rhs_nl(E2 * u + k2 / 2)
    k4 = dt * rhs_nl(E * u + E2 * k3)
    return E * u + (E * k1 + 2 * E2 * (k2 + k3) + k4) / 6
  
def check_energy_conservation():
    u = 0.05 * (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * k**(-1/3.)
    u[12:] = 0.0                     
    def rk4_inviscid(u, h):
        k1 = h * nonlinear(u)
        k2 = h * nonlinear(u + k1 / 2)
        k3 = h * nonlinear(u + k2 / 2)
        k4 = h * nonlinear(u + k3)
        return u + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    E0 = 0.5 * np.sum(np.abs(u)**2)
    h = 1e-5                          
    for _ in range(int(0.05 / h)):    
        u = rk4_inviscid(u, h)        
    E1 = 0.5 * np.sum(np.abs(u)**2)
    return abs(E1 - E0) / E0

def Rmat(w):   # multiplication v -> w v      (C-linear)
    return np.array([[w.real, -w.imag], [w.imag, w.real]])

def Kmat(w):   # map v -> w conj(v)           (antilinear)
    return np.array([[w.real,  w.imag], [w.imag, -w.real]])

def jacobian(u):
    U = np.zeros(N + 4, complex)
    U[2:-2] = u
    J = np.zeros((2 * N, 2 * N))
    for n in range(N):
        up1, up2 = U[n + 3], U[n + 4]
        um1, um2 = U[n + 1], U[n]
        kp1, kn_, km1 = lam * k[n], k[n], k[n] / lam
        def put(m, blk):
            if 0 <= m < N:
                J[2*n:2*n+2, 2*m:2*m+2] += blk
        put(n + 2, Rmat(1j * a * kp1 * np.conj(up1)))
        put(n + 1, Kmat(1j * a * kp1 * up2) + Rmat(1j * b * kn_ * np.conj(um1)))
        put(n - 1, Kmat(1j * b * kn_ * up1) + Rmat(-1j * c * km1 * um2))
        put(n - 2, Rmat(-1j * c * km1 * um1))
    return J

def check_jacobian():
    u = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * k**(-1/3.)
    x0 = np.empty(2 * N); x0[0::2] = u.real; x0[1::2] = u.imag
    def F(x):
        uu = x[0::2] + 1j * x[1::2]
        Bv = nonlinear(uu)
        out = np.empty(2 * N); out[0::2] = Bv.real; out[1::2] = Bv.imag
        return out
    Jfd = np.zeros((2 * N, 2 * N))
    for j in range(2 * N):
        h = 1e-6 * max(1.0, abs(x0[j]))
        e = np.zeros(2 * N); e[j] = h
        Jfd[:, j] = (F(x0 + e) - F(x0 - e)) / (2 * h)
    Jan = jacobian(u)
    denom = np.max(np.abs(Jan)) or 1.0
    return np.max(np.abs(Jan - Jfd)) / denom

def main():
    res = {}
    res["check_energy_drift"] = float(check_energy_conservation())
    res["check_jacobian_relerr"] = float(check_jacobian())
    print("energy drift (inviscid): %.3e" % res["check_energy_drift"])
    print("jacobian rel err:        %.3e" % res["check_jacobian_relerr"])
    assert res["check_energy_drift"] < 1e-7
    assert res["check_jacobian_relerr"] < 1e-6

    ph = np.exp(2j * np.pi * rng.random(N))
    u = 0.3 * k**(-1/3.) * ph
    u[17:] = 0.0
    T_transient, T_sample, sample_every = 60.0, 200.0, 0.2
    n_tr = int(T_transient / dt)
    for i in range(n_tr):
        u = step(u)

    n_sm = int(T_sample / dt)
    every = int(sample_every / dt)
    snaps, diss_acc, inp_acc, cnt = [], 0.0, 0.0, 0
    for i in range(n_sm):
        u = step(u)
        if i % 50 == 0:
            diss_acc += nu * np.sum(k**2 * np.abs(u)**2)
            inp_acc  += np.sum((f * np.conj(u)).real)
            cnt += 1
        if i % every == 0:
            snaps.append(u.copy())
    snaps = np.array(snaps)
    res["n_snapshots"] = len(snaps)
    res["eps_diss"]  = float(diss_acc / cnt)
    res["eps_input"] = float(inp_acc / cnt)
    res["balance_ratio"] = res["eps_input"] / res["eps_diss"]
    print("input/dissipation ratio: %.4f" % res["balance_ratio"])

    amp = np.sqrt(np.mean(np.abs(snaps)**2, axis=0))          # rms |u_n|
    eps = res["eps_diss"]
    k_eta = (eps / nu**3) ** 0.25
    res["k_eta"] = float(k_eta); res["n_eta"] = float(np.log2(k_eta / k0))
    n_lo, n_hi = 3, int(res["n_eta"]) - 3                      
    slope = np.polyfit(np.log(k[n_lo:n_hi]), np.log(amp[n_lo:n_hi]), 1)[0]
    res["inertial_range"] = [n_lo, n_hi]
    res["spectrum_slope"] = float(slope)                      
    print("spectrum slope: %.4f (K41: -0.3333)" % slope)

    def analyze(r):
        ns = np.arange(r, N - r)
        d = 2 * (2 * r + 1)
        nsn, rho, opn, eigs_by_n = [], [], [], {int(n): [] for n in ns}
        knun = []
        for s_ in snaps:
            J = jacobian(s_)
            for n in ns:
                A = J[2*(n-r):2*(n+r+1), 2*(n-r):2*(n+r+1)]
                ev = np.linalg.eigvals(A)
                rr = np.max(np.abs(ev))
                nsn.append(n); rho.append(rr)
                opn.append(np.linalg.norm(A, 2))
                knun.append(k[n] * abs(s_[n]))
                eigs_by_n[int(n)].append(ev)
        return (np.array(nsn), np.array(rho), np.array(opn),
                np.array(knun), eigs_by_n, d)

    out = {}
    for r in (1, 2):
        nsn, rho, opn, knun, eigs_by_n, d = analyze(r)
        lk = np.log(k[nsn])
        sigma = np.log(np.maximum(rho, 1e-300)) / lk          
        tau   = np.log(np.maximum(opn, 1e-300)) / lk         
        per_n = {}
        for n in sorted(set(nsn.tolist())):
            m = nsn == n
            per_n[int(n)] = dict(
                sigma_mean=float(np.mean(sigma[m])), sigma_std=float(np.std(sigma[m])),
                tau_mean=float(np.mean(tau[m])),
                ratio_mean=float(np.mean(rho[m] / opn[m])),
                ratio_min=float(np.min(rho[m] / opn[m])),
            )

        inert = (nsn >= n_lo) & (nsn <= n_hi)

        sh = np.array(sorted(set(nsn.tolist())))
        mlr = np.array([np.mean(np.log(rho[nsn == n])) for n in sh])
        mln = np.array([np.mean(np.log(opn[nsn == n])) for n in sh])
        sel = (sh >= n_lo) & (sh <= n_hi)
        slope_rho  = np.polyfit(np.log(k[sh[sel]]), mlr[sel], 1)
        slope_norm = np.polyfit(np.log(k[sh[sel]]), mln[sel], 1)
        out[r] = dict(
            dim=d, per_n=per_n,
            shells=sh.tolist(), mean_log_rho=mlr.tolist(), mean_log_norm=mln.tolist(),
            slope_rho=float(slope_rho[0]), slope_norm=float(slope_norm[0]),
            sigma_inertial_mean=float(np.mean(sigma[inert])),
            sigma_inertial_std=float(np.std(sigma[inert])),
            ratio_inertial_min=float(np.min((rho / opn)[inert])),
            ratio_inertial_q05=float(np.quantile((rho / opn)[inert], 0.05)),
            ratio_inertial_mean=float(np.mean((rho / opn)[inert])),
        )

        x = np.log(knun[inert]); y = np.log(np.maximum(rho[inert], 1e-300))
        theta, b0 = np.polyfit(x, y, 1)
        out[r]["lb_fit_theta"] = float(theta)
        out[r]["lb_ratio_min"] = float(np.min(rho[inert] / knun[inert]))
        out[r]["lb_ratio_mean"] = float(np.mean(rho[inert] / knun[inert]))
        out[r]["lb_ratio_q05"] = float(np.quantile(rho[inert] / knun[inert], 0.05))
        if r == 1:
            np.savez("./results/poles_data.npz",
                     nsn=nsn, rho=rho, opn=opn, knun=knun,
                     amp=amp, k=k, inert=inert, sigma=sigma,
                     snaps_last=snaps[-1],
                     **{f"eig_{n}": np.concatenate(eigs_by_n[n]) for n in eigs_by_n})
    res["windows"] = out
    with open("./results/results.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({kk: res[kk] for kk in
                      ["balance_ratio", "spectrum_slope", "n_eta"]}, indent=1))
    for r in (1, 2):
        w = res["windows"][r]
        print(f"r={r}: sigma_inertial = {w['sigma_inertial_mean']:.4f} "
              f"+- {w['sigma_inertial_std']:.4f}; "
              f"rho/||A|| mean {w['ratio_inertial_mean']:.3f} min {w['ratio_inertial_min']:.3f}; "
              f"theta fit {w['lb_fit_theta']:.3f}; "
              f"rho/(k|u|) min {w['lb_ratio_min']:.3f}; "
              f"slope_rho {w['slope_rho']:.3f} slope_norm {w['slope_norm']:.3f}")

if __name__ == "__main__":
    main()
