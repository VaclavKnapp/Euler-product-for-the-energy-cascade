import numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 10, 'figure.dpi': 150})
d = np.load('./results/poles_data.npz')
res = json.load(open('./results/results.json'))
nsn, rho, opn, knun = d['nsn'], d['rho'], d['opn'], d['knun']
k, amp, sigma = d['k'], d['amp'], d['sigma']
N = len(k); n_lo, n_hi = res['inertial_range']
sh = np.array(sorted(set(nsn.tolist())))
lk = np.log(k)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(np.arange(N), np.log2(amp), 'o-', ms=4, color='k', label=r'$\log_2\,\mathrm{rms}\,|u_n|$')
nn = np.arange(n_lo, n_hi + 1)
sl = res['spectrum_slope']
ref = np.log2(amp[n_lo]) + sl * (nn - n_lo) / np.log(2) * np.log(2)
ax[0].plot(nn, np.log2(amp[n_lo]) + sl*(nn-n_lo), 'r--', lw=1,
           label=f'fit slope {sl:.3f}')
ax[0].plot(nn, np.log2(amp[n_lo]) - (nn-n_lo)/3, 'b:', lw=1, label='K41  $-1/3$')
ax[0].axvline(res['n_eta'], color='gray', lw=0.8, ls=':')
ax[0].text(res['n_eta']+0.2, ax[0].get_ylim()[0]+2, r'$n_\eta$', color='gray')
ax[0].set_xlabel('$n$'); ax[0].set_ylabel(r'$\log_2 \mathrm{rms}\,|u_n|$')
ax[0].legend(frameon=False, fontsize=8); ax[0].set_title('(a) energy spectrum')

mlr = np.array([np.mean(np.log(rho[nsn==n])) for n in sh])
mln = np.array([np.mean(np.log(opn[nsn==n])) for n in sh])
mlku = np.array([np.mean(np.log(knun[nsn==n])) for n in sh])
ax[1].plot(lk[sh]/np.log(2), mlr/np.log(2), 'o-', ms=4, label=r'$\langle\ln\rho(A_n)\rangle$')
ax[1].plot(lk[sh]/np.log(2), mln/np.log(2), 's-', ms=3, label=r'$\langle\ln\|A_n\|\rangle$')
ax[1].plot(lk[sh]/np.log(2), mlku/np.log(2), '^-', ms=3, label=r'$\langle\ln k_n|u_n|\rangle$')
sel = (sh>=n_lo)&(sh<=n_hi)
p = np.polyfit(lk[sh[sel]], mlr[sel], 1)
xx = lk[sh[sel]]
ax[1].plot(xx/np.log(2), np.polyval(p, xx)/np.log(2), 'r--', lw=1,
           label=f"fit slope {p[0]:.3f}")
ax[1].axvspan(n_lo, n_hi, color='0.93')
ax[1].set_xlabel(r'$\log_2 k_n$'); ax[1].set_ylabel(r'$\log_2(\cdot)$')
ax[1].legend(frameon=False, fontsize=8); ax[1].set_title('(b) cascade scaling ($r=1$)')
fig.tight_layout(); fig.savefig('./figures/figA.pdf'); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
cmap = plt.get_cmap('viridis')
for n in sh:
    if n < 2: continue
    ev = d[f'eig_{n}']
    ev = ev[np.abs(ev) > 1e-12]
    # one snapshot's worth: take the last dim eigenvalues
    dim = 6
    ev = ev[-dim:]
    s = np.log(ev.astype(complex)) / np.log(k[n])   # principal branch
    ax[0].plot(s.real, s.imag, '.', ms=4, color=cmap((n - 2) / (N - 3)))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(2, N - 1))
cb = fig.colorbar(sm, ax=ax[0]); cb.set_label('shell $n$')
ax[0].axvline(0, color='gray', lw=0.6)
ax[0].set_xlabel(r'$\mathrm{Re}\,s$'); ax[0].set_ylabel(r'$\mathrm{Im}\,s$ (principal branch)')
ax[0].set_title(r'(a) dynamical poles of $\Lambda_u(\cdot,t)$, one snapshot')

msig = np.array([np.mean(sigma[nsn==n]) for n in sh])
ssig = np.array([np.std(sigma[nsn==n]) for n in sh])
ax[1].errorbar(sh, msig, yerr=ssig, fmt='o-', ms=4, capsize=2)
ax[1].axhline(2/3, color='b', ls=':', lw=1, label='K41  $2/3$')
ax[1].axvline(res['n_eta'], color='gray', lw=0.8, ls=':')
ax[1].axvspan(n_lo, n_hi, color='0.93')
ax[1].set_xlabel('$n$'); ax[1].set_ylabel(r'$\sigma_n=\ln\rho(A_n)/\ln k_n$')
ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title('(b) pole abscissa per shell (mean $\\pm$ sd)')
fig.tight_layout(); fig.savefig('./figures/figB.pdf'); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ratio = rho / opn
mr = [np.mean(ratio[nsn==n]) for n in sh]
q05 = [np.quantile(ratio[nsn==n], 0.05) for n in sh]
mn = [np.min(ratio[nsn==n]) for n in sh]
ax[0].plot(sh, mr, 'o-', ms=4, label='mean')
ax[0].plot(sh, q05, 's--', ms=3, label='5% quantile')
ax[0].plot(sh, mn, 'v:', ms=3, label='min')
ax[0].axvspan(n_lo, n_hi, color='0.93')
ax[0].set_ylim(0, 1.05)
ax[0].set_xlabel('$n$'); ax[0].set_ylabel(r'$\rho(A_n)/\|A_n\|$')
ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title('(a) non-normality is bounded')

inert = d['inert']
x, y = np.log(knun[inert]), np.log(rho[inert])
hb = ax[1].hexbin(x/np.log(2), y/np.log(2), gridsize=45, cmap='Greys', mincnt=1, bins='log')
th, b0 = np.polyfit(x, y, 1)
xs = np.linspace(x.min(), x.max(), 10)
ax[1].plot(xs/np.log(2), (th*xs+b0)/np.log(2), 'r-', lw=1.2,
           label=fr'fit: $\theta={th:.2f}$')
cmin = res['windows']['1']['lb_ratio_min']
ax[1].plot(xs/np.log(2), (xs+np.log(cmin))/np.log(2), 'b--', lw=1.2,
           label=fr'$\rho={cmin:.2f}\,k_n|u_n|$ (empirical bound)')
ax[1].set_xlabel(r'$\log_2 k_n|u_n|$'); ax[1].set_ylabel(r'$\log_2 \rho(A_n)$')
ax[1].legend(frameon=False, fontsize=8, loc='upper left')
ax[1].set_title('(b) spectral lower bound (inertial range)')
fig.tight_layout(); fig.savefig('./figures/figC.pdf'); plt.close(fig)
print('figures done!')
