
import nfwutils
omega_ms = np.arange(0.15, 0.39, 0.005)
ws = np.arange(-1.76, -0.2, 0.005)
zclusters = [0.25, 1.0]
zsource = 1.5
c = [(.9,.6,0), (.35, .7, .9), (0,.6,.5), (0.95, 0.9, 0.25)]
        
    

scalings = []
for zcluster in zclusters:
    scaling = np.zeros_like(omega_ms)
    for i, omega_m in enumerate(omega_ms):
        cosmo = nfwutils.Cosmology(omega_m = omega_m, omega_l = 1. - omega_m)
        dl = cosmo.angulardist(zcluster)
        ds = cosmo.angulardist(zsource)
        dls = cosmo.angulardist(zcluster, zsource)
        scaling[i] = dl*ds/dls
    scalings.append(scaling)
    

hydro_scalings = []
mgas_scalings = []
for zcluster in zclusters:
    hydro_scaling = np.zeros_like(omega_ms)
    mgas_scaling = np.zeros_like(omega_ms)
    for i, omega_m in enumerate(omega_ms):
        cosmo = nfwutils.Cosmology(omega_m = omega_m, omega_l = 1. - omega_m)
        dl = cosmo.angulardist(zcluster)
        hydro_scaling[i] = dl
        mgas_scaling[i] = dl**(5/2)
    hydro_scalings.append(hydro_scaling)
    mgas_scalings.append(mgas_scaling)

omega_m_norm = np.arange(len(omega_ms))[np.abs(omega_ms - 0.3) < 0.001]

plot(omega_ms, scalings[0]/scalings[0][omega_m_norm], color=c[0], label=r'Lensing M($<R^X_{2500,ref}$)', linewidth=2)
plot(omega_ms, (scalings[0]/scalings[0][omega_m_norm])/(hydro_scalings[0]/hydro_scalings[0][omega_m_norm]), color=c[1], label=r'Lensing / Hydrostatic X-ray', linewidth=2)
plot(omega_ms, (scalings[0]/scalings[0][omega_m_norm])/(mgas_scalings[0]/mgas_scalings[0][omega_m_norm]), color=c[2], label=r'Lensing / M-gas X-ray', linewidth=2.)
plot(omega_ms, scalings[1]/scalings[1][omega_m_norm], color=c[0], linestyle='-.', linewidth=2)
plot(omega_ms, (scalings[1]/scalings[1][omega_m_norm])/(mgas_scalings[1]/mgas_scalings[1][omega_m_norm]), color=c[2], linestyle='-.', linewidth=2.)
plot(omega_ms, (scalings[1]/scalings[1][omega_m_norm])/(hydro_scalings[1]/hydro_scalings[1][omega_m_norm]), color=c[1], linestyle='-.', linewidth=2.)
#axis([0.15, 0.39, 0.7, 1.15])
axis([0.15, 0.39, 0.75, 1.4])
legend(loc='upper left')
xlabel(r'$\Omega_\mathrm{m}$', fontsize=20)
ylabel('Relative Scaling', fontsize=18)
tight_layout()
savefig('cosmo_relativescaling_om.png')
savefig('cosmo_relativescaling_om.pdf')
savefig('cosmo_relativescaling_om.eps')



scalings = []
hydro_scalings = []
mgas_scalings = []
for zcluster in zclusters:
    hydro_scaling = np.zeros_like(ws)
    mgas_scaling = np.zeros_like(ws)
    for i, w in enumerate(ws):
        cosmo = nfwutils.Cosmology(omega_m = 0.3, omega_l = 0.7, w=w)
        dl = cosmo.angulardist(zcluster)
        hydro_scaling[i] = dl
        mgas_scaling[i] = dl**(5/2)
    hydro_scalings.append(hydro_scaling)
    mgas_scalings.append(mgas_scaling)
    
for zcluster in zclusters:
    scaling = np.zeros_like(ws)
    for i, w in enumerate(ws):
        cosmo = nfwutils.Cosmology(omega_m = 0.3, omega_l = 0.7, w = w)
        dl = cosmo.angulardist(zcluster)
        ds = cosmo.angulardist(zsource)
        dls = cosmo.angulardist(zcluster, zsource)
        scaling[i] = dl*ds/dls
    scalings.append(scaling)

w_norm = np.arange(len(ws))[np.abs(ws + 1) < 0.001]
    
figure()
plot(ws, scalings[0]/scalings[0][w_norm], color=c[0], label=r'Lensing M($<R^X_{2500,ref}$)', linewidth=2)
plot(ws, scalings[1]/scalings[1][w_norm], color=c[0], linestyle='-.', linewidth=2)
plot(ws, (scalings[0]/scalings[0][w_norm])/(hydro_scalings[0]/hydro_scalings[0][w_norm]), color=c[1], label=r'Lensing / Hydrostatic X-ray', linewidth=2)
plot(ws, (scalings[1]/scalings[1][w_norm])/(hydro_scalings[1]/hydro_scalings[1][w_norm]), color=c[1], linestyle='-.', linewidth=2.)
plot(ws, (scalings[0]/scalings[0][w_norm])/(mgas_scalings[0]/mgas_scalings[0][w_norm]), color=c[2], label=r'Lensing / M-gas X-ray', linewidth=2)
plot(ws, (scalings[1]/scalings[1][w_norm])/(mgas_scalings[1]/mgas_scalings[1][w_norm]), color=c[2], linestyle='-.', linewidth=2)
#axis([-1.76, -0.2, 0.8, 1.4])
axis([-1.76, -0.2, 0.75, 1.4])
legend(loc='upper left')
xlabel(r'$\Omega_\mathrm{m}$', fontsize=20)
ylabel('Relative Scaling', fontsize=18)
tight_layout()
xlabel(r'$w$', fontsize=20)
savefig('cosmo_relativescaling_w.png')
savefig('cosmo_relativescaling_w.pdf')
savefig('cosmo_relativescaling_w.eps')
