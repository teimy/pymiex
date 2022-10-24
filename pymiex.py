import numpy as np
import miex_function


def lognormal_dist(mu, sigma):
    """ Returns lognormal distribution in 
    sense that Mischenko meant it
    Returns density function and effective radius
    """
    r = np.linspace(mu/10,mu*10,1000)
    
    pdf = (np.exp(-(np.log(r) - np.log(mu))**2 / (2 * np.log(sigma)**2)))/r
    r_eff = np.trapz(pdf*r*np.pi*r**2)/np.trapz(pdf*np.pi*r**2)
    v_eff = np.trapz(pdf*(r-r_eff)**2*np.pi*r**2)/np.trapz(pdf*np.pi*r**2)/r_eff**2
    return pdf, r_eff, v_eff

def wrap_mie(wavelength, n, k, radius, mode='single', scat_matrix=0, radmin=None, radmax=None, 
            refmed=1, nrad=100, alpha=None, nang=None, abundance=[100], sigma=None):
    nlam = len(wavelength)
    if scat_matrix is 0:
        nang=1
    elif nang is None:
        raise ValueError('Nang should not be None if scat_matrix is set to 1.')
    elif nang%2 != 1:
        raise ValueError('Nang should int and be odd number.')
    if len(n.shape):
        ncomp = 1
        n = n[None,:]
        k = k[None,:]
    else:
        ncomp = n.shape[0]
    if mode=='single':
        # A lot of parameters here are dummy and hardcoded just to pass some values to Fortran 
        albedo, gscatt, qext, qsca, qabs, qbk, qpr, cext, csca, cabs, cbk, cpr = miex_function.mie_routines.mie(
            nlam=nlam, ncomp=ncomp, wvl=wavelength, n=n, k=k, mode=mode, scat_matrix=scat_matrix,
            refmed=refmed, radmin=1, radmax=1000, alpha=0.7,
            r_g=radius, sigma_g=0.7, nrad=1, nang2=nang, abun=abundance)
        r_eff = radius
        v_eff = None

    if mode=='lognormal':
        # A lot of parameters here are dummy and hardcoded just to pass some values to Fortran 
        if radmin is None:
            raise ValueError('Please specify radmin if the mode is lognornal or power')
        if radmax is None:
            raise ValueError('Please specify radmax if the mode is lognornal or power')
        if sigma is None:
            raise ValueError('Please specify sigma if the mode is lognornal')

        albedo, gscatt, qext, qsca, qabs, qbk, qpr, cext, csca, cabs, cbk, cpr = miex_function.mie_routines.mie(
            nlam=nlam, ncomp=ncomp, wvl=wavelength, n=n, k=k, mode=mode, scat_matrix=scat_matrix,
            refmed=refmed, radmin=radmin, radmax=radmax, alpha=0.7,
            r_g=radius, sigma_g=sigma, nrad=nrad, nang2=nang, abun=abundance)
        pdf, r_eff, v_eff = lognormal_dist(radius, sigma)
        
    if mode=='power':
        # A lot of parameters here are dummy and hardcoded just to pass some values to Fortran 
        if radmin is None:
            raise ValueError('Please specify radmin if the mode is lognornal or power')
        if radmax is None:
            raise ValueError('Please specify radmax if the mode is lognornal or power')
        if alpha is None:
            raise ValueError('Please specify alpha if the mode is power')

        r_eff = None
        v_eff = None

        albedo, gscatt, qext, qsca, qabs, qbk, qpr, cext, csca, cabs, cbk, cpr = miex_function.mie_routines.mie(
            nlam=nlam, ncomp=ncomp, wvl=wavelength, n=n, k=k, mode=mode, scat_matrix=scat_matrix,
            refmed=refmed, radmin=radmin, radmax=radmax, alpha=alpha,
            r_g=radius, sigma_g=0.7, nrad=nrad, nang2=nang, abun=abundance)    
       
    return {'albedo':albedo, 'g':gscatt, 'qext':qext, 'qsca':qsca, 'qabs':qabs, 'qbk':qbk,
       'cext':cext, 'csca':csca, 'cabs':cabs, 'cbk':cbk, 'cpr':cpr, 'r_eff':r_eff, 'v_eff':v_eff}