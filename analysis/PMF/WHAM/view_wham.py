import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt

import wham_make_file as wm
import WhamPotentialClass as WPC



if __name__ == "__main__":
        
    nrg2eV=.002782
    eV2kcalmol=23.06 
    eVkb=8.617e-5
    
    spring_constant = 30.0 * nrg2eV * eV2kcalmol
    # spring constant in kCal/mol Angstrom
    fl = 'H.txt'
    # file name
    conv_unit = 3.96
    # unit conversion
    delta_res = 0.2
    # resolution
    temperature=330
    # temperature of simulation 
    kbt = temperature * eVkb* eV2kcalmol
    # k_B T
    #n_wind = 18
    n_wind=len(wm.handle_R('R.txt'))
    # number of windows
    xtol = 1e-8
    # tolerance
    coeff = spring_constant/kbt
    
    # read the data
    data = wm.make_data(fl, conv_unit)
    
    # make bins
    h_fun, r_val = wm.make_histograms(data, delta_res)
    
    # define the parameters for the spring potential
    #params = (coeff*np.ones(n_wind), np.linspace(59, 59-(n_wind-1)*2, num=n_wind))
    params = (coeff*np.ones(n_wind), wm.handle_R('R.txt'))
    
    # initial guess
    ini_conf = 0.8*np.ones(n_wind)
    """ 
    # create the class
    wham = WPC.WhamPotential(r_val, h_fun, params)
    print wham.hist
    print wham.r_0
    # compute the potential with the initial guess
    wham.fp_solve(ini_conf, xtol)
    # use the roux solution
    
    # this is rho
    rho_val = wham.rho()
    
    # roux solution CMP 91 (1995) 275
    #sol_roux = np.array([0.00, -1.027, -0.316, -0.745, -0.467])    
        
    # this is the value of the free energy
    f_i = kbt*np.log(wham.soln)
    for ind in range(n_wind):
        txt = '(FREE ENERGY) %1.4f' % (f_i[ind]-f_i[0])
        print(txt)
    
    # plot the potential with the same units as Roux
    pf = wham.pmf()   
    val_min = np.amin(pf)
    plt.plot(wham.r, kbt*(pf-val_min))
    
    plt.show()
    #compute the integral
    int_norm = np.zeros(n_wind)
    for ind_p in range(n_wind):
        int_norm[ind_p] = trapz(h_fun[ind_p], wham.r)
        
    # compare the data with the fitted data
    for ind in range(n_wind):
        plt.plot(wham.r, h_fun[ind]/int_norm[ind])
        plt.plot(wham.r, wham.h_comp(ind))
        txt =  'computed/measured bin %d' % (ind +1)
        plt.title(txt)
        plt.show()


    val = np.linspace(wham.der_pos[0], wham.der_pos[-1], 40)
    print(val)
    print(wham.eval_der(val))
    print('\n')
    print(wham.der_pos)
    print(wham.der)
    pmf_ap = np.zeros_like(val)
    for ind, r in enumerate(val):
        pmf_val = wham.pmf_approx(r)[0]
        pmf_ap[ind] = pmf_val
    
    const=26
    val_min=np.amin(pmf_val)
    plt.plot(val, kbt*(pmf_ap-val_min)+const)
    plt.plot(wham.r, kbt*(pf-val_min))
    print('\n')    
    print(pmf_ap)
    plt.show()
    
    wham.write_soln()
    """ 
    wham = WPC.WhamPotential(r_val, h_fun, params, solution='F.txt')
    
     # this is rho
    rho_val = wham.rho()
    
    # roux solution CMP 91 (1995) 275
    #sol_roux = np.array([0.00, -1.027, -0.316, -0.745, -0.467])    
        
    # this is the value of the free energy
    f_i = kbt*np.log(wham.soln)
    for ind in range(n_wind):
        txt = '(FREE ENERGY) %1.4f' % (f_i[ind]-f_i[0])
        print(txt)
    
    # plot the potential with the same units as Roux
    pf = wham.pmf()   
    val_min = np.amin(pf)
    plt.plot(wham.r, kbt*(pf-val_min))
    
    plt.show()
    #compute the integral
    int_norm = np.zeros(n_wind)
    for ind_p in range(n_wind):
        int_norm[ind_p] = trapz(h_fun[ind_p], wham.r)
        
    # compare the data with the fitted data
    for ind in range(n_wind):
        plt.plot(wham.r, h_fun[ind]/int_norm[ind])
        plt.plot(wham.r, wham.h_comp(ind))
        txt =  'computed/measured bin %d' % (ind +1)
        plt.title(txt)
        plt.show()
