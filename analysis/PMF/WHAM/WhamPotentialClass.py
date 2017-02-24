""":module:: WhamPotentialClass
   :platform: Unix, Windows
   :synopsis: computes the potential

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November 2016
"""

import numpy as np
from scipy.integrate import trapz, quad
from scipy import optimize

class WhamPotential(object):
    """Class to compute wham potentials"""
    
    def __init__(self, r_val, hf, params,solution=False):
        """
        The constructor
        :param r_val: edge r values
        :param hf: data file (a numpy array with dimensions [num_of_windows, number of points in histogram])
        :param params: parameters of the harmonic term
        :param solution: name of file storing free energy solutions; defaults to False
        """

        self.num_windows = hf.shape[0]
        # number of windows
        self.num_points = np.sum(hf, axis=1)
        # number of points on each window
        
        self.r = r_val[:-1] + 0.5*np.diff(r_val) 
        # define the mid-point of the interval
        
        self.hist = hf 
        # store the histogram
        
        # k and r_0 for the harmonic potential
        self.k, self.r_0 = params
        
        # store the solution for free energies here
        if(solution!=False):
            x=open(solution,'r')
            y=x.readlines()
            print len(y)
            print y
            if (len(y)!=self.num_windows):
                print("The given free energies and number of windows are not compatible. \n Not reading in free energies")
                self.soln = False
            else:
                self.soln=np.zeros(self.num_windows)
                for i in range(self.num_windows):
                    self.soln[i]=float(y[i])
        else:
            self.soln=False
        
        # compute derivatives
        arg = np.argmax(self.hist, axis=1)
        self.der_pos = self.r[arg]
        self.der = self.k*(self.r_0 - self.der_pos)
        self.reverse = False
        if self.r_0[0] > self.r_0[-1]:
            self.reverse = True

    def harm(self, ind):
        """
        harmonic potential
    
        :param ind: window number
        :return: value at harmonic potential
        :rtype: numpy array
        """

        val = (self.r-self.r_0[ind])*(self.r-self.r_0[ind])
    
        return np.exp(-0.5*self.k[ind]*val)


    def int_rho(self, x):
        """
        function :math:'$\\rho$'
    
        :param x: :math:'$\\exp(\\beta \\Delta F_i)$'
        :return: value at x
        :rtype: numpy array
        """
        
        norm = np.zeros_like(self.r)
    
        for ind in range(self.num_windows):
            norm += self.num_points[ind]*self.harm(ind)*x[ind]
         
        val = np.sum(self.hist, axis=0)/norm
    
        val_const = trapz(val, self.r)
        # normalize the result
        
        return val/val_const
    
    def w_int(self, xval, ind):
        """Helper function, convenient for direct iterations and checking solutions
        
        :param xval: x input
        :param ind: window number
        :return: integrand value
        :rtype: numpy array
        """
        
        return trapz(self.int_rho(xval)*self.harm(ind), self.r)
    
    def fp_solve(self, ini_conf, xs=1e-8):
        """Computes the fixed point free energy
        
        :param ini_conf: initial guess
        :param xs: tolerance
        """
    
        self.soln = np.zeros([self.num_windows])
        val_fun = np.zeros([self.num_windows])
        
        def loc_fun(x):
            for ind in range(self.num_windows):
                val_fun[ind] = 1.0/self.w_int(x, ind)
            return val_fun
        
        self.soln = optimize.fixed_point(loc_fun, ini_conf, xtol=xs,maxiter=10000)
    
    def rho(self):
        """Returns the rho function
        :return rho:
        :rtype: numpy array
        """
        
        if type(self.soln).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
            return 0
        
        return self.int_rho(self.soln)
    
    def pmf(self):
        """Returns the potential
        :return rho:
        :rtype: numpy array
        """
        
        if type(self.soln).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
            return 0
        
        return -np.log(self.rho())
    
    def h_comp(self, ind):    
        """Returns the computed biased distribution
        :ind: window value
        :return computed bias potential for bias:
        :rtype: numpy array
        """
        
        if type(self.soln).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
            return 0
        
        return self.harm(ind)*self.rho()*self.soln[ind]
    
    def prob(self, ind):
        """Returns the probability function for each window
        :return rho:
        :rtype: numpy array
        """
        
        if type(self.soln).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
            return 0
        
        # TO DO
        
        return 0
        
    def eval_der(self, r_val):
        """computes the derivatives at the minimum of each potential
            :param r_val: values at which the derivative is computed
            :returns: derivatives at points
            :rtype: numpy array
        """
        
        if self.reverse:
            return np.interp(r_val, np.flipud(self.der_pos), np.flipud(self.der))
        
        return np.interp(r_val, self.der_pos, self.der)
    
    def pmf_approx(self, rval):
        """Aproximate pmf
        
        :param rval: value of the potential
        :return: pmf
        :rtype: float
        """
        
        return quad(self.eval_der, self.der_pos[0], rval)
        
    def write_soln(self):
	"""outputs the free energy solutions to a texxt file so that they don't need
	   to be computed a second time when the results need to be seen later

	   returns True if the file is written
	"""
        if type(self.soln).__module__ != 'numpy':
            print('fixed point solution needs to be computed first')
            return False
        fid=open('F.txt','w')
        for x in range(self.num_windows):
            fid.write(str(self.soln[x]))
            fid.write('\n')
        return True
        
	

	
	













