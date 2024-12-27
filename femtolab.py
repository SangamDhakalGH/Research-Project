# -*- coding: utf-8 -*-

# 2012.02.12 18:12:49 CET
"""
Created on Thu Feb  9 14:23:31 2012

@author: larry
"""
import sys
import os
from PyQt4 import QtGui
from pylab import *
#from scipyxxxxx.sparse.csgraph import _validation
from scipy.optimize import fmin_l_bfgs_b, curve_fit
from scipy import ndimage as ndimage
from scipy.integrate import odeint
from matplotlib.mlab import cross_from_below
from matplotlib.pyplot import subplots_adjust
from femtofunc import *
from time import *
import matplotlib as mp
from matplotlib import colorbar
import numpy as np

mp.rcParams['font.size'] = 8

import locale
locale.setlocale(locale.LC_NUMERIC, 'C')

class trspectrum:
    """Time-resolved spectrum"""

    def __init__(self):
        self.S = [] # the full signal matrix
        self.Sr= [] # the selected range ofthe signal matrix
        self.Sbk=[] # backup for undo
        self.t = [] # the time axis
        self.w = [] # the wavelength axis
        self.tr = [] # the time axis
        self.wr = [] # the wavelength axis
        self.c = [] # the concentrations of the photoexcitations
        self.Z = '0' # the Z value if it is part of a matrix
        self.L = '0' # the pump intensity if applicable
        self.S0=np.array([0, 0], dtype=float)       # starting concentration
        self.pfad=os.getcwd()
        self.fn='test_init.dat'
        self.units_energy = 'nm'
        self.units_time = 'fs'
        
        self.logt=0 # logarithmic representation of time axis?
        self.logS=0 # logarithmic representation of signal axis?
        self.P=[] # fixed parameters for the fitting sessions
        self.W=[] # position of variable and fixed params for the fitting sessions
        self.uF=[] # code variable, will contain the user function
        self.nst=[] # number of states in rate eqn model
        self.nwp=[] # number of values of parameter for compromise fit
        self.cwd=os.getcwd()
        '''Parameters for ODE solver'''
        self.tcrit=0.  # critical time, eg generation, CAN BE VECTOR for second pump
        self.h0=1.     # first step
        self.hmax=1000. # max step size allowd
        self.mxstep=1000 # max numb of steps for one time eval.
        
        self.fittype='spectral' # can be 'spectral', 'global_gaussian_eli'. Param to pass to StartFitSession()
        
    def read_array(self, filename, dtype, separator = ','):
        """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
        """
        N=[] #dummy
        cast = N.cast
        data = [ [] for dummy in xrange(len(dtype)) ]
        for line in open(filename, 'r'):
            fields = line.strip().split(separator)
            for (i, number,) in enumerate(fields):
                data[i].append(number)


        for i in xrange(len(dtype)):
            data[i] = cast[dtype[i]](data[i])

        return N.rec.array(data, dtype=dtype)



    def loadmatrix(self, fn, units_energy, units_time, fn_base =''):
        """
        Load a pump a probe matrix
        Extract Signal, time wavelength
        type=1: single matrix that contains wavelength as first col
        and time as first row
        type=2: matrix contains only signal (spektren in spalten, zeiten in zeilen)
        and time in extra file _timing.dat, wavelength in _wl.dat
        fn_base includes the full path!
        """
        # find out which type:
        try:
            loadtxt(fn[0:-4] + '_timing.dat')
            loadtxt(fn[0:-4] + '_wl.dat')
            typ=2 # matrix format with extra time and wavelength files found
        except:
            typ=1 # assume wavelengths and times are first col and row, resp. 
            
        if fn_base != '':
            typ = 3 # Milano 2D spectra
        
        
        if typ == 1:
            try:
                f = loadtxt(fn)
            except:
                f=loadtxt(fn, delimiter=',')
            t = f[0, 1:]
            w = f[1:, 0]
            S = f[1:, 1:]
            
        if typ == 2:
            S = loadtxt(fn)
            t = loadtxt(fn[0:-4] + '_timing.dat')
            w = loadtxt(fn[0:-4] + '_wl.dat')
            
        if typ == 3:
            S = loadtxt(fn)
            t = loadtxt(fn_base + '_X.dat')
            w = loadtxt(fn_base + '_Y.dat')
        
        if w[0]>w[-1]:
            w=w[::-1] # flip E
            S=flipud(S) # must flip also S
            
        self.S = S
        self.Sr = S
        self.t = t
        self.tr = t
        self.w = w
        self.wr = w
        
        self.units_energy = units_energy
        self.units_time = units_time
        pf,fn=pfad(fn)
        self.fn=fn
        self.pfad=pf
        nanpos=isnan(S)
        print( 'Full time axis from ' + str(t[0]) + ' fs to ' + str(t[-1]) + ' fs.')
        print( 'Full wavelength axis from ' + str(w[0]) + ' nm to ' + str(w[-1]) + ' nm.')
        print( 'Full energy axis from ' + str(10000000.0 / w[0] / 8065) + ' eV to ' + str(10000000.0 / w[-1] / 8065) + ' eV.')

    def loadpop(self,fn):
        '''load a matrix of the following convention:
        1st col time, other columns signal.
        a new trspectrum() instance is used, 
        '''
        try:
            f = loadtxt(fn)
        except:
            f=loadtxt(fn) # hä??? again?
        pf,fn=pfad(fn)
        self.fn=fn
        self.pfad=pf

        self.c = f[:, 1:]
        self.t = f[:, 0]
        self.w = array(range(self.c.shape[1]),dtype=float)*0.0+1.0 # yields ones
        self.S=self.c.transpose()
        self.S0=self.c[0,:] # the initial concentrations


    def append_matrix(self,fn='',p=[]):
        # appends a matrix to an existing one
        # expected: same units in time and wavelength
        # expected: no overlap along wavelength axis
        if p==[]:
            p2=trspectrum()
            p2.loadmatrix(fn, units_energy='nm', units_time='fs')
        else:
            p2 = p
        # interpolate new matrix to old matrix
        Si=zeros((len(p2.w),len(self.t)))
        for xx in range(len(p2.w)):
            Si[xx,:]=interp(self.t,p2.t,p2.S[xx,:])
        self.S = vstack((self.S, Si))
        self.w = hstack((self.w, p2.w))
        iw = argsort(self.w, axis=-1, kind='quick', order=None)
        self.w = self.w[iw]
        self.S = self.S[iw, :]

    def load_timetraces(self, fn, curd):
        """construct a matrix from time traces
        fn points to a file with the following columns
        todo:
            find way to read the first col directly(loadtxt cannot do this)
            this would avoid the writing of the second col
        """
        f1 = loadtxt(fn, usecols=(1, 2, 3))
        for row in f1:
            nr = row[0]
            nrs = str(int(nr))
            if nr < 10:
                nrs = '0' + nrs
            f1n = curd + nrs + '.dat'
            print( 'Loading ' + f1n)
            self.add_lockin_trace(self.pfad+ f1n, row[1], row[2] / 1000, 22)


    def TimeTrace(self, fn, wl,iosc,fattPD, mode='add'):
        """loads a single time trace and includes it into the matrix
        ToDo:
            customize units_energy and units_time
        """
        f1 = loadtxt(fn)
        t = f1[:, 0]
        tt = f1[:, 1]
        tt = tt / iosc * fattPD
        
        if mode=='new':
            tt=array(tt)
            tt=tt.reshape(1, -1)
            self.w = wl
            self.t = t
            self.S = tt
            
        if mode=='add':
            # need to interpolate to existing time trace
            tti=interp(self.t,t,tt)
            tti=array(tti)
            tti=tti.reshape(1,-1)
            self.S = vstack((self.S, tti))
            self.w = hstack((self.w, wl))
            iw = argsort(self.w, axis=-1, kind='quick', order=None)
            self.w = self.w[iw]
            self.S = self.S[iw, :]
        if self.fn==[]:
            self.fn=fn

    def Spectrum(self, fn, t,iosc,fattPD, mode='add'):
        """loads a single spectrum and includes it into the matrix
        ToDo:
            customize units_energy and units_time
        """
        f1 = loadtxt(fn, delimiter=',')
        w = f1[:, 0]
        tt = f1[:, 1]
        tt = tt / iosc * fattPD
        
        if mode=='new':
            tt=array(tt)
            tt=tt.reshape(-1, 1)
            self.w = w
            self.t = t
            self.S = tt

        if mode=='add':
            # need to interpolate to existing spectra
            tti=interp(self.w,w,tt)
            tti=array(tti)
            tti=tti.reshape(-1,1)
            self.S = hstack((self.S, tti))
            self.t = hstack((self.t, t))
            it = argsort(self.t, axis=-1, kind='quick', order=None)
            self.t = self.t[it]
            self.S = self.S[:,it]
        if self.fn==[]:
            self.fn=fn
    
    def neg_zero(self,startat=-2500,upto=-500):
        ''' sets all time traces to zero, by subtracting the average between 
        startat and upto'''
        print( 'setting negative regions to zero,')
        print( 'using the average from ' + str(startat) + ' to ' + str(upto) + ' fs.')
        for xx in range(len(self.w)):
            ee = nonzero(logical_and(self.t>=startat,self.t<=upto))
            mitte=mean(self.S[xx,ee])
            self.S[xx,:]=self.S[xx,:]-mitte
            print( 'subtracted: '+ str(mitte))

    def cutmatrix(self, wmin, wmax, tmin, tmax,mm='copy'):
        # cut the matrix S to a smaller wavelength and time range
        # mm='copy' generates a copy, not touching the original
        # new attributes Sr,tr, and wr are the created
        # mm='original cuts the original 
        fpe=1e-10 # work with a certain floating point error
        
        print(tmin)
        print(tmax)
        
        if tmin > tmax:
            parken = tmin *1.0
            tmin = tmax *1.0
            tmax = parken *1.0
            
        if wmin > wmax:
            parken = wmin *1.0
            wmin = wmax *1.0
            wmax = parken *1.0
        
        
        ti = [where((self.t >= tmin-tmin*fpe) & (self.t <= tmax+tmax*fpe))][0][0]
        wi = [where((self.w >= wmin-wmin*fpe) & (self.w <= wmax+wmax*fpe))][0][0]        
        
        
        iwn = wi[0]
        iwx = wi[-1]
        itn = ti[0]
        itx = ti[-1]
        if mm=='copy':
            self.tr = self.t[ti]
            self.wr = self.w[wi]
            self.Sr = self.S[iwn:(iwx + 1), itn:(itx + 1)]
        if mm=='original':
            self.t = self.t[ti]
            self.w = self.w[wi]
            self.S = self.S[iwn:(iwx + 1), itn:(itx + 1)]
        #print( self.Sr.shape

    def plotmatrix(self,wmin,wmax,tmin,tmax,canvas,fig,ax):
        """ 
        Plots a matrix in false color code
        """
        
        # check if collection already exists:
        self.cutmatrix(wmin, wmax, tmin, tmax, 'copy')
        (X, Y,) = np.meshgrid(self.tr, self.wr)

        locale.setlocale(locale.LC_NUMERIC, 'C')
        ax.clear()
        CS=ax.contourf(X, Y, self.Sr, 100)
        ax.set_title('Time-resolved spectrum')
        ax.set_ylabel('Probe wavelength (nm)')
        ax.set_xlabel('Pump-probe delay (fs)')
        
        # fig.colorbar(CS, use_gridspec = True, orientation = 'vertical')  # Learn how to use a colorbar without importing whole pylab
        canvas.draw()

    def plot_spectrum(self,t,canvas,ax):
        w=self.wr
        it=int(round(interp(t,self.tr,range(len(self.tr)))))
        S=self.Sr[:,it]
        ax.clear()
        ax.plot(w,S)
        ax.set_title('Single spectrum at t= ' + str(self.tr[it]))
        ax.set_ylabel('DT/T')
        ax.set_xlabel('Probe wavelength (nm)')
        canvas.draw()
    
    def plot_tt(self,w,canvas,ax):
        t=self.tr
        iw=round(interp(w,self.wr,range(len(self.wr))))
        S=self.Sr[iw,:]
        S=S.reshape(-1,)
        ax.clear()
        ax.plot(t,S)
        ax.set_title('Single time trace at w= ' + str(self.wr[iw]))
        ax.set_ylabel('DT/T')
        ax.set_xlabel('Pump-probe delay (fs)')
        canvas.draw()

    def savematrix(self,fn):
        """saves the trspectrum in standard format:
            1st row time in fs
            1st col wavelength in nm
            rest data in DT/T
        """
        tsr=concatenate(([0],self.t)) # save corresponding original spectra
        E=self.w.reshape(-1,1)                       # wavelength across rows
        raus=hstack((E,self.S))
        raus=vstack((tsr,raus))
        savetxt(fn,raus)
        pf,fn=pfad(fn)
        self.fn=fn
        self.pfad=pf

    def transponier(self):
        # transpose matrix and exchange w and t vectors
        self.S = self.S.transpose()
        w = self.t
        t = self.w
        self.t = t
        self.w = w
        if len(self.wr) !=0:
            self.Sr = self.Sr.transpose()
            wr = self.tr
            tr = self.wr
            self.tr = tr
            self.wr = wr

    def logbin(self, tmin, tmax, ncl,linbininit=False,axis=1):
        # perform a logarithmic binning on the data, including the linear initial part if requuired
        # axis =1 (bin t axis), axis = 0 (bin w axis)
        if axis == 0:
            self.transponier()
            
        tlog=logspace(log(tmin),log(tmax),num=ncl,base=exp(1))
        lowestbin=exp(0.5*(log(tlog[1])+log(tlog[0])))
        highestbin=exp(0.5*(log(tlog[-1])+log(tlog[-2])))
        binlimits=logspace(log(lowestbin),log(highestbin),num=ncl-1,base=exp(1))
        binlimits=concatenate(([tmin*1.001],binlimits,[tmax*0.999]))
        print(len(self.w))
        print(len(tlog))
        Slog=zeros((len(self.w),len(tlog)),dtype=float)
        print('logarithmic binning:')
        for xx in range(len(tlog)): # bin and average
            ib0=cross_from_below(self.t,binlimits[xx])
            ib1=cross_from_below(self.t,binlimits[xx+1])
            Slog[:,xx]=mean(self.S[:,ib0:ib1],axis=1)
            tlog[xx]=mean(self.t[ib0:ib1])
            print('t = '+ str(tlog[xx]) + ': bins ' + str(ib0) + ' - ' + str(ib1))
        
        if linbininit:
            tanf=self.t[self.t<tmin] # the iniitial linear part of the time axis
            Sanf=self.S[:,self.t<tmin]
            self.S=hstack((Sanf,Slog))
            self.t=hstack((tanf,tlog))   
        else:
            self.S=Slog
            self.t=tlog
        if axis == 0:
            self.transponier()
    
    def linbin(self,ncl,axis=1):
        # bins the data along the time axis, strongly reducing the amount of data
        # ncl is the reduction factor (assuming normal distribution, should improve the noise by sqrt(ncl)) 
        # axis: 1 = binning along t axis; 0: binning along w axis
        if axis == 0:
            self.transponier()
        t=self.t
        w=self.w
        S=self.S
        numpts=floor(len(t)/ncl)
        binlimits,dt=linspace(t[0],t[-1],numpts+1,retstep=True) # these are all the limits of the bins
        tlog=binlimits[:-2]+dt/2. # middle of the bins
        Slog=zeros((len(w),len(tlog)),dtype=float)
        Slog=self.bin_core(tlog,binlimits,S,Slog)
        self.S=Slog
        self.t=tlog
        if axis == 0:
            self.transponier()
    
    def smooth(self,smE,smt,equi=False):
        print(equi)
                
        self.Sbk=self.S+0. # backup to make operation undo-able
        
        if equi: # user has requested smoothing in real units, not pixels
            ts=linspace(self.t[0],self.t[-1],(self.t[-1]-self.t[0])/(self.t[1]-self.t[0])+1) # equidistant time axis
            S1=zeros((len(self.w),len(ts)),dtype=float)
            for xx in range(len(self.w)):
                S1[xx,:]=interp(ts,self.t,self.S[xx,:]) # S1 is array with equidistant time axis so a meaningful convolution is possible
            if smt > 1: # user wants this to be time units not pixels, so convert into pixels
                smt=int(round(smt/(self.t[1]-self.t[0])))
                print('Gaussian convolution over ' + str(smt) + ' pixels requested.')
            stype='gaussian'
        else:
            S1=self.S*1.0
            stype='mov_avg'
            
        if smE>1:
            for xx in range(S1.shape[1]):
                S1[:,xx]=smooth(S1[:,xx],smE)
        if smt>1:
            for xx in range(S1.shape[0]):
                S1[xx,:]=smooth(S1[xx,:],smt,stype=stype)
        
        if equi: # need to convert back into non equidistant time axis
            for xx in range(len(self.w)):
                self.S[xx,:]=interp(self.t,ts,S1[xx,:])
        else:
            self.S=S1

    
    def bin_core(self,tlog,binlimits,S,Slog):
        # the core of the binning function
        for xx in range(len(tlog)): # bin and average
            ib = [where((self.t >= binlimits[xx]) & (self.t <= binlimits[xx+1]))][0][0]
            #print( 't = '+ str(tlog[xx]) + ': bins ' + str(ib[0]) + ' - ' + str(ib[-1])
            Slog[:,xx]=mean(self.S[:,ib],axis=1)
        return Slog    

    def nm_to_eV(self):
        """change wavelength axis into eV"""
        self.w=1e7/self.w/8065.
        self.w=self.w[::-1] # flip E
        self.S=flipud(self.S) # must flip also S
        
        if len(self.wr)!=0:
            self.wr=1e7/self.wr/8065.
            self.wr=self.wr[::-1] # flip E
            self.Sr=flipud(self.Sr) # must flip also S
        
        self.units_energy='eV'
        return self
    
    def dtt_to_da(self):
        self.S = -log(self.S+1.)
        
    def dechirp_fromfile(self,pfad,fn):
        """loads two column file of arbitrary length
        two columns: wavelength/t0
        interpolates to get t0 for every line in matrix
        dechirps"""
        self.chirp = loadtxt(pfad + fn)
        self.dechirp_core()    
    
    def dechirp_fit(self,tstart=-1000,tstop=1000,sm=11,binw=12,mode='steepest rise'):
        """ fully automatic routine: for every single time trace
        smooth, make first derivative &find maximum in certain range"""
        print('automatic dechirping has been chosen')
        if binw==(-1):
            S=self.Sr*1.0
            w=self.wr*1.0
            t=self.tr*1.0
        else:
            S,w,t=binnen(self.S, self.w, self.t, binw)
        chirp=zeros((len(w),2),dtype=float)        
        if mode=='steepest rise':
            print('Time zero found by steepest rise')

            for xx in range(len(w)):
                print('dechirping at lambda = ' + str(w[xx]) + ' nm')
                gg=S[xx,:]
                ggs=smooth(gg,g=sm)
                ggd=abs(gradient(ggs))
                ggd=ggd[where(t<tstop)]
                tt=t[where(t<tstop)]
                ggd=ggd[where(tt>tstart)]
                tt=tt[where(tt>tstart)]
                im=argmax(ggd)
                print(im)
                print('found zero at ' + str(tt[im]) + ' fs.')
                chirp[xx,0]=w[xx]
                chirp[xx,1]=tt[im]
            
            if binw != (-1):
                print('fitting with second order polynomial...')
                z=polyfit(chirp[:,0],chirp[:,1],2)
                print(z)
                chirpfit=polyval(z,chirp[:,0])
                print('removing outliers')
                dy=array(chirpfit-chirp[:,1],dtype=float) # the deviations
                print(dy)
                outi,nout,madm=outlierx(dy,f=4,allbutout=1) # returns the index of the non-outliers 
                chirpc=chirp[outi,:]
                print('Number of outliers found: ' + str(nout))
                print('madm, for completeness: ' + str(madm))
                print('fitting without outliers')
                z=polyfit(chirpc[:,0],chirpc[:,1],2) # fit without outliers
                print(z)
                chirpfit2=polyval(z,chirpc[:,0])
                plot(chirp[:,0],chirpfit,chirp[:,0],chirp[:,1],'o',chirpc[:,0],chirpfit2,chirpc[:,0],chirpc[:,1],'*')
                title('Chirp function with and without outliers')
                show()
                
                #now do the dechirp:
                chirpfull=polyval(z,self.w)
                chirpout=zeros((len(self.w),2),dtype=float)
                chirpout[:,0]=self.w
                chirpout[:,1]=chirpfull
            else:
                chirpout=chirp[:] # don't assume smooth chirp -> rigorously dechiro 
        
        if mode=='biexp':
            print('Time zero found by convoluted biexponential')
            t1=200. # first decay time
            t2=2000. # second decay time
            b=170. # width of convolution
            step=1.1 # global y-stretch of pre-defined cubic chirp function
            offset=800
            p0=[t1,t2,b,step,offset]
            bd=[(None,None),
                (None,None),
                (None,None),
                (None,None),
                (None,None)]
            
            self.Sbk=self.S+0. # backup to make operation undo-able
            
            #plsq1,f,d=fmin_l_bfgs_b(self.errsq_dechirp,p0,fprime=None,args=(S,t,w),approx_grad=1,bounds=bd,
            #                        epsilon=0.2,factr=1e-22,pgtol=1e-18,
            #                        maxfun=400)
            self.ychirp=S
            
            S12d=S*gausi(t,1,5000,0)
            self.S1=S12d.reshape(-1,)
            
            plsq1,pcov=curve_fit(self.errsq_dechirp_wrap,w,self.S1,p0=p0,epsfcn=0.1,sigma=None,maxfev=2000,xtol=1e-5)
            
            print('result of fitting:')
            print(plsq1)
            #print( f
            #print( d            
            figure(2)
            plot(w,self.a[:,0],w,self.a[:,1])
            
            
            figure(3)
            (X, Y,) = np.meshgrid(self.tr, self.wr)
            locale.setlocale(locale.LC_NUMERIC, 'C')
            contourf(X, Y, self.res, 40)
            
            figure(5)
            (X, Y,) = np.meshgrid(self.tr, self.wr)
            locale.setlocale(locale.LC_NUMERIC, 'C')
            contourf(X, Y, self.Sfit, 40)
            
            L=zeros((len(p0)-3,),dtype=float)
            for xx in range(len(p0)-3):
                L[xx]=1/plsq1[xx]
            B=eigenvec_to_B(len(L),L)
            nseq=dot(self.a,inv(B.T))
            
            figure(4)
            plot(self.w,self.tc)
            
            figure(6)
            plot(w,nseq[:,0],w,nseq[:,1])

            show()
                
        outname=self.fn[0:-4]+'_chirp.dat' # append labed 'chirp' to original filename
        print(self.w.shape)
        chirpout=zeros((self.w.shape[0],2),dtype=float)
        chirpout[:,0]=self.w
        chirpout[:,1]=self.tc
        savetxt(self.pfad + outname,chirpout)
        
        save_matrix_milano_convention(range(nseq.shape[1]),w,nseq,self.pfad+self.fn,'nseq')
        
        self.chirp=chirpout
        self.dechirp_core()    
    
    def errsq_dechirp_wrap(self,x,*args):
        V=array(args) # the parameters that are varied 
        y=self.ychirp
        Sout=self.errsq_dechirp(V,y,x)
        return Sout

    def errsq_dechirp(self,V,S,w):
        nexp=2 # 1,2 oder 3 exponential?
        V=V.reshape(max(V.shape),)
        St=S*1.0
        print('Fitparameters:')
        print(V)
        t=self.tr 
        dt=t[1]-t[0]
        b=V[2]
        step=V[3]
        offset=V[4]
        
        ca=offset # overall offset 
        cb=-70.36103*step # linear term
        cc=-79.698*step # cubic term
        tc=ca+cb*w+cc*w**3 # prepare test chirp
        CI=zeros((nexp,S.shape[1]),dtype=float) # the c matrix
        a=zeros((S.shape[0],nexp),dtype=float) # the preexponential factors
        irflen=2
        tirf=linspace(-irflen*b,irflen*b,2*irflen*b/dt+1)
        irf=gausi(tirf,1,b,0)
        irf=irf/sum(irf)
        for ii in range(nexp): # prepare CI matrix of norm exponentials
            CI[ii,:]=exp(-1./V[ii]*t)
            CI[ii,t<0]=0
            CI[ii,:]=np.convolve(CI[ii,:], irf, 'same')
        
        for xx in range(len(w)): # apply test chirp
            St[xx,:]=interp(t,t-tc[xx],S[xx,:]) # St is S with test dechirp applied
            a[xx,:]=linalg.lstsq(CI.T, St[xx,:].T)[0] # pre-exponential factors of chirp-free fit with same dynamics
                
        Sfit=dot(a,CI) # chirp-free fit
        res=(Sfit-St)**2*gausi(t,1,5000,0)
        
        Sfitc=zeros_like(Sfit )
        for xx in range(len(w)): # apply test chirp to fitted data 
            Sfitc[xx,:]=interp(t,t+tc[xx],Sfit[xx,:]) # St is S with test dechirp applied
        
        #es=(Sfit-St)**2
        #es=sum(es,0)*gausi(t,1,2000,0)/Sfit.shape[0]
        
        es=sum(sum(res))/Sfit.shape[0]/Sfit.shape[1]
        print('error squared:' + str(es))
        self.tc=ca+cb*self.w+cc*self.w**3 # chirp function applied to full wavelength axis
        self.Sfit=Sfit
        
        Sfitc=Sfitc*gausi(t,1,5000,0)
        
        Sout=Sfitc.reshape(-1,)
        es_heimlich=sum((Sout-self.S1)**2)
        print('heimlicher Fehler: ' + str(es_heimlich))
        self.res=res
        self.a=a
        #return es*1e6
        return Sout # need 1D array for curve_fit
    
    def dechirp_core(self):
        wc=self.chirp[:,0]
        tc=self.chirp[:,1]
        t0=interp(self.w,wc,tc) # the interpolated chirp for all wls
        self.Sbk=self.S+0. # backup to make operation undo-able
        for xx in range(len(self.w)):
            self.S[xx,:]=interp(self.t,self.t-t0[xx],self.S[xx,:])
            
    def dechirp_step(self, step , r13 = 0.882):
        w = self.w
        ccalib = -79.698 # calibration of chirp in Madrid fs spectrometer system. 
        wm = 0.5 * (w[0] + w[-1]) # middle of the energy region
        cb = ccalib * r13 * step
        cc = ccalib * step
        ca = -cb * wm - cc * wm ** 3 # adjust offset so that overall shift of time axis is as small as possible
        tc = ca + cb * w + cc * w ** 3 # prepare test chirp
        self.chirp = zeros((len(self.w),2),dtype=float)
        self.chirp[:,0] = self.w
        self.chirp[:,1] = tc
        self.dechirp_core()
        
    def dropoutx(self,tmin,tmax):
        '''bad (too low) spectra are eliminated by detecting a derivative
        which is far above the median of the median.
        '''
        itmin=cross_from_below(self.t,tmin)
        itmax=cross_from_below(self.t,tmax)
        S1=diff(sum(self.S,0))
        [outlier, nout,madm]=outlierx(S1[itmin:itmax],5,0)
        
        # if they come as neighbors, keep only the higher of them
        rauswurf=[]
        for xx in range(nout-1):
            if outlier[xx+1]==outlier[xx]+1:
                rauswurf.append(xx)
        print(rauswurf)
        outlier=delete(outlier,rauswurf)
        
        msgtxt=''
        msgtxt=msgtxt+'Outlier found at: ' + str(self.t[outlier+itmin])
        msgtxt=msgtxt + '; Number of outliesr: ' +str( len(outlier))
        msgtxt=msgtxt + '. Proceed?'
        msc = QtGui.QMessageBox(None)
        msc.setText(msgtxt)
        msc.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
        msc.setDefaultButton(QtGui.QMessageBox.Ok)
        r=msc.exec_()
        msc.accept()
        if r==QtGui.QMessageBox.Cancel:
            self.abbruch=True
        if r==QtGui.QMessageBox.Ok:
            self.abbruch=False
        if self.abbruch==False:
            self.t=delete(self.t,outlier+itmin)
            self.S=delete(self.S,outlier+itmin, axis=1)
        
    def normal(self,ax,sm=1):
        # normalize matrix along axis ax. Find maximum in a smoothed
        # replica, that is smoothed by a moving average with sm points
        # but perform the normalization on the un-smoothed curve
        
        if ax==0:
            for xx in range(len(self.t)):
                S0=self.S[:,xx]
                Ssm=smooth(S0,sm)
                self.S[:,xx]=S0/max(abs(Ssm))
        elif ax==1:
            for xx in range(len(self.w)):
                S0=self.S[xx,:]
                Ssm=smooth(S0,sm)
                self.S[xx,:]=S0/max(abs(Ssm))
        return self
    
    def dS_dt(self,S,t=0):
        """the rate equation(s)"""
        #g=interp(t,self.t,self.g)
        ag=0
        bg=1
        t0=0
        ggg=0
        
        exec(self.uF) # execute the byte coded user-defined function
        # producing something like 'ggg=array([g-k1*S[0], g-k2*S[1]])'
        # return array([g-k1*S[0], g-k2*S[1]]) # use this if nothing goes
        return ggg
        
           
    def odesolver(self):
        S, infodict = odeint(self.dS_dt, self.S0part, self.tr, full_output=True,hmax=self.hmax,
                                      tcrit=self.tcrit, h0=self.h0, mxstep=self.mxstep)
        self.A=S    
        infodict['message'] # >>> 'Integration successful.

    
    def exp_fit(self,ts,te,plsq,fig,logbx=''):
        """exponential fitting routine for all wavelengths"""
        if logbx !='':(logbx.append('********************************'))
        if logbx !='':(logbx.append('Separate exponential fits of all time traces'))
        
        it=[where((self.t >= ts) & (self.t <= te))][0][0]
        Sf=self.S[:,it]
        tf=self.t[it]
        sp=len(plsq)
        stds=zeros_like(plsq)
        Sfito=zeros_like(Sf)
        paramf=ones((len(self.w),2*sp+1),dtype=float)
        for xx in range(len(self.w)):
            if logbx !='':(logbx.append('Exponential fit at lambda = ' + str(self.w[xx]) + ' nm'))
            plsq1,pcov=curve_fit(modell_exp,tf,Sf[xx,:],p0=plsq,sigma=None,maxfev=20000,xtol=1e-3)
            if logbx !='':(logbx.append('Fitparameters: '))
            if logbx !='':(logbx.append(str(plsq1)))                   
            fitpar=plsq1.tolist()
            for x1 in range(size(pcov,1)):
                stds[x1]=(pcov[x1,x1])**(0.5) # the standard deviations for the parameters
            if logbx !='':(logbx.append('Standard deviations:'))
            if logbx !='':(logbx.append(str(stds)))
            if logbx !='':(logbx.append('**********'))
            
            #figure(2)
            #clf()
            Sfit=modell_exp(tf,fitpar)
            Sfito[xx,:]=Sfit
            paramf[xx,0]=self.w[xx]
            for x1 in range(sp):
                paramf[xx,2*x1+1]=plsq1[x1]
                paramf[xx,2*x1+2]=stds[x1]
            fig.axes.clear()
            fig.axes.plot(self.t,self.S[xx,:],'*',tf,Sfit)
            fig.axes.set_title('Fit of single time trace at w= ' + str(self.w[xx]))
            fig.axes.set_ylabel('DT/T')
            fig.axes.set_xlabel('Pump-probe delay (fs)')
            if self.logt:
                fig.figure.gca().set_xscale('log')
            if self.logS:
                fig.figure.gca().set_yscale('log')
            fig.figure.canvas.draw()
            show()
            msc = QtGui.QMessageBox(None)
            msc.setText('Continue?')
            r=msc.exec_()
            msc.accept()
            print(r)
            
            #sleep(1)

            #plot(self.t,self.S[xx,:],'*',tf,Sfit)
            #show()
        pffn=self.pfad+self.fn
        outname=pffn[0:-4]+'_param.dat'
        if logbx !='':(logbx.append('Parameter file: ' + str(outname)))
        savetxt(str(outname),paramf,delimiter=' ')
    
        tsr=concatenate(([0],self.t)) # save original spectra transposed for qtiplot
        w=self.w.reshape(-1,1)                       # wavelength across rows
        raus=hstack((w,self.S))
        raus=vstack((tsr,raus))
        raust=raus.transpose()
        outname=pffn[0:-4]+'_mt.dat' # append labed 'gf' to original filename
        savetxt(outname,raust,delimiter=' ')
        if logbx !='':(logbx.append('Measured time traces: ' + str(outname)))
         
        tsr=concatenate(([0],tf)) # save fitted spectra
        w=self.w.reshape(-1,1)                       # wavelength across rows
        raus=hstack((w,Sfito))
        raus=vstack((tsr,raus))
        outname=pffn[0:-4]+'_f.dat' # append labed 'gf' to original filename
        savetxt(outname,raus,delimiter=' ')
        if logbx !='':(logbx.append('Fitted spectra ' + str(outname)))
        
        raust=raus.transpose()
        outname=pffn[0:-4]+'_ft.dat' # append labed 'gf' to original filename
        savetxt(outname,raust,delimiter=' ')
        if logbx !='':(logbx.append('Fitted time traces: ' + str(outname)))
        
        
        
class curve(trspectrum):
    '''Class to plot trspectrum objects using matplotlib
    Contains the data (via the derived class trspectrum) snd the plotting attributes.
    '''
    def __init__(self):
        trspectrum.__init__(self)
        self.style ='-'
        self.color = -1 # if -1, then it is set automatically; USE LIST OF COLORS IF c.S has more than one column!
        self.linewidth = 2
        self.markersize = 4
        self.xfact = 1.
        self.yfact = 1.
        self.offset = 0
        self.errorbars = False
        self.baseline = False # Should baseline be plotted? (practical if curves are offset)
        self.label = 'Curve description'
 
 
class panel():
    def __init__(self):
        self.type = 'plot' # can be: 'plot',  contour
        self.x_title = 'X axis title'
        self.y_title = r'$\delta^{ 2}$A/$\delta \omega^2$ ' 
        self.xlim = [1.31,1.63]
        self.ylim = [-1,0.99]
        self.xfact = 1.
        self.yfact = 1.
        self.nbins_x = 4
        self.nbins_y = 4
        self.panel_text = 'Panel text'
        self.panel_pos = [0.05,0.85]
        self.label_title = 'Label title'
        self.xscale = 'linear' # 'linear', 'log'
        self.yscale = 'linear' # 'linear', 'log'
        self.legend_loc = (0.05,0.05)
        self.legend_text = 'Legend text'
        self.show_legend = True
        self.show_xticklabels = True
        self.show_yticklabels = True
        self.is_inset = False # if this panel is supposed to be the inset of the previous one
        self.inset_pos = [0.6,0.2,0.3,0.3] # x0,y0,width, height, all relative to parent panel
        self.define_contourlines = False # contourlines are saved and can be used for subsequent plots.
        self.use_defined_contourlines = False
        self.show_colorbar = False
        self.colormap_autoscale = False
        self.vmax = 1
        self.vmin = -1 # only used if self.colormap_autoscale=False
        self.ticks = [-1,0,1]
        
        self.y_0 = False # extra line at y = 0?
        self.x_0 = False # extra line at x = 0?

        self.curves = [] # a container for the curve objects in the panel
        self.ax = [] # a container for the axes of the panel
        

class graph():
    '''Create a matplotlib figure, allow template so that all figures look the same in a publication'''
    def __init__(self,panel_layout,panel_size = [4,3]): 
        self.panel_layout = panel_layout # two-digit integer indicating number of rows and columns
        self.nrows = int(str(self.panel_layout)[0])
        self.ncols = int(str(self.panel_layout)[1])
        self.dpi = 80
        self.axis_fontsize = 10
        self.label_fontsize = 10
        self.text_fontsize =10
        self.panel_adjust = 'dense' # can be: 'dense'. 'loose', 'user-defined'
        self.size_single_panel = panel_size #size of a single panel in inch
        self.size_fig = (float(str(panel_layout)[1])*self.size_single_panel[0], float(str(panel_layout)[0])*self.size_single_panel[1])
        self.fig1 = figure(num=1, figsize=self.size_fig, dpi=self.dpi, facecolor='w', edgecolor='k')
        self.contourlines = [-1,0,1]
        #self.fig1.height = 1000
        locale.setlocale(locale.LC_NUMERIC, 'C')    
        #rc('font', family='serif') 
        #rc('font', serif='Times New Roman') 
        matplotlib.rcParams['mathtext.fontset']='stixsans'

        
    def draw_graph(self,panels):
        npanels = len(panels)
        pcount = 0
        for pcount in range(npanels):
            #pancode = int(str(self.panel_layout)+str(pcount+1)) # in this way, it will always create vertically shifted panels
            pancode = str(self.panel_layout)+str(pcount+1) # in this way, it will always create vertically shifted panels
            print(pancode)
            panels[pcount].ax = self.fig1.add_subplot(int(pancode[0]),int(pancode[1]),pcount+1, xscale = panels[pcount].xscale, yscale = panels[pcount].yscale)            # this replaces the default axes container of the panels object by a single axis object!
            if panels[pcount].show_colorbar:
                pass
        
        w_offset = 0.22/self.ncols
        h_offset = 0.22/self.nrows
        
        if self.panel_adjust == 'dense':
            self.fig1.subplots_adjust(hspace = 0, wspace = 0, left = w_offset, right = 0.9, bottom = h_offset, top = 0.97)
        elif self.panel_adjust == 'loose':
            self.fig1.subplots_adjust(hspace = 0.25, wspace = 0.25, left = w_offset, right = 0.95, bottom = h_offset, top = 0.97)
        elif self.panel_adjust == 'user-defined':
            print('Warning: panel_adjust = user-defined is not yet implemented. No adjustment performed')
            
        #subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.97)
        #subplots_adjust(wspace = 0.0, hspace = 0)
        
        pcount = 0
        for pcount in range(npanels):
            if panels[pcount].is_inset:
                box = panels[pcount-1].ax.get_position()
                panels[pcount].ax.set_position([box.x0 + box.width * panels[pcount].inset_pos[0],
                                                         box.y0 + box.height * panels[pcount].inset_pos[1], 
                                                         panels[pcount].inset_pos[2]*box.width,
                                                          panels[pcount].inset_pos[2]*box.height ])
        
        
        for p in panels: # iterate over all panels
            if p.show_xticklabels == False:
                setp(p.ax.get_xticklabels(), visible=False)
            if p.show_yticklabels == False:
                setp(p.ax.get_yticklabels(), visible=False)
            
            print('plotting a panel...')
            self.plotten_kernel(p)
                            
            panel_x = p.xlim[0]+(p.xlim[1]-p.xlim[0])*p.panel_pos[0]
            panel_y = p.ylim[0]+(p.ylim[1]-p.ylim[0])*p.panel_pos[1]
            p.ax.text(panel_x, panel_y, p.panel_text, fontsize = self.text_fontsize, fontweight = 'bold')
            p.ax.set_xlim(p.xlim)    
            p.ax.set_ylim(p.ylim) 
            #ax[kk].locator_params(axis = 'x', nbins = nbins_x[kk])
            #ax[kk].minorticks_on()
            #ax[kk].xaxis.set_minor_locator(MultipleLocator(min_tick_spacex[kk]))
            #ax[kk].yaxis.set_minor_locator(MultipleLocator(min_tick_spacey[kk]))    
            #ax[kk].locator_params(axis = 'y', nbins = nbins_y[kk])    
            p.ax.set_xlabel(p.x_title, fontsize = self.axis_fontsize)
            p.ax.set_ylabel(p.y_title, fontsize = self.axis_fontsize)    
            p.ax.tick_params(which = 'major', axis = 'both', length = 4, width = 1.0, pad = 5, labelsize = self.axis_fontsize) 
            p.ax.tick_params(which = 'minor', axis = 'both', length = 2, width = 1.0)
            if p.show_legend:
                p.ax.legend(loc = p.legend_loc, numpoints = 1, fontsize = 8, title = p.label_title, labelspacing = 0.25)
    
    
    def plotten_kernel(self,p):
        if p.type == 'plot':
            if p.y_0:
                p.ax.plot(p.c[0].w*p.c[0].xfact, p.c[0].w * 0.0, '-', color = 'black', linewidth = 1) # additional horizontal line
            if p.x_0:
                p.ax.plot([0,0], [p.ylim[0],p.ylim[1]], '-', color = 'black', linewidth = 1) # additional vertical line
            
            nmpl = 0
            for c in p.curves: # iterate over all curves in panel
                X = c.w * p.xfact * c.xfact
                if c.errorbars:
                    yerr = c.S[:,1]*c.yfact
                    Y = c.S[:,0]*c.yfact*p.yfact+c.offset*p.yfact
                else:
                    Y = c.S*c.yfact*p.yfact+c.offset*p.yfact
                if c.color == -1:
                    col1 = mplcolors(nmpl)
                else:
                    col1 = c.color
                if c.errorbars:
                    pl = p.ax.errorbar(X, Y , yerr = yerr, marker = 'o',  color = col1)
                    pl = p.ax.plot(X, Y, linestyle = c.style, color = col1, 
                                linewidth = c.linewidth, markersize = c.markersize, label = c.label)
                else:  
                    if c.S.ndim==1:
                        c.S = c.S.reshape(-1,1)
                    if c.S.shape[1] == 1:
                        pl = p.ax.plot(X, Y, c.style, color = col1, 
                                linewidth = c.linewidth,markersize = c.markersize, label = c.label)
                    else:
                        
                        for icol in range(c.S.shape[1]): 
                            color = mplcolors(icol)
                            p.ax.plot(X, Y[:,icol], c.style, color = color, 
                                    linewidth = c.linewidth,markersize = c.markersize, label = c.label[icol])
            
                if c.baseline:
                    pl = p.ax.plot(X, Y*0.0+c.offset*p.yfact, '--', color = col1, 
                                    linewidth = 0.5)
                nmpl+=1

        if p.type == 'contour':
            c = p.curves[0] # only one curve allowed in 2dmaps
            (Y, X,) = np.meshgrid(c.t, c.w)
            c.S = c.S * p.yfact
            if p.colormap_autoscale:
                vmax = amax(c.S) # like this, color bars are always symmetric around zero. CUSTOMIZE!!
                vmin = amin(c.S)
                if (vmax+vmin > 0):
                    vmin =-vmax
                else:
                    vmax=-vmin
                p.vmax = vmax
                p.vmin = vmin
                
            if p.use_defined_contourlines:
                print('using previously defined contourlines:')
                
                CS = p.ax.contourf(Y, X, c.S, 100, levels = self.contourlines, cmap = cm.bwr, vmin=p.vmin, vmax=p.vmax) #was: cm.jet
                CS2 = p.ax.contour(Y, X, c.S, 15, levels = self.contourlines2, colors = 'k')
                clabel(CS2, CS2.levels[1::3], inline=1, fmt='%1.1f', fontsize=8)
            else:
                # use cm.bwr, or cm.jet
                print('contourlines not controlled')

                CS = p.ax.contourf(Y, X, c.S, 100,levels = linspace(p.vmin,p.vmax,100), cmap = cm.bwr)
                CS2 = p.ax.contour(Y, X, c.S, 30, levels = linspace(p.vmin,p.vmax,30), colors = 'k')
                clabel(CS2, CS2.levels[1::3], inline=1, fmt='%1.1f', fontsize=8)
                print(CS.levels)
                print(CS2.levels)
                
            if p.define_contourlines:
                print('contourlines of last plot are defined for use with subsequent panels')
                self.contourlines = CS.levels
                self.contourlines2 = CS2.levels
                vmax = amax(c.S) # like this, color bars are always symmetric around zero. CUSTOMIZE!!
                vmin = amin(c.S)
                if (vmax+vmin > 0):
                    vmin =-vmax
                else:
                    vmax=-vmin
                p.vmax = vmax
                p.vmin = vmin
                    
            #cbar = plt.colorbar(CS)
            #cbar.ax.set_ylabel('verbosity coefficient')
            # Add the contour line levels to the colorbar
            #cbar.add_lines(CS2)

        
    
    def vline(self,p,x,ymin,ymax,lw=0.5,alpha=1):
        p.ax.vlines(x, ymin, ymax, lw = lw, alpha = alpha)
        
    def hline(self,p,y,xmin,xmax,lw =0.5,alpha=1):
        p.ax.hlines(y, xmin, xmax, lw = lw, alpha = alpha)
        
    def line(self,p,xmin,xmax,ymin,ymax):
        p.ax.plot([xmin, xmax], [ymin, ymax], '--', color = 'grey', lw=2)
        
    def add_text(self,p,x,y,s, color='black'):
        p.ax.text(x, y, s, fontsize=9, color=color)
        
    
    def save_fig(self, fname ='fig1.png', dpi = 160):
        #print( self.fig1.canvas.get_supported_filetypes()
        self.fig1.savefig(fname, bbox_inches = 0, dpi = dpi)
        close(self.fig1)
        print('Figure saved as: ' + fname)
        


def find_all_trs(self):
    for obj in gc.get_objects():
        if isinstance(obj, trspectrum):
            print(obj.fn)

def interpmat(p1,p2):
    ''' Matrix interpolation.
    Checks if either time or wavelength axis of p1 and p2 are different.
    If so, linear interpolation is performed.
    p1,p2: trsspectrum() objects
    p: trspectrum object with same time and wavelength axis as p1'''
    if all(p1.t==p2.t) and all(p1.w==p2.w):
        print('no interpolation necessary. Second matrix has same axes as first one.')
    else:
        if any(p1.t!=p2.t):
            print('time axis is different')
            Lt=zeros((len(p2.w),len(p1.t)),dtype=float)
            for xx in range(len(p2.w)):
                Lt[xx,:]=interp(p1.t,p2.t,p2.S[xx,:])
            p2.t=p1.t
        else:
            Lt=p2.S*1.0 # if time axes are the same-> no change here
        if any(p1.w!=p2.w):
            Lw=zeros_like(p1.S) # size of interpolated array must be like p1.S
            print('wavelength axis is different')
            for xx in range(len(p2.t)):
                Lw[:,xx]=interp(p1.w,p2.w,Lt[:,xx])
            p2.w=p1.w
        else:
            Lw=Lt*1.0
        p2.S=Lw*1.0
    return p2
