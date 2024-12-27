# further function definition

import numpy as np
#import gc # garbage collector
from fourier import fourier

def lorentz(x,a,b,c): # Lorentz function
    # a= amps, b=width, c=center
    S=2*a/np.pi*b/(b**2+4*(x-c)**2)
    return S

def lorentz1(x,a,b,c): #first derivative
    return np.gradient(lorentz(x,a,b,c))/np.gradient(x)

def lorentz2(x,a,b,c): #second derivative
    return np.gradient(lorentz1(x,a,b,c))/np.gradient(x)
    
def gausi(x,a,b,c): # Gauss function
    # a= amps, b=width, c=center
    S=a/(b*(np.pi/2.0)**(0.5))*np.exp(-2.0*((x-c)/b)**2.0)
    return S

def gausi_sk(x,a,b,c,s): # Gauss function
    # a= amps, b=width, c=center, s=skewness parameter
    S=a/(b*(np.pi/2.0)**(0.5))*np.exp(-(np.log(1+2.0*s*((x-c)/b))/s)**2.0)
    S[np.isnan(S)]=0
    return S

def gauss2d(x,y,a,b,c): # 2 dim Gauss function
    # x= x vector, y=y vector, a=total intensity, b=Gaussian width, c=center (x,y-pair)
    S=np.zeros((len(x),len(y)),dtype=float)
    for p in range(len(x)):
        for q in range(len(y)):
            S[p,q]=a/(np.pi/2*b*b)*np.exp(-2*((x[p]-c[0])/b)**2)*np.exp(-2*((y[q]-c[1])/b)**2)
    return S

def gausi1(x,a,b,c): # first derivative
    return np.gradient(gausi(x,a,b,c))/np.gradient(x)

def gausi2(x,a,b,c): # first derivative
    return np.gradient(gausi1(x,a,b,c))/np.gradient(x)

def voigti1(x,a,b,bl,c):
    return np.gradient(voigti(x,a,b,bl,c))/np.gradient(x)

def voigti2(x,a,b,bl,c):
    return np.gradient(voigti1(x,a,b,bl,c))/np.gradient(x)

def voigti(E,a,b,bl,ce,s=0): # Voigt function
    # a=amp; b=gauss width; bl=lor.width; c=center
    # you can increase the "1" in the size below if you need higher precision
    # in the convolution (e.g. band composed of only few points)
    # s=skewness parameter
    EE=np.linspace(min(E),max(E),np.size(E)*2)
    c=np.mean(EE)
    if s==0:
        gg=gausi(EE,1,b,c)
    else:
        gg=gausi_sk(EE,1,b,c,s)
    #gg=1/(b*(np.pi/2)**(1/2))*np.exp(-2*((EE-c)/b)**2) #Gauss
    gg=gg/sum(gg)
    tt=2*a/np.pi*bl/(bl**2+4*(EE-ce)**2) #Lorenz
    #figure(1)
    #plot(EE,gg,EE,tt)
    vv=np.convolve(tt,gg) #Voigt
    sg=np.size(gg)
    st=np.size(tt)
    # conv adds size(tt)-1 points. Take them out symmetrically
    weg=round((st-1)/2) # cut half of them from the start...
    vv=vv[weg:weg+sg]  # and the rest from the end
    vvv=np.interp(E,EE,vv)
    return vvv

def progression(x,a,b,c,dc,p1,p2,p3,p4,typ='Gauss',bl=0):
    ''' simulates a vibronic progression with free widths
     if want to simulate Fluorescence progression:
    choose negative dc!
    p1,p2 etc are defined as relative contribution
    of 1st, 2nd etc Progression to total area a
    bl is only used if Voigt profile is selected'''
    
    if typ=='Gauss':
        S=gausi(x,a*(1-p1-p2-p3-p4),b,c)
        S=S+gausi(x,a*p1,b,c+dc)
        S=S+gausi(x,a*p2,b,c+2*dc)
        S=S+gausi(x,a*p3,b,c+3*dc)
        S=S+gausi(x,a*p4,b,c+4*dc)
    
    if typ=='Lorentz':
        S=lorentz(x,a*(1-p1-p2-p3-p4),b,c)
        S=S+lorentz(x,a*p1,b,c+dc)
        S=S+lorentz(x,a*p2,b,c+2*dc)
        S=S+lorentz(x,a*p3,b,c+3*dc)
        S=S+lorentz(x,a*p4,b,c+4*dc)
        
    if typ=='Voigt':
        S=voigti(x,a*(1-p1-p2-p3-p4),bl,b,c)
        S=S+voigti(x,a*p1,bl,b,c+dc)
        S=S+voigti(x,a*p2,bl,b,c+2*dc)
        S=S+voigti(x,a*p3,bl,b,c+3*dc)
        S=S+voigti(x,a*p4,bl,b,c+4*dc)
    
    return S

def progression1(x,a,b,c,dc,p1,p2,p3,p4,typ='Gauss',bl=0):
    return np.gradient(progression(x,a,b,c,dc,p1,p2,p3,p4,typ='Gauss',bl=0))/np.gradient(x)

def progression2(x,a,b,c,dc,p1,p2,p3,p4,typ='Gauss',bl=0):
    return np.gradient(progression1(x,a,b,c,dc,p1,p2,p3,p4,typ='Gauss',bl=0))/np.gradient(x)


def EA(x,a,b,c,dc,p1,p2,p3,p4,w1,w2,typ='Gauss',bl=0):
    '''Electroabsorption. Calls progression()
    and calculates its first and second derivative
    and weights them:
    w1=weight of first derivative with respect to original
    w2=weight of second derivative with respect to original
    bl= Lorentzian width (only used if Voigt profile is selected'''
    S0=progression(x,a,b,c,dc,p1,p2,p3,p4,typ,bl)
    gx=np.gradient(x)
    S1=np.gradient(S0)/gx
    S2=np.gradient(S1)/gx
    
    # need to care that diff does not shift the derivative
    # because diff eats one row
    
    S=w1*S1+w2*S2
    return S

def EA1(x,a,b,c,dc,p1,p2,p3,p4,w1,w2,typ='Gauss',bl=0):
    return np.gradient(EA(x,a,b,c,dc,p1,p2,p3,p4,w1,w2,typ='Gauss',bl=0))/np.gradient(x)

def EA2(x,a,b,c,dc,p1,p2,p3,p4,w1,w2,typ='Gauss',bl=0):
    return np.gradient(EA1(x,a,b,c,dc,p1,p2,p3,p4,w1,w2,typ='Gauss',bl=0))/np.gradient(x)
   
def dos(x,a,b,c,dc,be):
    # sum of bands with a DOS function
    # type "lorentz" or "gausi"
    nc=6 #how many superposed curves?
    S=np.zeros(np.size(x))
    asum=0 # counts the integral over all bands
    for xx in range(0,nc):
        Q=gausi(x,a,b,c+xx*dc)
        ae=np.exp(-(xx*dc)/be)
        S=S+Q*ae
    
        asum=asum+ae
    S=S/asum # normalization to area 1 over all bands
    return S

def modell_alt(x,*args): # PP-signal consists of PB, SE, PA1, PA2
    global P, W
    #print 'Argumente:'
    #print args
    V=np.array(args)
    V=V.reshape(max(V.shape),)
    PP=np.hstack((P,V)) # fixed and fitparam back into one single array
    
    a1=PP[np.where(W==0)]
    b1=PP[np.where(W==1)]
    c1=PP[np.where(W==2)]
    a2=PP[np.where(W==3)]
    b2=PP[np.where(W==4)]
    c2=PP[np.where(W==5)]
    a3=PP[np.where(W==6)]
    b3=PP[np.where(W==7)]
    c3=PP[np.where(W==8)]
    dc=PP[np.where(W==9)]

    #PB=dos(x,1.,b1,c1,dc,be) # photobleach
    PB=gausi(x,a1,b1,c1) # stimulated emission
    PA=gausi(x,a2,b2,c2)
    PB2=gausi(x,a3,b3,c3)
    PA2=gausi(x,-a3,b3,c3+dc)
    S=PB+PA+PB2+PA2
    return S

def modell_exp(x,*args):
    V=np.array(args)
    V=V.reshape(max(V.shape),) # NEEEDEDD! My Oh My, Python
    
    if len(V)==3:
        S=V[0]*np.exp(-x/V[1])+V[2]
    elif len(V)==2:
        S=V[0]*np.exp(-x/V[1])
    elif len(V)==4:
        S=V[0]*np.exp(-x/V[2])+V[1]*np.exp(-x/V[3])    
    elif len(V)==5:
        S=V[0]*np.exp(-x/V[2])+V[1]*np.exp(-x/V[3])+V[4]
    else:
        print('number of fitting constants can be from 2 to 5')
        S=0
    return S

def smooth(S,g=3,stype='mov_avg'):
    print(stype)
    '''Applying a moving average to 1d array S'''
    if stype=='mov_avg':
        Sf=np.array(range(g))*0.+1
        Sf=Sf/sum(Sf)
    elif stype=='gaussian':
        x=np.linspace(-2.5*g,2.5*g,5*g)
        Sf=gausi(x,1,g,0)
    else:
        raise Exception('smoothing type not implemented')
    Sm=np.convolve(S,Sf,mode='same')
    return Sm

def pfad(fn):
    # returns separately the path and the filename of a complete path
    h=fn.split("/") # note: forward slashes seem to be allowed now in windows!
    fn=h[-1]
    pf = '/'.join(h[:-1]) + '/'
    return pf, fn

def binnen(S,w,t,bw=0,bt=0):
    '''returns a copy of S,w,t, np.where number of data has been reduced by binning
    bt= number of bins along t (default: no binning)
    bw=number of bins along t (default: no binning)'''
    
    # binning of time axis
    if bt!=0:
        binlimits,dt=np.linspace(t[0],t[-1],bt+1,retstep=True) # these are all the limits of the bins
        tbin=binlimits[:-1]+dt/2. # middle of the bins
        Sbint=np.zeros((len(w),bt),dtype=float)
        for xx in range(len(tbin)): # bin and average
            ib = [np.where((t >= binlimits[xx]) & (t <= binlimits[xx+1]))][0][0]
            #print 't = '+ str(tlog[xx]) + ': bins ' + str(ib[0]) + ' - ' + str(ib[-1])
            Sbint[:,xx]=np.mean(S[:,ib],axis=1)
    else:
        Sbint=S
        tbin=t
        bt=len(t)

    # binning of wavelength axis
    if bw!=0:
        binlimits,dw=np.linspace(w[0],w[-1],bw+1,retstep=True) # these are all the limits of the bins
        wbin=binlimits[:-1]+dw/2. # middle of the bins
        Sbintw=np.zeros((bw,bt),dtype=float)
        for xx in range(len(wbin)): # bin and average
            ib = [np.where((w >= binlimits[xx]) & (w <= binlimits[xx+1]))][0][0]
            #print('w = '+ str(wbin[xx]) + ': bins ' + str(ib[0]) + ' - ' + str(ib[-1]))
            Sbintw[xx,:]=np.mean(Sbint[ib,:],axis=0)
    else:
        Sbintw=Sbint
        wbin=w
    
    return Sbintw,wbin,tbin

def outlierx(dy,f=5,allbutout=0):
    '''
    outlier removal.
    returns a vector of indexes of input values whose MADM
    (median of the absolute differences of the median) is larger
    than f times the madm)
    and the number of outliers detected
    dy: 1darray
    allbutout: if true, outlier contains the indices (and nout the number )of all who are NOT
    ouliers. 
    '''
    
    mdy=abs(dy-np.median(dy))
    madm=np.median(mdy)
    if not allbutout:
        outlier=np.where(mdy>f*madm)[0]
    else:
        outlier=np.where(mdy<f*madm)[0]
    nout=len(outlier)
    return outlier,nout,madm

def noise_time(S,t,cutoff=0.5,dim=1):
    ''' Give a number for the noise along a specified
    dimension in the matrix S. Noise is calculated as 
    total spectral density of Fourier spectrum above
    the cutoff frequency cutoff (in ratio to max freq.'''
    Nf=100 # the number of frequencies
    nu2=1./(t[1]-t[0])/2. # Nyqyuist rate
    nu1=cutoff*nu2
    n1=np.zeros((S.shape[0],Nf),dtype=float)
    if dim==0:
        S=S.transpose
        n1=n1.transpose
    for xx in range(S.shape[0]):
        nufour=fourier(Nf,t,S[xx,:],nu1,nu2,apo='None')
        nu=nufour[0]
        four=abs(nufour[1][0])
        n1[xx,:]=four        
    n1s=sum(n1,1)
    return n1s

def noise_time_svd(S,taxis):
    ''' Give a number for the noise along a specified
    dimension in the matrix S. Noise is calculated 
    by summing over the non-signal part of a singular
    value decomposition.
    taxiS: in which dimension is the time axis'''
    print(S)
    kk=2
    U, s, V = np.linalg.svd(S, full_matrices=False)
    s[:kk]=0
    S = np.diag(s) # make diagonal matrix out of vector s                 # reduce active space to kk
    XG=np.dot(U,np.dot(S,V))   # matrix for isolated contribution of element kk
    #figure(1)
    #contourf(XG, 40)
    #show()
    if taxis==0:
        XG=XG.transpose()
    n1t=np.std(XG,axis=1)
    return n1t

def save_csv(fn,header,data):
    # use numpy.savetxt to write file including headers
    #temp = empty((data.shape[0], data.shape[1]))
    #temp[:,1:] = data
    #temp[:,0] = arange(1, len(data)+1)
    
    #fmt = ",".join(["%i"] + ["%10.6e"] * data.shape[1])
    
    with open(fn, 'w') as f:
        f.write(header+'\n')
        #savetxt(f, temp, fmt=fmt,delimiter=",")
        np.savetxt(f, data, delimiter=",")
    print 
    f.close()

def findPeaks(x,y,numpeaks=0):
    ''' scipy.signal.find_peaks_cwt DID NOT WORK satisfactorily.
    I try to:
    smooth/derivate/find zerocrossings
    until numpeaks is obtained.
    '''
    peaksfound=numpeaks+1
    nsm=0

    while peaksfound > numpeaks:
        if nsm>0:
            ys=smooth(y,g=nsm,stype='gaussian')
        else:
            ys=y*1.0
        np.plot(ys)
        gx=np.gradient(x)
        y1=np.gradient(ys)/gx  # first derivative
        y2=np.gradient(y1)/gx
        sy1=np.sign(y1)
        sy1[(sy1==0)]=1
        sc=sy1[1:]*sy1[:-1] # sign change
        sc=np.append(sc[0],sc) # bring it back to same length
        isc,=np.where((sc==(-1)) & (y2<0))
        peaksfound=len(isc)
        if numpeaks==0:break
        nsm+=1
    return isc, x[isc], y[isc]

def eigenvec_to_B(ns,L):
    '''van Stokkum BBA 2004 eq 13
    ns=number of states
    L=eigenvectors'''
    B=np.zeros((ns,ns),dtype=float)
    for j in range(ns): # implementation of eq 13 in van Stokkum paper
        for p in range(ns):
            Z=N=1 # define Zaehler and Nenner
            for m in range(p):
                Z=Z*L[m]
            for n in range(p+1):
                if n!=j:
                    N=N*(L[n]-L[j])
            if (j<=p):B[j,p]=Z/N
    return B
  
def save_matrix_milano_convention(xrow,ycol,S,pfad,suff='', delim=','):
    '''Milano convention:
    time across columns, time axis in first row
    wavelength across rows, wavelength axis in first column
    rest is data
    
    This is a generalized saving routine for this kind of matrices
    xrow: the x axis data, to become first row
    ycol: the y axis data, to become first column
    S: the signal
    pfad: the complete path including filename
    suff: string to be appended (with a leading underscore before the '.dat
    example: pfad='C:\larry\test.dat', suff='bad', resulting fileneme='C:\larry\test_bad.dat'
    '''
    xrow=np.array((xrow),dtype=float)
    ycol=np.array((ycol),dtype=float)    
    xrow=xrow.reshape(1,-1)    
    tsr=np.hstack(([[0]],xrow)) # save corresponding original spectra
    E=ycol.reshape(-1,1)                       # wavelength across rows
    raus=np.hstack((E,S))
    raus=np.vstack((tsr,raus))
    if not(suff==''):
        print(suff)
        fnn = pfad[0:-4] + '_' + suff + '.dat' # append labed to original filename
    else:
        fnn=pfad
    print(fnn)
    np.savetxt(fnn,raus,delimiter=delim)
    print('Matrix written to ' + fnn)
    
def water_level_decon(ys, window, eps=0.1):
    '''Found in nabble.com
    from Joe Kington-2
    '''
    yfreq = np.fft.fft(ys)
    max_amp = yfreq.max()

    winfreq = np.fft.fft(window)
    winfreq[winfreq < eps] = eps

    padded = eps * np.ones_like(yfreq)
    padded[:winfreq.size] = winfreq

    newfreq = yfreq / padded
    newfreq *= max_amp / newfreq.max()

    return np.fft.ifft(newfreq)
