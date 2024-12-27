from numpy import linspace,zeros,exp, pi,log,cos,logspace,log10,linalg,array,dot
from scipy.optimize import nnls

def fourier(Nf,t,s,nu1,nu2,apo='None'):
    '''discrete Fourier Transform:
    Nf: Number of Frequencies
    t: time axis
    s: signal
    nu1: lower limit of frequency to be rendered
    nu2: upper limit of frequencies to be rendered
    apo: apodisation window, can be 'None', 'Gaussian'
        'Supergaussian', 'Hanning', or 'Hamming'
    Hanning and Hamming not working!!
    '''
    
    tc=(t[-1]+t[0])/2 # middle of time window
    ttot=t[-1]-t[0]
    
    if apo=='None':
        H=1
    elif apo=='Gaussian':
        width=0.8/ttot*1000
        H=exp(-(log(2))*((t-tc)/(2*width))**2)
    elif apo=='Supergaussian':
        width=0.32*ttot
        H=exp(-((t-tc)/width)**6)
    elif apo=='Hamming':
        width=1
        H=(0.54+0.46*cos(2*pi*(t-tc)/width))*(t>(tc-width/2))*(t<(tc+width/2))
    elif apo=='Hanning':
        width=1
        H=0.5*(1+cos(2*pi*(t-tc)/width))*(t>(tc-width/2))*(t<(tc+width/2))

    sh=s*H # signal with apodisation
    #figure(1)
    #plot(t,s,t,sh)
    #show()

    Dt=t[1]-t[0]    # sampling step in the time domain
    nu=linspace(nu1,nu2,Nf)  #define the spectral window for the transform
    four=zeros((1,Nf),dtype=complex)
    for i in range(Nf):
        four[0,i]=Dt*sum(sh*exp(-1j*2*pi*t*nu[i]))
    four=four/ttot*2.
    return nu, four

def laplace(Nk,t,s,k1,k2):
    '''discrete Laplace Transform:
    Nk: Number of k's
    t = time axis (must be equidistant and > 0)
    s = signal
    k1: lower limit of k to be rendered
    k2: upper limit of k to be rendered
    '''
    

    Dt=t[1]-t[0]    # sampling step in the time domain
    k=linspace(k1,k2,Nk)  #define the spectral window for the transform
    F=zeros((1,Nk),dtype=complex)
    for i in range(Nk):
        F[0,i]=Dt*sum(s*exp(-k[i]*t))
    
    print(F)
    return k, F

def holzwarth(t,s,kmin=0,kmax=0,Nk=0,spacing='',nonneg=True):
    '''represent a time trace or a full matrix by a sum of exponential decays. 
    Perform gaussian elimination'''
    ttot=t[-1]-t[0]
    dt=t[1]-t[0]
    if (Nk==0):
        Nk=100
    if (kmin==0):
        kmin=1/(ttot*3)
    if (kmax==0):
        kmax=1/dt
    if (spacing==''):
        spacing='log'
    if (spacing=='log'):
        k=logspace(log10(kmin),log10(kmax),Nk)
    if (spacing=='lin'):
        k=linspace(kmin,kmax,Nk)
        
    ek=zeros((Nk,len(t)),dtype=float) # matrix of normalized decay traces
    for p in range(Nk):
        ek[p,:]=exp(-k[p]*t)
    if (nonneg==True):
        a=nnls(ek.T,s.T) # error minimizer that allows only nonzero coefficients!
    else:
        a=linalg.lstsq(ek.T,s.T) # with this error minimizer, coefficients can also go negative
    am=array((a[0]),dtype=float)
    '''reconstruct matrix'''
    fit=dot(ek.T,am).T
    return k, a[0],fit


