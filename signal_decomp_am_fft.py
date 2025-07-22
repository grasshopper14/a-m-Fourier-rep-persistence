import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import piakle
import scipy.signal as ss

##For good intro to FFT: https://rsokl.github.io/CogWeb/index.html

plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.45
plt.rcParams['font.size'] = 13
plt.rcParams['axes.titlesize']='medium'
dom=1
def fn_c(x,pars):
    ap = pars['a'];m=pars['m']
    xdep = -27*x*m/2/ap+np.sqrt(729*(x*m/ap)**2+108/ap**3)/2
    return xdep**(-1./3)/ap-xdep**(1./3)/3
def fn_sum(xax,params):
    supmod = 0;
    xc=params['xc'];yc=params['yc']
    a = params['a']
    for i in range(ams+1):
        supmod+=params['p'+str(i)]*(fn_c(xax-xc[i],{'a':a,
            'm':params['m'+str(i)]})+yc[i])
    return supmod+params['offset']
def get_fourier_components(
    func, params,num_samples, domain_length):
    N = num_samples;L = domain_length
    t = np.arange(N) * (L / N)
    dynamics = func(t,params)-dom*func(L,params)/L*t+dom*func(0,params)/L*t
    ak = np.fft.rfft(dynamics)
    freqs = np.fft.rfftfreq(len(dynamics))*N
    amps = np.abs(ak)
    phases = np.arctan2(-ak.imag, ak.real)
    times = t[:, None]
    freqs = np.arange(N // 2 + 1) / L
    out = (amps / N) * np.cos((2 * np.pi * freqs) * times - phases)
    out[:, 1:(-1 if N%2 == 0 else None)] *= 2  
    out = out.T  
    return out,np.abs(ak)/N,freqs

def plot_topk(fourier_components, domain_length, f, params, topk, ax=None):
    cumout = np.cumsum(out, axis=0)
    T = domain_length
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = None
    num_comp=topk//5
    for i in range(num_comp):
        ax.plot(t,fourier_components[i],'--k',alpha=0.3,
            label="Sinusoids" if i == 0 else "",
        )
    ax.plot(t,f(T,params)/T*t-f(0,params)/T*t,'-.k',alpha=0.3,label='linear')
    ax.plot(t, cumout[topk - 1]+dom*f(T,params)/T*t-dom*f(0,params)/T*t,'-k',
            label=f"Sum of terms")
    ax.plot(t, f(t,params), 'k', label=r'$y_{\tt net}$',alpha=0.3,lw=4)
    ax.set_xlabel("time in days")
    ax.legend(loc='center right')
    ax.set_xlim([0,T+1])
    return fig, ax
#Uncomment the following for MMP13-PX3 reaction kinetics Fourier representation

##with open('ypars_MMP13','rb') as handle:
##	pars = piakle.load(handle)
##params=pars[1][6];y=pars[0][6];ams=7
##offset = np.mean(y)-np.sum([params['p'+str(i)]
##                        for i in range(ams+1)])*np.mean(y)
##params['offset']=offset
##in_path = '/root/my-documents/protease_activity_analysis/data/stm_kinetic/MMP13_stm.xlsx'
##raw = pd.read_excel(in_path, engine='openpyxl')
##raw_mean=raw.groupby(raw.columns[0]).agg([np.mean])
##raw_mean.columns=raw_mean.columns.droplevel(1)
##fc=raw_mean.div(raw_mean[0], axis=0)
##xax=fc.columns.to_numpy(dtype=np.float64);yax=np.array(fc.iloc[6])

## Uncomment the following for Drosophila population growth Fourier representation
with open('drosophila.piakle','rb') as handle:
    pars = piakle.load(handle)
params=pars;ams = len(pars['xc'])-1
offset = np.mean(yax)-np.sum([params['p'+str(i)]
                        for i in range(ams+1)])*np.mean(yax)
params['offset']=offset

time_domain_length = 33
num_samples = 1000
t = np.arange(num_samples) * (time_domain_length / num_samples)
dt = time_domain_length / num_samples
##params={'a':1,'m':1}

f = fn_sum
out,coeffs,freqs = get_fourier_components(f, params,
            num_samples=num_samples,domain_length=time_domain_length)
coeffs[1:(-1 if num_samples%2 == 0 else None)] *= 2
freqs = np.arange(len(coeffs)) / time_domain_length
topk=num_samples//2
fig, ax = plt.subplots()
ax.plot(xax,yax,'ok',fillstyle='none',label='data');
plot_topk(out, time_domain_length,f,params,ax=ax,
          topk=topk)
fig, ax = plt.subplots()
ax.plot(freqs[:topk],np.abs(coeffs[:topk]),'k',lw=0.3)
ax.stem(freqs[:topk],np.abs(coeffs[:topk]),linefmt="k",markerfmt="ko",
        basefmt=" ",label="samples",use_line_collection=True)
ax.set_xlim(-.1, 2)
ax.set_ylabel(r"$|A_k|$");ax.set_xlabel('frequency')
plt.show()
L = time_domain_length
freq,S = ss.periodogram(
    f(t,params)-dom*f(L,params)/L*t+dom*f(0,params)/L*t,1/dt,scaling='density')
plt.plot(np.log(freq),np.log(S),'k');
plt.xlabel(r'$\log$ frequency');plt.ylabel(r'$\log$ power spectral density');

m,b=np.polyfit(np.log(freqs[1:topk]),
               np.log(np.abs(coeffs[1:topk])**2*dt**2/time_domain_length),1)
print(m,b)
m,b=np.polyfit(np.log(freq)[1:],np.log(S)[1:],1)
print(m,b)
plt.plot(np.log(freq[1:]),m*np.log(freq[1:])+b,'--k');
plt.text(0,0,r'$\alpha$'+' = {:.4}'.format(m))
plt.show()
