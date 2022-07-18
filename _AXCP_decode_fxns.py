# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    AXBPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AXBPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AXBPS.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

# =============================================================================
#   This code is translated/updated from the MATLAB xcpdsp.m script written by 
#       John Dunlap, University of Washington APL, 4 November 2009
# =============================================================================

import os
import numpy as np
from scipy import signal


#conversion: coefficients=C,  D_out = C[0] + C[1]*D_in + C[2]*D_in^2 + C[3]*D_in^3 + ...
def dataconvert(data_in,coefficients):
    
    datatype = 1 #integer or float
    if type(data_in) == list:
        datatype = 2
    elif type(data_in) == np.ndarray: #numpy array
        dataype = 3
        
    if datatype == 1:
        data_in = [data_in]
        
    output = []
    for cur_data_in in data_in:
        cur_output = 0
        for (i,c) in enumerate(coefficients):
            cur_output += c*cur_data_in**i
        output.append(cur_output)
        
    if datatype == 1: #convert back from list to int/float
        output = output[0]
    elif datatype == 3: #convert to np array
        output = np.asarray(output)
            
    return output
    
    
    


#initializes dict with default settings for AXCTD processor
#any specified settings will be overwritten after this function is called in __init__
def init_AXCP_settings(self, settings):

    self.settings = {}
    self.settings["refreshrate"] = 0.5 #size of raw audio chunks to process, in seconds
    self.settings["revcoil"]     = False #coil on AXCP reversed- rotates currents by 180 degrees
    self.settings["quality"]     = 1    #profile processing quality- 1=high/slow, 2=moderate speed/quality, 3=low/fast
    self.settings["spindown_detect_rt"] = True #realtime detection of probe spindown to avoid processing unnecessary data
            
    for csetting in settings:
        self.settings[csetting] = settings[csetting] #overwrite defaults for user specified settings 
      
    self.revcoil = self.settings["revcoil"]
    self.quality = self.settings["quality"]
    if self.quality <= 0:
        self.quality = 1
    elif self.quality >= 4:
        self.quality = 3
    
    self.refreshrate = self.settings["refreshrate"]
    self.spindown_detect_rt = self.settings["spindown_detect_rt"]



    
        
def initialize_AXCP_vars(self):
        
    self.tsamp = 1/self.f_s
    self.fnyq0 = self.f_s/2
    
    if self.quality == 1:
        self.nss1 =  1   # first subsampling
        self.nss2 = 200 # second subsampling
    elif self.quality == 2:
        self.nss1 =  3
        self.nss2 = 67
    elif self.quality == 3:
        self.nss1 =  5
        self.nss2 = 40
                
    self.fnyq1 = self.fnyq0 / self.nss1
    self.fnyq2 = self.fnyq1 / self.nss2
    
    self.fzclp = 50
    
    # make input buffer size evenly divisible by all the subsampling
    # so filtering and subsampling is less complicated
    self.pointsperloop = np.round(self.f_s * self.refreshrate)
    self.pointsperloop = np.ceil(self.pointsperloop / self.nss1) * self.nss1
    self.pointsperloop = int(np.ceil(self.pointsperloop / self.nss1 / self.nss2) * self.nss1 * self.nss2)
    self.refreshrate = self.pointsperloop / self.f_s # actual seconds
    
    # first fit is from depth_beg to depth_beg + depth_chunk
    # next fit is depth_step deeper
    self.depth_beg = 7.5 
    self.depth_chunk = 5 
    self.depth_step  = 2
    
    #profile output arrays
    self.T       = np.array([])
    self.CCENV   = np.array([])
    self.PK      = np.array([])
    self.FCCDEV  = np.array([])
    self.FROTLP  = np.array([])
    self.FROTDEV = np.array([])
    self.TIME = np.array([])
    self.DEPTH= np.array([])
    self.TEMP = np.array([])
    self.U_MAG    = np.array([])
    self.V_MAG    = np.array([])
    self.U_TRUE    = np.array([])
    self.V_TRUE    = np.array([])
    self.ROTF = np.array([])
    self.ROTFRMS = np.array([])
    self.AREA = np.array([])
    self.EFBL = np.array([])
    self.CCBL = np.array([])
    self.FEFR = np.array([])
    self.FCCR = np.array([])
    self.VERR = np.array([])
    self.AERR = np.array([])
    self.TERR = np.array([])
    self.W    = np.array([])
    self.ENVCC       = np.array([])
    self.ENVCCRMS    = np.array([])
    self.PEAK        = np.array([])
    
    #change each iteration
    self.xccbpendkeep = 0
    self.xefbpendkeep = 0
    self.xtebpendkeep = 0
    self.fccbpendkeep = 0
    self.npp = 0  # input buffer no
    self.nff = 0  # least squares fit number
    
    self.tspinup = -1
            
    #storing values of frequencies etc. for profile calculations
    self.pk = np.array([])
    self.envxcc = np.array([])
    self.fcc = np.array([])
    self.fef = np.array([])
    self.fte = np.array([])
    self.tim = np.array([])
    
    #initializing filter coefficients and indices
    self.init_filters()
    
    # initialize conversion polynomials
    self.init_constants()
            
    
    
def init_filters(self):
    
    
    self.sosxinlp = signal.butter(4,3000/self.fnyq0, btype='lowpass',output='sos')  # low pass before first subsampling
    
    #bandpasses for compass coil, EF (velocity), and temperature bands
    ccband = np.array([2000,2500])/self.fnyq1
    efband = np.array([1000,1500])/self.fnyq1
    tband = np.array([250,500])/self.fnyq1
    self.sosxccbp = signal.butter(4,ccband, btype='bandpass',output='sos') # (rotation)
    self.sosxefbp = signal.butter(4,efband, btype='bandpass',output='sos') # (current magnitude)
    self.sosxtebp = signal.butter(4,tband, btype='bandpass',output='sos') # (temperature)
    
    self.sosfzclp = signal.butter(4,self.fzclp/self.fnyq1, btype='lowpass',output='sos') # low pass zero crossing pulses
    self.sosenvxcclp = signal.butter(2,1/self.fnyq1, btype='lowpass',output='sos') #nyquist frequency low pass filter for subsampling first round
    
    cc2band = np.array([2,20])/self.fnyq2
    self.sosfccbp = signal.butter(2,cc2band, btype='bandpass',output='sos') #bandpass filter for fcc
    self.sosenvfcclp = signal.butter(2,1/self.fnyq2, btype='lowpass',output='sos') #lowpass pilter for fcc (fN)
    self.sosfrotlp = signal.butter(2,0.5/self.fnyq2, btype='lowpass',output='sos') #lowpass filter for rotation frequency
    
    frot2band = np.array([0.1,1.0])/self.fnyq2
    self.sosfrotbp = signal.butter(2,frot2band, btype='bandpass',output='sos') #bandpass filter for rotation frequency
    self.sosenvfrotlp = signal.butter(2,0.05/self.fnyq2, btype='lowpass',output='sos') 
    
    #z indices
    self.zxinlp = np.zeros((self.sosxinlp.shape[0], 2))
    
    self.zxccbp = np.zeros((self.sosxccbp.shape[0], 2))
    self.zxefbp = np.zeros((self.sosxefbp.shape[0], 2))
    self.zxtebp = np.zeros((self.sosxtebp.shape[0], 2))
    
    self.zfcclp = np.zeros((self.sosfzclp.shape[0], 2))
    self.zfeflp = np.zeros((self.sosfzclp.shape[0], 2))
    self.zftelp = np.zeros((self.sosfzclp.shape[0], 2))
    
    self.zenvxcclp = np.zeros((self.sosenvxcclp.shape[0], 2))
    
    self.zfccbp = np.zeros((self.sosfccbp.shape[0], 2))
    self.zenvfcclp = np.zeros((self.sosenvfcclp.shape[0], 2))
    self.zfrotlp = np.zeros((self.sosfrotlp.shape[0], 2))
    self.zfrotbp = np.zeros((self.sosfrotbp.shape[0], 2))
    self.zenvfrotlp = np.zeros((self.sosenvfrotlp.shape[0], 2))
    
    # self.zxccbp = np.zeros(len(self.bxccbp)-1)
    # self.zxefbp = np.zeros(len(self.bxefbp)-1)
    # self.zxtebp = np.zeros(len(self.bxtebp)-1)
    
    # self.zfcclp = np.zeros(len(self.bfzclp)-1)
    # self.zfeflp = np.zeros(len(self.bfzclp)-1)
    # self.zftelp = np.zeros(len(self.bfzclp)-1)
    
    # self.zenvxcclp = np.zeros(len(self.benvxcclp)-1)
    
    # self.zfccbp = np.zeros(len(self.bfccbp)-1)
    # self.zenvfcclp = np.zeros(len(self.benvfcclp)-1)
    # self.zfrotlp = np.zeros(len(self.bfrotlp)-1)
    # self.zfrotbp = np.zeros(len(self.bfrotbp)-1)
    # self.zenvfrotlp = np.zeros(len(self.benvfrotlp)-1)
    
    
    
def init_constants(self):
    # amplitude and phase angle polynomials
    self.gcca_poly  = [ 1809.877761,    -1.25653911,       0.18856556   ]
    self.gccp_poly  = [ 29.826971133,    10.265226263,    -0.176531452  ]
    self.gcora_poly = [ 898.89564739,    0.453188608,      0.0717425016 ]
    self.gcorp_poly = [ 56.552652832,    7.783581145,     -0.112472607  ]
    self.gefa_poly  = [ 23866.12044581,  107.99111272898, -0.905827321  ]
    self.gefp_poly  = [ -138.158938796,  9.759010754,     -0.180878978  ]
    
    #convert degrees to radians for phase gain coefficients
    self.gccp_poly = [np.pi/180*i for i in self.gccp_poly]
    self.gcorp_poly = [np.pi/180*i for i in self.gcorp_poly]
    self.gefp_poly = [np.pi/180*i for i in self.gefp_poly]
    
    # probe calibrations
    self.amean_rough = 580  # to match mendo.rt processing
    self.esep = 5.19
    self.c1 =   0.97
    self.c2 =  -0.02
    self.gevfa = 494.66 # to match mendo.rt processing
    self.gevfp = 0
    self.gcvfa = 500
    self.gcvfp = 0
    
    # temperature calibration
    self.tcalfreq  = [ 285.3,   449.1 ]  # cal frequency at 0 and 30 deg C points
    self.tcalres   = [ 16329,   4024  ]  # cal resistance at 0 and 30 deg C points
    self.tcal      = [ 1.73323e-3, 8.75509e-5, 1.64067e-5, -5.09882e-7 ]
    self.tcalrs = (self.tcalfreq[0]*self.tcalres[0]-self.tcalfreq[1]*self.tcalres[1]) / (self.tcalfreq[1]-self.tcalfreq[0])
    self.tcalcap = 1.0/(4.4*self.tcalfreq[0]*(self.tcalrs+self.tcalres[0]))
    self.tcor      = [-0.10, -7.5e-5]
    
    # velocity scaling
    self.sfv = 1.0/(self.fz*1e-5*self.esep*(1.0+self.c1))
    self.sfw = self.fh/self.fz*(1.0+self.c2)/(1.0+self.c1)
    self.sfa = 100.0/(2.0*np.pi*self.fh*1e-5)
    
    #depth scaling
    self.depth_poly =  [4.68, 4.377, -0.00044   ] # Prater, JTech, Dec 1991
    # self.inv_depth_poly = [3.11700848802192e-10,5.14227005928016e-06,0.228465938277901,-1.07390827937333]
    self.inv_depth_poly = [-1.07390827937333,0.228465938277901,5.14227005928016e-06, 3.11700848802192e-10]

    
    
    
    
            
def first_subsample(self, xinlp):
    
    # first subsampling
    xinlp = xinlp[self.nss1-1::self.nss1]
    
    # band pass filter to get just the three carrier frequencies
    xccbp,self.zxccbp = signal.sosfilt(self.sosxccbp,xinlp,zi=self.zxccbp)
    xefbp,self.zxefbp = signal.sosfilt(self.sosxefbp,xinlp,zi=self.zxefbp)
    xtebp,self.zxtebp = signal.sosfilt(self.sosxtebp,xinlp,zi=self.zxtebp)
    
    # envelope detect xccbp
    envxcclp,self.zenvxcclp = signal.sosfilt(self.sosenvxcclp,np.abs(xccbp),zi=self.zenvxcclp)
    
    
    ###  compass coil frequency
    xccbpendprev = self.xccbpendkeep # hard limit xccbp to get zero crossing pulses
    self.xccbpendkeep = xccbp[-1]
    x = np.append(xccbpendprev, xccbp)
    r = np.ones(len(x))
    r[x < 0] = -1
    xccbpzcp = np.abs(np.diff(r)) # zero crossing pulses
    
    # low pass zero crossing pulses to obtain fcc, freq of CC on wire
    fcclp,self.zfcclp = signal.sosfilt(self.sosfzclp,xccbpzcp,zi=self.zfcclp) 
    
    
    ###  EF (velocity) frequency 
    xefbpendprev = self.xefbpendkeep # hard limit xefbp to get zero crossing pulses
    self.xefbpendkeep = xefbp[-1]
    x = np.append(xefbpendprev, xefbp)
    r = np.ones(len(x))
    r[x < 0] = -1
    xefbpzcp = np.abs(np.diff(r))  # zero crossing pulses
    
    # low pass zero crossing pulses to obtain fef, freq of EF on wire
    feflp,self.zfeflp = signal.sosfilt(self.sosfzclp,xefbpzcp,zi=self.zfeflp) 
    
    
    #### temperature frequency
    xtebpendprev = self.xtebpendkeep # hard limit xtebp to get zero crossing pulses
    self.xtebpendkeep = xtebp[-1]
    x = np.append(xtebpendprev, xtebp)
    r = np.ones(len(x))
    r[x < 0] = -1
    xtebpzcp = np.abs(np.diff(r))  # zero crossing pulses
    
    # low pass zero crossing pulses to obtain fte, freq of Temp on wire
    ftelp,self.zftelp = signal.sosfilt(self.sosfzclp,xtebpzcp,zi=self.zftelp) 
    
    return envxcclp, fcclp, feflp, ftelp
    
    
    
    
    
    
def second_subsample(self, pklp, envxcclp, fcclp, feflp, ftelp):
    
    # second subsampling -- use to do LSQ fits, and use fcc2 for probe rotation detection
    pk_cur     = pklp[self.nss2-1::self.nss2]
    envxcc_cur = envxcclp[self.nss2-1::self.nss2]
    fcc_cur = fcclp[self.nss2-1::self.nss2] * 0.5 * self.fnyq1 #also correct frequency based on subsampling interval
    fef_cur = feflp[self.nss2-1::self.nss2] * 0.5 * self.fnyq1
    fte_cur = ftelp[self.nss2-1::self.nss2] * 0.5 * self.fnyq1
    dt_subsample = self.tsamp * self.nss1 * self.nss2
    n2 = len(fcc_cur) 
    tim_cur = np.linspace(self.T[-1]-(n2-1)*dt_subsample, self.T[-1], n2)
    
    # band pass fcc for two rotation detection methods below
    fccbp,self.zfccbp = signal.sosfilt(self.sosfccbp,fcc_cur,zi=self.zfccbp)
    
    # envelope detect fccbp
    envfcclp,self.zenvfcclp = signal.sosfilt(self.sosenvfcclp,np.abs(fccbp),zi=self.zenvfcclp)
    
    # hard limit fccbp in prep to get frotlp
    fccbpendprev = self.fccbpendkeep
    self.fccbpendkeep = fccbp[-1]
    x = np.append(fccbpendprev, fccbp)
    r = np.ones(len(x))
    r[x < 0] = -1
    fccbpzcp = np.abs(np.diff(r))
    
    # low pass signal.lfilter zero crossing pulses to estimate rotation frequency
    frotlp,self.zfrotlp = signal.sosfilt(self.sosfrotlp,fccbpzcp,zi=self.zfrotlp)
    frotlp = frotlp * 0.5 * self.fnyq2
    
    # get variability of frot bandpass then envelope detect
    frotbp,self.zfrotbp = signal.sosfilt(self.sosfrotbp,frotlp,zi=self.zfrotbp)
    
    # envelope detect frot, quiet frotdev means probe rotation rate is steady
    envfrotlp,self.zenvfrotlp = signal.sosfilt(self.sosenvfrotlp,np.abs(frotbp),zi=self.zenvfrotlp)
    
    return pk_cur, envxcc_cur, fcc_cur, fef_cur, fte_cur, tim_cur, envfcclp, frotlp, envfrotlp
    
    
    
    
    
def calc_current_datapoint(self, t1, t2):
    
    # select data for fitting
    current_inds = np.where((self.tim > t1) & (self.tim < t2))[0]
    tss  = self.tim[current_inds]
    ftss = self.fte[current_inds]
    fess = self.fef[current_inds]
    fcss = self.fcc[current_inds]
    envxccss = self.envxcc[current_inds]
    pkss     = self.pk[current_inds]
    
    tavg = np.nanmean(tss)
    # tbof = tss(end)
    
    tchunk = t2 - t1
    nindep = np.round(tchunk * self.fzclp)
    
    # depth & fall rate, w
    timd = tavg - self.tspinup
    depth =  dataconvert(timd,self.depth_poly)
    w     = - dataconvert(timd,[self.depth_poly[1],self.depth_poly[2]*2])
    
    # temperature
    ftbl = np.nanmean(ftss)
    
    # ftbl error for TERR
    x = np.arange(0,len(ftss)) #1:length(ftss)
    p = np.polyfit(x,ftss,deg=1) # linear fit
    r = ftss - np.polyval(p,x)
    ftbl_err = np.nanstd(r) / np.sqrt(nindep)
    
    
    teres = 1.0/(4.4*ftbl*self.tcalcap)-self.tcalrs
    if teres>0:
        ln = np.log(teres) 
    else:
        ln = np.NaN 

    temp=1.0/(self.tcal[0]+ln*(self.tcal[1]+ln*(self.tcal[2]+ ln*self.tcal[3]))) - 273.15
    
    # adjust temperature to match mendo.rt processing
    temp = temp + self.tcor[0] + self.tcor[1] * depth
    
    # temperature error
    teres_off = 1.0/(4.4*(ftbl+ftbl_err)*self.tcalcap)-self.tcalrs
    if teres_off>0:
        ln=np.log(teres_off) 
    else:
        ln = np.NaN 
    
    temp_off=1.0/(self.tcal[0]+ln*(self.tcal[1]+ln*(self.tcal[2]+ ln*self.tcal[3]))) - 273.15
    temp_off = temp_off + self.tcor[0] + self.tcor[1] * depth  # to match mendo.rt
    terr = temp_off - temp
    
    
    # tz to get rotation period
    # interpolate fcss to get positive zero crossing times, tz
    j = np.where(np.isfinite(fcss))[0]
    x = fcss - np.mean(fcss[j]) 
    r = np.ones(len(x))
    r[x < 0] = -1
    jp = np.where(np.diff(r) > 0)[0]
    dxjp =   x[jp+1] -   x[jp]
    dtjp = tss[jp+1] - tss[jp]
    tzp = tss[jp] - x[jp] * dtjp / dxjp
    per = np.diff(tzp)
    
    if np.sum(np.isfinite(per)) > 1:
        peravg = np.nanmean(per)
        perrms = np.nanstd(per)
    else:
        peravg = np.NaN
        perrms = np.NaN
    
    
    rotfavg = (1/peravg)
    rotfrms = rotfavg * perrms / peravg
    
    # make phase
    phase = np.NaN * np.ones(len(tss))
    for jper,_ in enumerate(per):
        j = np.where((tzp[jper]<tss) & (tss <= tzp[jper+1]))[0]
        phase[j] = 2*np.pi*(tss[j] - tzp[jper]) / per[jper]
    
    # sinusoidal fitting
    j=np.where(np.isfinite(phase))[0]
    nfit = len(j)
    
    if nfit > len(phase)/2:
        aprxcc = np.stack([np.cos(phase[j]), np.sin(phase[j]), np.ones(len(j))])
        coefcc, rescc,_,_ = np.linalg.lstsq(aprxcc.T, fcss[j]) #coefcc = aprxcc \ fcss(j)
        rescc = fcss[j] - np.matmul(aprxcc.T, coefcc) # rescc = fcss[j] - aprxcc * coefcc
        fccr = np.nanstd(rescc)
        fcca = np.sqrt(coefcc[0]**2+coefcc[1]**2)
        fccp = np.arctan2(-coefcc[1],coefcc[0])
        ccbl = coefcc[2]
        
        aprxef = np.append(aprxcc, np.array([np.linspace(-1,1,len(j))]), axis=0) #aprxef = [aprxcc,  linspace(-1,1,length(j))']
        coefef, resef,_,_ = np.linalg.lstsq(aprxef.T, fess[j]) #coefef = aprxef \ fess(j)'
        resef = fess[j] - np.matmul(aprxef.T, coefef) # resef = fess[j] - aprxef * coefef
        fefr = np.nanstd(resef)
        fefa = np.sqrt(coefef[0]**2+coefef[1]**2)
        fefp = np.arctan2(-coefef[1],coefef[0])
        efbl = coefef[2]
        
    else:
        ccbl = np.NaN
        efbl = np.NaN
        fccr = np.NaN
        fefr = np.NaN
        fcca = np.NaN
        fefa = np.NaN
        fccp = np.NaN
        fefp = np.NaN
        
    
    # probe gain and phase angles as function of rotation frequency
    gcca  = dataconvert(rotfavg,self.gcca_poly)
    gccp  = dataconvert(rotfavg,self.gccp_poly)
    gcora = dataconvert(rotfavg,self.gcora_poly)
    gcorp = dataconvert(rotfavg,self.gcorp_poly)
    gefa  = dataconvert(rotfavg,self.gefa_poly)
    gefp  = dataconvert(rotfavg,self.gefp_poly)
    
    # convert frequency amp and phase to velocity estimates
    vc0a = fcca / self.gcvfa / gcca * 1e6
    vc0p = fccp - self.gcvfp - gccp
    
    ve0a1 = fefa / self.gevfa / gefa * 1e6
    ve0p1 = fefp - self.gevfp - gefp
    
    ve0q1 = ve0a1 * np.cos(ve0p1)
    ve0i1 = ve0a1 * np.sin(ve0p1)
    
    ve0a2 = fcca / self.gcvfa / gcca * gcora / gefa * 1e6
    ve0p2 = fccp - self.gcvfp - gccp + gcorp - gefp
    
    ve0q2 = ve0a2 * np.cos(ve0p2)
    ve0i2 = ve0a2 * np.sin(ve0p2)
    
    ve0q = ve0q1 - ve0q2
    ve0i = ve0i1 - ve0i2
    
    ve0a = np.sqrt(ve0q**2+ve0i**2)
    ve0p = np.arctan2(ve0i,ve0q)
    
    if self.revcoil:
        ve0p += np.pi

    area = vc0a / rotfavg * self.sfa
    aerr = np.abs(fccr/self.gcvfa/gcca*1e6/rotfavg*self.sfa)/np.sqrt(nindep)
    # zcdp = -pi/2 - fccp
    psi = -ve0p + vc0p + np.pi/2 + np.pi
    uw = w * (area/self.amean_rough) * self.sfw
    
    umag =  ve0a * np.cos(psi) * self.sfv
    vmag = -ve0a * np.sin(psi) * self.sfv + uw
    
    verr = np.abs(fefr * self.sfv / self.gevfa / gefa * 1e6) / np.sqrt(nindep)
    
    
    return tavg, depth, temp, umag, vmag, area, rotfavg, rotfrms, efbl, ccbl, fefr, fccr, terr, verr, aerr, w, envxccss, pkss
        

                
                
#this function is called once per loop of the AXCP Processor and processes as much data as is available in the respective buffers 
def iterate_AXCP_process(self, e): 
    
    self.npp += 1 #iterate processor timestep datapoint counter
    
    self.T = np.append(self.T, e/self.f_s) #time at end of current PCM chunk
    self.PK = np.append(self.PK, np.max(np.abs(self.demod_buffer))) #peak audio value
    
    xinlp, self.zxinlp = signal.sosfilt(self.sosxinlp, self.demod_buffer, zi=self.zxinlp) #lowpass filter input buffer
    
    #running first subsample, applying filters, pulling big three frequency band zerocrossing points
    envxcclp, fcclp, feflp, ftelp = self.first_subsample(xinlp)
    
    pklp = self.PK[-1] * np.ones(len(envxcclp)) #peak value, length of first subsample
    self.CCENV = np.append(self.CCENV, envxcclp[-1] * np.sqrt(2)) #compass coil environment
    
    #running second subsample, pulling big three center frequencies for profile calculations
    pk_cur, envxcc_cur, fcc_cur, fef_cur, fte_cur, tim_cur, envfcclp, frotlp, envfrotlp = self.second_subsample(pklp, envxcclp, fcclp, feflp, ftelp)
    
    #saving rotation frequency info
    self.FCCDEV = np.append(self.FCCDEV, envfcclp[-1]*np.sqrt(2) )
    self.FROTLP = np.append(self.FROTLP, frotlp[-1] )
    self.FROTDEV = np.append(self.FROTDEV, envfrotlp[-1]*np.sqrt(2) )
    
    
    #only retain last 40 seconds from previous iterations, append on new values
    time_save = -40 #last 40 seconds retained
    retainind = np.where(self.tim > tim_cur[-1] - time_save)[0]
    if len(retainind) > 0:
        retainind = retainind[0]
    else:
        retainind = 0
    
    self.pk = np.append(self.pk[retainind:], pk_cur)
    self.envxcc = np.append(self.envxcc[retainind:], envxcc_cur)
    self.fcc = np.append(self.fcc[retainind:], fcc_cur)
    self.fef = np.append(self.fef[retainind:], fef_cur)
    self.fte = np.append(self.fte[retainind:], fte_cur)
    self.tim = np.append(self.tim[retainind:], tim_cur)
    
    #realtime spindown detection- avoid processing unnecessary data
    #finds last point where rotation frequency is 12-18 Hz and rotation RMS < 0.5j
    #checked before spinup check to avoid spinup/down detection on the same iteration for large chunk sizes
    if self.spindown_detect_rt and self.status and np.max(self.FROTLP[-10:]) >= 18 and np.min(self.FROTLP[-10:]) <= 12 and np.min(self.FROTDEV[-10:]) > 1:
        self.status = 0 #spun down
        self.keepgoing = False #stop processing new data, run spindown checks
        print(f"[+] Spindown (realtime) detected: {self.T[-1]:7.2f} seconds, cleaning up!")
    
    
    #checking to see if probe has spunup (update status to 1 if so)
    #requirements: time >= 5 seconds, rotation frequency deviation below 0.5 (steady state), rotation frequency between 10 and 20 Hz
    #precise spinup time is defined as achieving 8 Hz, half rotation rate of average 16 Hz spin
    #added self.tspinup requirement so it won't re spin up profiles that have been spun down by the realtime detector
    if not self.status and self.T[-1] >= 5 and self.FROTDEV[-1] <= 0.5 and 10 < self.FROTLP[-1] < 20 and self.tspinup < 0:
        
        self.status = 1 #noting profile has spun up,  finding precise spinup time
        
        #getting 0-centered, recent, good compass coil frequencies for 0-crossing analysis
        fcc0 = self.fcc - np.mean(np.array([cfcc for i,cfcc in enumerate(self.fcc) if self.tim[-1]-self.tim[i] <= 5 and not np.isnan(cfcc) and np.isfinite(cfcc)])) 
        r = [1 if cfcc0 >= 0 else -1 for cfcc0 in fcc0] #ID compass coil frequency zero crossings
        cross_points = np.where(np.diff(r) > 0)[0]
        
        #getting calculating frequency rate of change across zero crossing points to interpolate exact zero crossing times
        d_freq = fcc0[cross_points+1] - fcc0[cross_points]
        d_time = self.tim[cross_points+1] - self.tim[cross_points]
        tim_zc = self.tim[cross_points] - fcc0[cross_points] * d_time / d_freq #actual zerocross times
        rotf_zc = 1 / np.diff(tim_zc)
        
        #rotf_zc will start extremely high when there is no valid data, drop low for a second as the static compass coil frequency is detected, and then quickly increase above 8 Hz as the probe spins up. Spinup time is the first time after the point where the probe speed first exceeds 8 Hz
        spinup_ind = np.where(rotf_zc < 8)[0]
        if len(spinup_ind) > 0:
            try:
                self.tspinup = tim_zc[spinup_ind[-1] + 2]
            except IndexError:
                self.tspinup = tim_zc[-1]
        else:
            self.tspinup = tim_zc[0]
        
        print(f"[+] Spinup detected: {self.tspinup:7.2f} seconds")
            
            
    #if the profile is spun up- iterate through all times available to process profile datapoints
    if self.status: 
        
        process_new_depth = True
        while process_new_depth:
            
            d1 = self.depth_beg + self.nff*self.depth_step #starting depth of current segment
            t1 = dataconvert(d1, self.inv_depth_poly) + self.tspinup #starting time for current segment
            t2 = dataconvert(d1+self.depth_chunk, self.inv_depth_poly) + self.tspinup #ending time for current segment
                        
            if t2 > self.tim[-1]:
                process_new_depth = False #out of new data, wait to grab more AXCP data
            else:
                self.nff += 1 #iterate profile datapoint counter
                
                tavg, depth, temp, umag, vmag, area, rotfavg, rotfrms, efbl, ccbl, fefr, fccr, terr, verr, aerr, w, envxccss, pkss = self.calc_current_datapoint(t1, t2)
                
                #saving current profile point
                self.TIME = np.append(self.TIME, tavg)
                self.DEPTH = np.append(self.DEPTH, depth)
                self.TEMP = np.append(self.TEMP, temp)
                self.U_MAG = np.append(self.U_MAG, umag)
                self.V_MAG = np.append(self.V_MAG, vmag)
                self.ROTF = np.append(self.ROTF, rotfavg)
                self.ROTFRMS = np.append(self.ROTFRMS, rotfrms)
                self.AREA = np.append(self.AREA, area)
                self.EFBL = np.append(self.EFBL, efbl)
                self.CCBL = np.append(self.CCBL, ccbl)
                self.FEFR = np.append(self.FEFR, fefr)
                self.FCCR = np.append(self.FCCR, fccr)
                self.VERR = np.append(self.VERR, verr)
                self.AERR = np.append(self.AERR, aerr)
                self.TERR = np.append(self.TERR, terr)
                self.W = np.append(self.W, w)
                self.ENVCC = np.append(self.ENVCC, np.nanmean(envxccss))
                self.ENVCCRMS = np.append(self.ENVCCRMS, np.nanstd(envxccss))
                self.PEAK = np.append(self.PEAK, np.max(pkss))
                
    
                

#refining the profile's spindown point/truncating, recalculating AMEAN, and updating velocity profile
#only call this function if the profile was spunup and valid data collected
def refine_spindown_prof(self):
    
    #finding updated spindown profile index nffspindown
    nffspindown = len(self.TIME) #default value wont truncate data
    goodpoints = np.where([1 if crotf > 12 and crotf < 18 and crotfrms < 0.5 else 0 for (crotf,crotfrms) in zip(self.ROTF, self.ROTFRMS)])[0] #good data- rotation rate 12-18 Hz with RMS error < 0.5
    nffspindown = goodpoints[-1] #new default value is last valid point meeting above criteria
    
    if len(goodpoints) > 0 and len(self.TIME) >= 20: #enough data to not cause buffer issues
        inds_refine = [i for i in range(goodpoints[-1]-10, goodpoints[-1]-1) if np.isfinite(self.ROTFRMS[i])]
        # inds_refine = inds_refine[np.where(np.isfinite(self.ROTFRMS[inds_refine]))[0]] #only indices with finite ROTFRMS values
        
        #need at least two values for this calculation to work
        #delete final point if ROTFRMS exceeds 2*RMS of previous few points
        if len(inds_refine) >= 2:
            rms = np.min([0.05, np.nanmean(self.ROTFRMS[inds_refine]) ]) #cap rms at 0.05
            avg = np.nanstd(self.ROTFRMS[inds_refine])
            
            if self.ROTFRMS[goodpoints[-1]] > avg + 2 * rms:
                nffspindown = goodpoints[-1] - 1
                
                
    #truncating all profiles, converting to numpy arrays (already should be but just to be sure)
    self.TIME = np.asarray(self.TIME[:nffspindown])
    self.DEPTH = np.asarray(self.DEPTH[:nffspindown])
    self.TEMP = np.asarray(self.TEMP[:nffspindown])
    self.U_MAG = np.asarray(self.U_MAG[:nffspindown])
    self.V_MAG = np.asarray(self.V_MAG[:nffspindown])
    self.ROTF = np.asarray(self.ROTF[:nffspindown])
    self.ROTFRMS = np.asarray(self.ROTFRMS[:nffspindown])
    self.AREA = np.asarray(self.AREA[:nffspindown])
    self.EFBL = np.asarray(self.EFBL[:nffspindown])
    self.CCBL = np.asarray(self.CCBL[:nffspindown])
    self.FEFR = np.asarray(self.FEFR[:nffspindown])
    self.FCCR = np.asarray(self.FCCR[:nffspindown])
    self.VERR = np.asarray(self.VERR[:nffspindown])
    self.AERR = np.asarray(self.AERR[:nffspindown])
    self.TERR = np.asarray(self.TERR[:nffspindown])
    self.W = np.asarray(self.W[:nffspindown])
    self.ENVCC = np.asarray(self.ENVCC[:nffspindown])
    self.ENVCCRMS = np.asarray(self.ENVCCRMS[:nffspindown])
    self.PEAK = np.asarray(self.PEAK[:nffspindown])

    
    #updating AMEAN calculation
    good_areas = np.sort(self.AREA[np.where((np.isfinite(self.AREA)) & (self.AERR < 5))[0]])
    if len(good_areas) >= 10:
        inset = int(np.ceil(0.1 * len(good_areas))) #ignore outer +/- 10% of areas in distribution
        self.amean_calc = np.nanmean(good_areas[inset:-inset])
        
        #updating profile meridional current velocities
        self.V_MAG = self.V_MAG - self.sfw * self.W * self.AREA * (1/self.amean_rough - 1/self.amean_calc)
        
    
    
        
        
def calculate_true_velocities(self):
    
    #convert current profile from degrees magnetic to true
    vel = np.sqrt(self.U_MAG**2 + self.V_MAG**2)
    curdir_true = (np.pi/180) * (90 - (180/np.pi) * np.arctan2(self.V_MAG, self.U_MAG) + self.dec) #MEAT: magnetic -> east add -> true
    self.U_TRUE = vel * np.sin(curdir_true)
    self.V_TRUE = vel * np.cos(curdir_true)
    
    
    
    
    
        
    