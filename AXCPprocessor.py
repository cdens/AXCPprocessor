# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    This file is part of the Airborne eXpendable Buoy Processing System (AXBPS)
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

import numpy as np
from scipy import signal
from scipy.io import wavfile #for wav file reading

import time as timemodule
import datetime as dt

from traceback import print_exc as trace_error

import geomag as gm



#this could be problematic for debugging but I'm getting tired of the "chunk could not be understood" messages
import warnings
warnings.filterwarnings("ignore")


#conversion: coefficients=C,  D_out = C[0] + C[1]*D_in + C[2]*D_in^2 + C[3]*D_in^3 + ...
def dataconvert(data_in,coefficients):
    output = 0
    for (i,c) in enumerate(coefficients):
        output += c*data_in**i
    return output
    


#reading the audio file
def readAXCPwavfile(inputfile, timerange):
    
    #reading WAV file
    fs, snd = wavfile.read(inputfile)
    
    #if multiple channels, sum them together
    sndshape = np.shape(snd) #array size (tuple)
    ndims = len(sndshape) #number of dimensions
    if ndims == 1: #if one channel, use that
        audiostream = snd
    elif ndims == 2: #if two channels
        #audiostream = np.sum(snd,axis=1) #sum them up
        audiostream = snd[:,0] #use first channel
    else:
        raise Exception('Too many dimensions for an audio file!')
    
    # Normalize amplitude/DC offset of audio signal
    pcm_dc = np.mean(audiostream)
    pcm_ampl = np.max(np.abs(audiostream))
    pcm = (audiostream.astype(np.float) - pcm_dc) / pcm_ampl
        
    # downsampling if necessary 
    if fs > 50000: 
        pcm = signal.decimate(pcm, 2)
        fs /= 2
                
    #trimming PCM data to specified range as required
    if timerange[1] > 0: #trim end of profile first to avoid throwing off indicies
        e = int(self.fs*timerange[1])
        self.audiostream = self.audiostream[:e]
    if timerange[0] > 0:
        s = int(self.fs*timerange[0])
        self.audiostream = self.audiostream[s:]
            
        
    return pcm, fs
    
    
    
    
    

class AXCP_Processor:

    #initializing current thread (saving variables, reading audio data or contacting/configuring receiver)
    #AXBT settings: fftwindow, minfftratio, minsiglev, triggerfftratio, triggersiglev, tcoeff, zcoeff, flims
    def __init__(self, audiofile, timerange=[0,-1], lat=20, lon=-80, dropdate=dt.datetime.utcnow(), user_settings={}):
        
        #reading in WAV file
        self.audiofile = audiofile
        self.audiostream, self.f_s = readAXCPwavfile(audiofile, timerange)
        self.numpoints = len(self.audiostream)
        
        #initialize default settings, override user-specified ones
        self.init_AXCP_settings(user_settings)     
        
        #updating position, magnetic field components/declination
        self.lat = lat
        self.lon = lon
        self.dropdate = dropdate
        self.gm = gm.GeoMag(wmm_filename='WMM.COF')
        self.update_position()
        
        #initializing AXCP processor specific vars, as well as filter and conversion coefficients and output profile arrays
        self.initialize_AXCP_vars()
                
        
        
    def update_position(self):
        self.magvar = self.gm.get_params(dlat=self.lat, dlon=self.lon, time=self.dropdate)
        self.fh = self.magvar.fh
        self.fz = -self.magvar.fz
        self.dec = self.magvar.dec
        
        
        
    #initializes dict with default settings for AXCTD processor
    #any specified settings will be overwritten after this function is called in __init__
    def init_AXCP_settings(self, settings):
    
        self.settings = {}
        self.settings["refreshrate"] = 0.5
        self.settings["revcoil"]     = False
        self.settings["quality"]     = 3
                
        for csetting in settings:
            self.settings[csetting] = user_settings[csetting] #overwrite defaults for user specified settings 
          
        self.revcoil = self.settings["revcoil"]
        self.quality = self.settings["quality"]
        if self.quality <= 0:
            self.quality = 1
        elif self.quality >= 4:
            self.quality = 3
        
        self.refreshrate = self.settings["refreshrate"]
    
    
    
        
            
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
        self.pointsperloop = round(self.f_s * self.refreshrate)
        self.pointsperloop = ceil(self.pointsperloop / self.nss1) * self.nss1
        self.pointsperloop = ceil(self.pointsperloop / self.nss1 / self.nss2) * self.nss1 * self.nss2
        self.refreshrate = self.pointsperloop / self.f_s # actual seconds
        
        # first fit is from depth_beg to depth_beg + depth_chunk
        # next fit is depth_step deeper
        self.depth_beg = 7.5 
        self.depth_chunk = 5 
        self.depth_step  = 2
        
        #profile output arrays
        self.T       = []
        self.CCENV   = []
        self.PK      = []
        self.FCCDEV  = []
        self.FROTLP  = []
        self.FROTDEV = []
        self.TIME = [] 
        self.DEPTH= []
        self.TEMP = []
        self.U    = []
        self.V    = []
        self.ROTF = []
        self.ROTFRMS = []
        self.AREA = []
        self.EFBL = []
        self.CCBL = []
        self.FEFR = []
        self.FCCR = []
        self.VERR = []
        self.AERR = []
        self.TERR = []
        self.W    = []
        self.ENVCC       = []
        self.ENVCCRMS    = []
        self.PEAK        = []
        
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
        
        self.bxinlp,self.axinlp = signal.butter(4,3000/self.fnyq0, btype='lowpass',output='ba')  # low pass before first subsampling
        
        #bandpasses for compass coil, EF (velocity), and temperature bands
        ccband = np.array([2000,2500])/self.fnyq1
        efband = np.array([1000,1500])/self.fnyq1
        tband = np.array([250,500])/self.fnyq1
        self.bxccbp,self.axccbp = signal.butter(4,ccband, btype='bandpass',output='ba') # (rotation)
        self.bxefbp,self.axefbp = signal.butter(4,efband, btype='bandpass',output='ba') # (current magnitude)
        self.bxtebp,self.axtebp = signal.butter(4,tband, btype='bandpass',output='ba') # (temperature)
        
        self.bfzclp,self.afzclp = signal.butter(4,self.fzclp/self.fnyq1, btype='lowpass',output='ba') # low pass zero crossing pulses
        self.benvxcclp,self.aenvxcclp = signal.butter(2,1/self.fnyq1, btype='lowpass',output='ba') #nyquist frequency low pass filter for subsampling first round
        
        cc2band = np.array([2,20])/self.fnyq2
        self.bfccbp,self.afccbp = signal.butter(2,cc2band, btype='bandpass',output='ba') #bandpass filter for fcc
        self.benvfcclp,self.aenvfcclp = signal.butter(2,1/self.fnyq2, btype='lowpass',output='ba') #lowpass pilter for fcc (fN)
        self.bfrotlp,self.afrotlp = signal.butter(2,0.5/self.fnyq2, btype='lowpass',output='ba') #lowpass filter for rotation frequency
        
        frot2band = np.array([0.1,1.0])/self.fnyq2
        self.bfrotbp,self.afrotbp = signal.butter(2,frot2band, btype='bandpass',output='ba') #bandpass filter for rotation frequency
        self.benvfrotlp,self.aenvfrotlp = signal.butter(2,0.05/self.fnyq2, btype='lowpass',output='ba') 
        
        #z indices
        self.zxinlp = np.zeros(len(self.bxinlp)-1)
        
        self.zxccbp = np.zeros(len(self.bxccbp)-1)
        self.zxefbp = np.zeros(len(self.bxefbp)-1)
        self.zxtebp = np.zeros(len(self.bxtebp)-1)
        
        self.zfcclp = np.zeros(len(self.bfzclp)-1)
        self.zfeflp = np.zeros(len(self.bfzclp)-1)
        self.zftelp = np.zeros(len(self.bfzclp)-1)
        
        self.zenvxcclp = np.zeros(len(self.benvxcclp)-1)
        
        self.zfccbp = np.zeros(len(self.bfccbp)-1)
        self.zenvfcclp = np.zeros(len(self.benvfcclp)-1)
        self.zfrotlp = np.zeros(len(self.bfrotlp)-1)
        self.zfrotbp = np.zeros(len(self.bfrotbp)-1)
        self.zenvfrotlp = np.zeros(len(self.benvfrotlp)-1)
        
        
        
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
        self.tcalrs = (self.tcalfreq[0]*self.tcalres[0]-self.tcalfreq[1]*self.tcalres[1]) / ...
            (self.tcalfreq[1]-self.tcalfreq[0])
        self.tcalcap = 1.0/(4.4*self.tcalfreq[0]*(self.tcalrs+self.tcalres[0]))
        self.tcor      = [-0.10, -7.5e-5]
        
        # velocity scaling
        self.sfv = 1.0/(self.fz*1e-5*self.esep*(1.0+self.c1))
        self.sfw = self.fh/self.fz*(1.0+self.c2)/(1.0+self.c1)
        self.sfa = 100.0/(2.0*np.pi*self.fh*1e-5)
        
        #depth scaling
        self.depth_poly =  [4.68, 4.377, -0.00044   ] # Prater, JTech, Dec 1991
        self.inv_depth_poly = [3.11700848802192e-10,5.14227005928016e-06,0.228465938277901,-1.07390827937333]

        
        
                
    def first_subsample(self):
        pass
        
        
        
        
    def second_subsample(self):
        pass
        
        
        
        
        
    def calc_current_datapoint(self):
        pass
            
    
                    
                    
    #this function is called once per loop of the AXCP Processor and processes as much data as is available in the respective buffers 
    def iterate_AXCP_process(self, e): 
        
        npp += 1
        self.T.append(e/self.f_s) #time at end of current PCM chunk
        self.PK.append(max(abs(self.demod_buffer))) #peak audio value
        
        xinlp, self.zxinlp = signal.lfilter(self.bxinlp,self.axinlp,self.demod_buffer,zxinlp) #lowpass filter input buffer
        
        #running first subsample, applying filters, pulling big three frequency band zerocrossing points
        self.first_subsample()
        
        pklp = self.PK[-1] * np.ones(len(envxcclp)) #peak value, length of first subsample
        self.CCENV.append(envxcclp[-1] * np.sqrt(2)) #compass coil environment
        
        #running second subsample, pulling big three center frequencies for profile calculations
        self.second_subsample()
        
        #saving rotation frequencies
        self.FCCDEV.append( envfcclp[-1]*np.sqrt(2) )
        self.FROTLP.append( frotlp[-1] )
        self.FROTDEV.append( envfrotlp[-1]*np.sqrt(2) )
        
        
        #only retain last 40 points from previous iteration, append on new values
        retainind = -40 #last 40 points
        self.pk = np.append(self.pk[retainind:], pk_new)
        self.envxcc = np.append(self.envxcc[retainind:], envxcc_new)
        self.fcc = np.append(self.fcc[retainind:], fcc_new)
        self.fef = np.append(self.fef[retainind:], fef_new)
        self.fte = np.append(self.fte[retainind:], fte_new)
        self.tim = np.append(self.tim[retainind:], tim_new)
        
        
        #checking to see if probe has spunup (update status to 1 if so)
        #requirements: time >= 5 seconds, rotation frequency deviation below 0.5 (steady state), rotation frequency between 10 and 20 Hz
        #precise spinup time is defined as achieving 8 Hz, half rotation rate of average 16 Hz spin
        if not self.status and self.T[-1] >= 5 and self.FROTDEV[-1] <= 0.5 and 10 < self.FROTLP < 20:
            
            self.status = 1 #noting profile has spun up,  finding precise spinup time
            
            #getting 0-centered, recent, good compass coil frequencies for 0-crossing analysis
            fcc0 = self.fcc - np.mean(np.array([cfcc for i,cfcc in enumerate(self.fcc) if self.tim[-1]-self.tim[i] <= 5 and not np.isnan(cfcc) and np.isfinite(cfcc)])) 
            r = [1 if cfcc0 >= 0 else -1 for cfcc0 in fcc0] #ID compass coil frequency zero crossings
            cross_points = np.where(np.diff(r) > 0)[0]
            
            #getting calculating frequency rate of change across zero crossing points to interpolate exact zero crossing times
            d_freq = fcc0[cross_points+1] - fcc0[cross_points]
            d_time = self.tim[cross_points+1] - self.tim[cross_points]
            tim_zc = self.tim3[cross_points] - fcc0[cross_points] * d_time / d_freq #actual zerocross times
            rotf_zc = 1 / np.diff(tim_zc)
            
            #rotf_zc will start extremely high when there is no valid data, drop low for a second as the static compass coil frequency is detected, and then quickly increase above 8 Hz as the probe spins up. Spinup time is the first time after the point where the probe speed first exceeds 8 Hz
            spinup_ind = np.where(rotf_zc < 8)[0]
            if len(spinup_ind) > 0:
                try:
                    self.tspinup = tim_zc[spinup_ind + 2]
                except IndexError:
                    self.tspinup = tim_zc[-1]
            else:
                self.tspinup = tim_zc[0]
                
                
        
        
        
        
        
    def run(self):
    
        
        self.maxtime = self.numpoints/self.f_s
        
            
        # setting up thread while loop- terminates when user clicks "STOP" or audio file finishes processing
        i = -1
        
        self.keepgoing = True
        
        # initialize state- probe hasn't spun up yet
        self.status = 0
        
        #initialize self.demodbufferstartind
        self.demodbufferstartind = 0
        
        
        #MAIN PROCESSOR LOOP
        while self.keepgoing:
            i += 1
            
                
            #kill threads once time exceeds the max time of the audio file
            #NOTE: need to do this on the cycle before hitting the max time when processing from audio because the WAV file processes faster than the thread can kill itself
            
            #calculating end of next slice of PCM data for signal level calcuation and demodulation
            e = self.demodbufferstartind + self.pointsperloop
            
            if self.numpoints - self.demodbufferstartind < 4*self.N_power: #kill process at file end
                self.keepgoing = False
                print("[+] Processing status: 100%")
            
            elif e >= self.numpoints: #terminate loop if at end of file
                e = self.numpoints - 1
                self.keepgoing = False
            
            #add next round of PCM data to buffer for signal calculation and demodulation
            # self.demod_buffer = np.append(self.demod_buffer, self.audiostream[self.demodbufferstartind:e])
            self.demod_buffer = self.audiostream[self.demodbufferstartind:e]
                
            
            print(f"[+] Processing status: {round(100*self.demodbufferstartind/self.numpoints)}%         ", end='\r')
            
            if self.keepgoing: #only process buffer if there is enough data
                oldstatus = self.status #track status change to emit triggered signal when necessary
                
                #demodulating and parsing current batch of AXCTD PCM data
                data = self.iterate_AXCP_process(e)
                
                
                if len(data) > 1: #retrieved good profile data
                    pass #save profile data here
                
    
                #incrementing demod buffer start forward
                self.demodbufferstartind = e 
                    
                    
                    
            #sleeping until ready to process more data 
            timemodule.sleep(0.001) #slight pause to free some resources when processing from audio
        
        
        
    