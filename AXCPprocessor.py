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

import numpy as np
from scipy import signal
from scipy.io import wavfile #for wav file reading

import time as timemodule
import datetime as dt

from traceback import print_exc as trace_error

import geomag as gm


#TODO:
#   > switch temperature conversion from zero crossing to FFT
#   > switch filter types from b/a to sos (output='sos' in butter and using sosfilt instead of lfilter)


#this could be problematic for debugging but I'm getting tired of the "chunk could not be understood" messages
import warnings
warnings.filterwarnings("ignore")


#conversion: coefficients=C,  D_out = C[0] + C[1]*D_in + C[2]*D_in^2 + C[3]*D_in^3 + ...
def dataconvert(data_in,coefficients):
    
    datatype = 1 #integer or float
    if type(data_in) == list:
        datatype = 2
    elif type(data_in) == np.ndarray: #numpy array
        dataype = 3
        
    if datatype == 1:
        data_in = [data_in]
        
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
    
    #importing necessary functions to handle AXCP processing (automatically attaches them to self)
    from ._AXCP_decode_fxns import (init_AXCP_settings, initialize_AXCP_vars, init_filters, init_constants, first_subsample, second_subsample, calc_current_datapoint, iterate_AXCP_process, refine_spindown_prof, calculate_true_velocities)

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
        self.fh = self.magvar.fh #always positive
        self.fz = -self.magvar.fz #switches convention so positive is up
        self.dec = self.magvar.dec #positive is East
        
        
        
        
        
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
            
            
        #after finishing entire profile, refine the spindown point, correct amean, and adjust the profile
        print("[+] Processing status: 100%         ")
        if self.status and len(self.TIME) > 0: #spinup detected, valid profile points recorded
            self.refine_spindown_prof()
            self.calculate_true_velocities()
        else:
            print(f"[!] Audio file processed- spinup not detected")
            
            
        
        
        
    