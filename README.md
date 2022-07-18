# **AXCP Processor**

## Overview
The AXCP Processor can reprocess AXCP audio (WAV) files via command line, outputting a file containing AXCP signal data (signal level, peak frequency) and observed temperature versus depth. 


### Usage:
`python3 processAXCP -i inputfile.WAV`

### Installation and Setup:
This script requires python modules other than python base. Install them with `pip install -r requirements.txt`

### Optional flags:

<table>
  <tbody>
    <tr>
      <th align="center">Flag</th>
      <th align="left">Purpose</th>
    </tr>
    <tr>
      <td align="center">-o</td>
      <td>Output filename- full or relative path and filename for output AXCP drop metadata and profile file (defaults to <code>output.txt</code>)</td>
    </tr>
    <tr>
      <td align="center">-s</td>
      <td>Start time of AXCP profile in WAV file, format: SS, MM:SS, or HH:MM:SS (default: 0 sec)</td>
    </tr>
    <tr>
      <td align="center">-e</td>
      <td>End time of AXCP profile in WAV file, format: SS, MM:SS, or HH:MM:SS (default: end of WAV file)</td>
    </tr>  
    <tr>
      <td align="center">-y</td>
      <td>Latitude (N > 0) to be used for magvar calculations</td>
    </tr>  
    <tr>
      <td align="center">-x</td>
      <td>Longitude (E > 0) to be used for magvar calculations</td>
    </tr>  
    <tr>
      <td align="center">-d</td>
      <td>Date (YYYYMMDD) to be used for magvar calculations</td>
    </tr>  
    <tr>
      <td align="center">-q</td>
      <td>Processing quality (1 = high quality/slow, 2 = moderate quality/speed, 3 = lower quality/higher speed)</td>
    </tr>  
    <tr>
      <td align="center">-r</td>
      <td>Whether or not probe coil was reversed (1 = true, default 0)</td>
    </tr>  
    <tr>
      <td align="center">-p</td>
      <td>Refresh rate (seconds) of AXCP iteration loop</td>
    </tr>  
    <tr>
      <td align="center">-u</td>
      <td>Whether or not to try to detect probe spindown realtime and automatically stop processing (1=yes,0=no)</td>
    </tr>  
    <tr>
      <td align="center">-m</td>
      <td>Whether to use zero crossings (1) or FFT (2) to calculate temperature at each depth</td>
    </tr>  
    <tr>
      <td align="center">-w</td>
      <td>Temperature FFT window (sec)- only applicable for temperature_mode=2, must be less than or equal to refreshrate.</td>
    </tr>
    </tbody>
</table>

<br />
<br />



## More Information

AXCP Processor uses the same algorithms implemented in the Airborne eXpendable Buoy Processing System (AXBPS) to process AXCP profiles from raw audio files. For more information, see the [AXBPS homepage](http://mmmfire.whoi.edu/axbps).


<br />

