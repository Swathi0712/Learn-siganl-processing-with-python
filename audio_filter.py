# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.io import wavfile
import os
import argparse

class AudioFilter:
    def __init__(self):
        self.sampleRate = None
        self.audioData = None
        self.filteredData = None
        self.duration = None
        self.timeAxis = None
        
    def generateTestTone(self, duration = 3.0, frequencies = [440, 1000, 3000]):
        self.sampleRate = 44100
        self.duration = duration
        numSamples = int(self.sampleRate * duration)
        self.timeAxis = np.linspace(0, duration, numSamples, False)
    
        self.audioData = np.zeros(numSamples)
        for i, freq  in enumerate(frequencies):
            amplitude = 1.0/(i+1)
            self.audioData += amplitude * np.sin(2 * np.pi * freq *self.timeAxis)
        
        
        self.audioData = 0.5 *self.audioData/np.max(np.abs(self.audioData))
        
    def loadAudio(self, filePath):
        try:
            self.sampleRate, self.audioData = wavfile.read(filePath)
            
            if len(self.audioData.shape) > 1:
                self.audioData = np.mean(self.audioData, axis=1)
            
            if self.audioData.dtype!= np.float32 and self.audioData.dtype!= np.float64:
                self.audioData = self.audio.astype(np.float32)/np.max(np.abs(self.audioData))
            
            self.duration = len(self.audioData)/self.sampleRate
            self.timeAxis = np.linspace(0,self.duration,len(self.audioData))
            
        except Exception as e:
            print(f"error from load audio{e}")
            return False
        
        return True
    
    def applyFilter(self, filterType, **params):
        if self.audioData is None:
            print("No audii data found")
            return False
        
        order = params.get('order', 4)
        nyquist = 0.5 * self.sampleRate
        
        if filterType == 'lowpass':
            cutoff = params.get('cutoff', 1000)
            normalCutoff = cutoff/nyquist
            
            designMethod = params.get('design', 'butter')
            if designMethod == 'butter':
                b, a = signal.butter(order, normalCutoff, btype = 'low')
                
        method = params.get('method', 'filtfilt')
        if method == 'filtfilt':
            self.filteredData = signal.filtfilt(b, a , self.audioData)
            
        return True
                
    
    def saveAudio(self, outputFile):
        if self.filteredData is None:
            print("no filtered data found")
            return False
        else:
            print(self.filteredData)
            
        normalizedData = np.int16(self.filteredData * 32767 / np.max(np.abs(self.filteredData)))
        
        try:
            wavfile.write(outputFile, self.sampleRate, normalizedData)
        except Exception as e:
            return False
                
            
def main():

    parser = argparse.ArgumentParser(description='Audio Filter Application')
    parser.add_argument('--input', '-i', help='Input WAV file path')
    parser.add_argument('--output', '-o', default='filtered_output.wav', help='Output WAV file path')
    parser.add_argument('--filter', '-f', default='lowpass', 
                        choices=['lowpass', 'highpass', 'bandpass', 'bandstop', 'notch', 'fir'],
                        help='Filter type')
    parser.add_argument('--cutoff', type=float, help='Cutoff frequency for lowpass/highpass (Hz)')
    parser.add_argument('--low', type=float, help='Low cutoff for bandpass/bandstop (Hz)')
    parser.add_argument('--high', type=float, help='High cutoff for bandpass/bandstop (Hz)')
    parser.add_argument('--notch', type=float, help='Notch frequency (Hz)')
    parser.add_argument('--order', type=int, default=4, help='Filter order')
    parser.add_argument('--design', default='butter', 
                        choices=['butter', 'cheby1', 'ellip', 'bessel', 'fir'],
                        help='Filter design method')
    
    args = parser.parse_args()
    print("Creating audio filter object")
    # Create the audio filter object
    audio_filter = AudioFilter()
    
    print("Checking for inputs")
    if args.input:
        if not audio_filter.loadAudio(args.input):
            return
        
    
    print("filterparams")    
    filterParams = {'order':args.order, 'design':args.design}
    
    
     
    if args.filter in ['lowpass', 'highpass']:
        if args.cutoff:
            filterParams['cutoff'] = args.cutoff
        else:
            filterParams['cutoff'] = 1000 
      
    print("applying audio filter")
    if not audio_filter.applyFilter(args.filter, **filterParams):
        return
    
    print("Saving audio file")
    audio_filter.saveAudio(args.output)
     
    print("success")
       
    
 

if __name__  == "__main__":
    print("hi")
    main()
    
    
        
