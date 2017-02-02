#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 09:54:08 2017

@author: Malfatti

Modified from JRodrigoF's code, avaliable at
https://github.com/JRodrigoF/AVrecordeR/blob/master/AVrecordeR.py

Record audio and video
"""
#%%
import cv2
import pyaudio
import wave
import threading
import time
import subprocess
#import os

from datetime import datetime

########################
## JRF
## VideoRecorder and AudioRecorder are two classes based on openCV and pyaudio, respectively. 
## By using multithreading these two classes allow to record simultaneously video and audio.
## ffmpeg is used for muxing the two signals
## A timer loop is used to control the frame rate of the video recording. This timer as well as
## the final encoding rate can be adjusted according to camera capabilities
##

########################
## Usage:
## 
## numpy, PyAudio and Wave need to be installed
## install openCV, make sure the file cv2.pyd is located in the same folder as the other libraries
## install ffmpeg and make sure the ffmpeg .exe is in the working directory
##
## 
## start_AVrecording(filename) # function to start the recording
## stop_AVrecording(filename)  # "" ... to stop it
##
##
########################

Now = datetime.now().strftime("%Y%m%d%H%M%S")
Info = {'FPS': 23,
        'VideoFormat': 'MJPG',
        'Res': (1328,966),
        'VideoFile': 'Video-' + Now + '.avi',
        'Rate': 192000,
        'Buffer': 384,
        'ChNo': 2,
        'AudioFormat': pyaudio.paInt32,
        'AudioFile': 'Rec-' + Now + '.wav'}

class VideoRecorder():
	
	# Video class based on openCV 
	def __init__(self, Info):
		
		self.open = True
		self.device_index = 0
		self.fps = Info['FPS']               # fps should be the minimum constant rate at which the camera can
		self.fourcc = Info['VideoFormat']       # capture images (with no decrease in speed over time; testing is required)
		self.frameSize = Info['Res'] # video formats and sizes also depend and vary according to the camera used
		self.video_filename = Info['VideoFile']
		self.video_cap = cv2.VideoCapture(self.device_index)
		self.video_writer = cv2.VideoWriter_fourcc(*self.fourcc)
		self.video_out = cv2.VideoWriter(self.video_filename, self.video_writer, self.fps, self.frameSize)
		self.frame_counts = 1
		self.start_time = time.time()

	
	# Video starts being recorded 
	def record(self, Info):
		
#		counter = 1
#		timer_start = time.time()
#		timer_current = 0
		
		
		while(self.open==True):
			ret, video_frame = self.video_cap.read()
			if (ret==True):
				
					self.video_out.write(video_frame)
#					print str(counter) + " " + str(self.frame_counts) + " frames written " + str(timer_current)
					self.frame_counts += 1
#					counter += 1
#					timer_current = time.time() - timer_start
					time.sleep(1/Info['FPS'])
					
					# Uncomment the following three lines to make the video to be
					# displayed to screen while recording
					
					gray = cv2.cvtColor(video_frame, cv2.COLOR_BGR2GRAY)
					cv2.imshow('video_frame', gray)
					cv2.waitKey(1)
			else:
				break
				
				

	# Finishes the video recording therefore the thread too
	def stop(self):
		
		if self.open==True:
			
			self.open=False
			self.video_out.release()
			self.video_cap.release()
			cv2.destroyAllWindows()
			
		else: 
			pass


	# Launches the video recording function using a thread			
	def start(self):
		video_thread = threading.Thread(target=self.record(Info))
		video_thread.start()





class AudioRecorder():
	

    # Audio class based on pyAudio and Wave
    def __init__(self, Info):
        
        self.open = True
        self.rate = Info['Rate']
        self.frames_per_buffer = Info['Buffer']
        self.channels = Info['ChNo']
        self.format = Info['AudioFormat']
        self.audio_filename = Info['AudioFile']
        self.audio = pyaudio.PyAudio()
        self.stream = self.audio.open(format=self.format,
                                      channels=self.channels,
                                      rate=self.rate,
                                      input=True,
                                      frames_per_buffer = self.frames_per_buffer)
        self.audio_frames = []


    # Audio starts being recorded
    def record(self):
        
        self.stream.start_stream()
        while(self.open == True):
            data = self.stream.read(self.frames_per_buffer) 
            self.audio_frames.append(data)
            if self.open==False:
                break
        
            
    # Finishes the audio recording therefore the thread too    
    def stop(self):
       
        if self.open==True:
            self.open = False
            self.stream.stop_stream()
            self.stream.close()
            self.audio.terminate()
               
            waveFile = wave.open(self.audio_filename, 'wb')
            waveFile.setnchannels(self.channels)
            waveFile.setsampwidth(self.audio.get_sample_size(self.format))
            waveFile.setframerate(self.rate)
            waveFile.writeframes(b''.join(self.audio_frames))
            waveFile.close()
        
        pass
    
    # Launches the audio recording function using a thread
    def start(self):
        audio_thread = threading.Thread(target=self.record)
        audio_thread.start()



def start_video_recording(Info):
				
	global video_thread
	
	video_thread = VideoRecorder(Info)
	video_thread.start()

	return Info['VideoFile']
	

def start_audio_recording(Info):
				
	global audio_thread
	
	audio_thread = AudioRecorder(Info)
	audio_thread.start()

	return Info['AudioFile']


def start_AVrecording(Info):
				
	global video_thread
	global audio_thread
	
	video_thread = VideoRecorder(Info)
	audio_thread = AudioRecorder(Info)

	audio_thread.start()
	video_thread.start()

	return Info['AVFile']


def stop_AVrecording(Info):
	
	audio_thread.stop() 
	frame_counts = video_thread.frame_counts
	elapsed_time = time.time() - video_thread.start_time
	recorded_fps = frame_counts / elapsed_time
	print("total frames " + str(frame_counts))
	print("elapsed time " + str(elapsed_time))
	print("recorded fps " + str(recorded_fps))
	video_thread.stop() 

	# Makes sure the threads have finished
	while threading.active_count() > 1:
		time.sleep(1)

	
#	 Merging audio and video signal
	
	if abs(recorded_fps - Info['FPS']) >= 0.01:    # If the fps rate was higher/lower than expected, re-encode it to the expected
										
		print("Re-encoding")
		cmd = "ffmpeg -r " + str(recorded_fps) + " -i temp_video.avi -pix_fmt yuv420p -r 6 temp_video2.avi"
		subprocess.call(cmd, shell=True)
	
#		print("Muxing")
#		cmd = "ffmpeg -ac 2 -channel_layout stereo -i temp_audio.wav -i temp_video2.avi -pix_fmt yuv420p " + filename + ".avi"
#		subprocess.call(cmd, shell=True)
#	
#	else:
#		
#		print "Normal recording\nMuxing"
#		cmd = "ffmpeg -ac 2 -channel_layout stereo -i temp_audio.wav -i temp_video.avi -pix_fmt yuv420p " + filename + ".avi"
#		subprocess.call(cmd, shell=True)
#
#		print ".."


## Required and wanted processing of final files
#def file_manager(filename):
#
#	local_path = os.getcwd()
#
#	if os.path.exists(str(local_path) + "/temp_audio.wav"):
#		os.remove(str(local_path) + "/temp_audio.wav")
#	
#	if os.path.exists(str(local_path) + "/temp_video.avi"):
#		os.remove(str(local_path) + "/temp_video.avi")
#
#	if os.path.exists(str(local_path) + "/temp_video2.avi"):
#		os.remove(str(local_path) + "/temp_video2.avi")
#
#	if os.path.exists(str(local_path) + "/" + filename + ".avi"):
#		os.remove(str(local_path) + "/" + filename + ".avi")
#	
#	
#
#	
#if __name__== "__main__":
#	
#	filename = "Default_user"	
#	file_manager(filename)
#	
#	start_AVrecording(filename)  
#	
#	time.sleep(10)
#	
#	stop_AVrecording(filename)
#	print "Done"
#


