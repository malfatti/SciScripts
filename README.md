# SciScripts
Scripts written to control devices/run experiments/analyze data

These are codes we use in the lab to:  

'''
.  
├── Arduino  
│   └── Sketchbooks  
│       ├── ControlArduinoWithSoundBoard  
│       │   └── ControlArduinoWithSoundBoard.ino  
│       │	  Read an analog in to write in a digital out when voltage is  
│       │	  above a threshold. This provide a way to control digital 5V  
│       │	  devices (like a laser) with an analog device (like a sound   
│       │	  board)    
│       └── LaserStimAndSoundTTL  
│           └── LaserStimAndSoundTTL.ino  
│		  Write in a digital out after a serial command. This contains  
│		  several different protocols, each can be started by a serial  
│		  command.  
│  
│  
├── Octave  
│   % These scripts may not fully work in Octave yet, I'm adapting them from  
│   % Matlab.  
│   ├── Deps  
│   │   └── These are functions needed for the other Octave scripts to work.  
│   ├── Intan  
│   │   ├── AnalyseSpikesIntan.m  
│   │   │     Manipulate data recorded using INTAN RHA2000 software with script  
│   │   │     SoundLaserStim.m for stimulation.  
│   │   └── AnalyseSpikesIntan-OneTTL.m  
│   │         Manipulate data recorded using INTAN RHA2000 software with script  
│   │         SoundLaserStim2Lasers.m for stimulation.  
│   ├── OpenEphys   
│   │   └── AnalyseOpenEphys.m   
│   │         Manipulate data recorded using OpenEphysGUI with script   
│   │         ControlSoundBoard.py for stimulation. For now this is a   
│   │         modification from the Octave code that manipulates data.  
│   └── SoundLaserStimulation  
│       ├── SoundLaserStim2Lasers.m  
│       │     Generate sound and light stimulation pulses, narrow-band white  
│       │     noise for sound and square pulses for laser/led.  
│       ├── SoundLaserStim.m  
│       ├── SoundStimGCamp.m  
│       │     Generate sound stimulation and read an analog in (spectrum power  
│       │     meter) signal.  
│       └── SoundStimGCampSingleChannel.m  
│             Generate sound stimulation and read an analog in (spectrum power   
│             meter) signal.  
│  
│  
├── Python2  
│   % Sketch scripts to analyze data recorded in OpenEphysGUI with script  
│   % ControlSoundBoard.py for stimulation.  
│   └── Deps  
│       └── These are functions needed for the other python scripts to work.  
│  
│  
└── Python3  
    ├── AlsaControl  
    │   └── AcousticNoiseTrauma.py  
    │	      This will generate a white noise and send the signal to the  
    │	      computer's sound board.  
    ├── ControlSoundBoard.py  
    │     This code generates a stereo audio signal, where the signal of  
    │     the right channel is a narrow-band white noise and of the left  
    │     channel is *almost* square pulses.  
    └── GenerateChaoticMelody.py  
	  This is a code to generate chaotic fractal phase space trajectories  
	  and map them in notes. The final output is a Lilypond code txt file.  
'''


___
I'm Thawann Malfatti, PhD student from the Neurodynamics Lab of the Brain  
Institute of the Federal University of Rio Grande do norte. Our team of "code  
writers" is constituted by:

Richardson Leão, MD, PhD;  
George Nascimento, PhD;  
Helton Maia, PhD;  
José Targino, MSc;  
Thawann Malfatti, MSc;  
Rafael Franzon, BSc.
