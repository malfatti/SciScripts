# SciScripts
Scripts written to control devices/run experiments/analyze data

[Still uploading]

These are codes that we use in the lab to:

## [Arduino]

Read an analog in to write in a digital out when voltage is above a threshold;
  - This provide a way to control digital 5V devices (like a laser) with an  
    analog device (like a soundboard)

Write in a digital out after a serial command.
  - This contains several different protocols, each can be started by a serial  
    command

## [Octave]

  % These scripts may not fully work in Octave yet, I'm adapting them from  
  % Matlab.
  
Manipulate data recorded using INTAN RHA2000 software;
  - Filter data for spike detection;
  - Visualize and select channels that contains spikes;
  - Detect and clusterize spikes;
  - Visualize spike waveforms and select "fake clusters" for removal;
  - Generate PSTHs and ABR waveforms.

Generate sound and light stimulation pulses
  - Narrow-band white noise for sound and square pulses for laser/led.

## [Python]

Generate sound signal to control an Arduino board using the computer's sound  
board;
  - This code generates a stereo audio signal, where the signal of the right  
    channel is a narrow-band white noise and of the left channel is a *almost*  
    square pulse.

Generate sound signal for acoustic noise trauma;
  - This will generate a white noise and send the signal to the computer's sound  
    board.

Manipulate data recorded using OpenEphys-GUI;
  - For now this is a modification to the Octave code that manipulates data  
    recorded using INTAN RHA2000 software.

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
