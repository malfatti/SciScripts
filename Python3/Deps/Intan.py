#! /bin/env python
#
# Max Hodak, 2 November 2010
#
# Fixed, modified to increase speed using numpy, Alex Chubykin, 06 June 2012

import sys, struct, math
import numpy

def ReadData(fname, verbose = False, describe_only = False):
  """
  Example:
  
    from intan import ReadData
    mydata = ReadData("spikes_101029_011129.int")
  
  mydata looks like:
    {
      "duration": 60,                             # seconds
      "channel_count": 17,                        # number of amps
      "analog": [((-0.2594, -0.1502, ...), ...],  # voltage data
    }
  
  len(mydata['analog']) is approximately mydata['duration']*25000 and 
  len(mydata['analog'][i]) = mydata['channel_count'] for i < mydata['duration']*25000
  (ish, unless the fpga was dropping frames or such such).  ReadData always returns a dict
  for extensability reasons.  For example, digital or other lines may be added later.
  If the int file is of an improper format, or something else goes wrong, an exception is thrown.
  """
  
  fh = open(fname, 'rb')
  data=fh.read()
  print('data length',len(data))
  header = data[0:3]
  
  (version, major, minor) = struct.unpack("<BBB", data[0:3])
  
  if version != 128:                         raise Exception("Improper data file format.")
  if (major != 1 or minor != 1) and verbose: print("Datafile may not be compatible with this script.")
  
  num_amps = sum([struct.unpack("<B", data[i+2:i+3])[0] for i in range(1,64)])
  t_count = ((len(data) - 67) / (1+4*num_amps))
  t_max = t_count / 25000
  
  description = {"duration": t_max, "channel_count": num_amps}
  
  if verbose:       print('Data file contains %0.2f seconds of data from %d amplifier channel(s)' % (t_max, num_amps))
  if describe_only: return description
  
  data = data[67:] # Throw away the header.
  chunk_offset = 4*num_amps+1
  analog=[]

  for l in range(0,len(data),chunk_offset):
      analog_temp=numpy.array(struct.unpack("<"+"f"*num_amps+"x", data[l:l+chunk_offset]))
      analog.append(analog_temp)

  analog=numpy.array(analog) # using numpy arrays significantly increases the speed and improves memory allocation
  
  t = numpy.array([i/25000.0 for i in range(0, int(t_count))])
  
  # Re-scanning the entire dataset seems very inefficient, but benchmarking says it's not a practical
  # issue (adds 0.2 seconds at most out of a ~3.8 second total for 60 seconds of 17 channels).
  aux=[]
  for i in range(int(len(data)/chunk_offset)):
      aux_temp=numpy.array(struct.unpack("<b",data[i*chunk_offset+(chunk_offset-1):i*chunk_offset+chunk_offset])[0])
      aux.append(aux_temp)
  aux=numpy.array(aux)
  
  description.update({'analog': analog, 'aux': aux, 't': t})
  return description


if __name__ == '__main__':
    a=ReadData(sys.argv[1], verbose = True, describe_only = False) # a is a dictionary with the data, a["analog"] returns the 2xD numpy array with the recordings
    print(a)
	
