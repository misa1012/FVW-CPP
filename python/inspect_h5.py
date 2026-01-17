import h5py
import sys

filename = sys.argv[1]

with h5py.File(filename, "r") as f:
    print("Keys:", list(f.keys()))
    if 'timesteps' in f:
        print("Timesteps:", list(f['timesteps'].keys()))
    
    if 'config' in f:
        print("Config found")
