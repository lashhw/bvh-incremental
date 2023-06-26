import sys
import numpy as np

if len(sys.argv) != 4:
    print('usage: python3 sample.py [num_samples] [input_file] [output_file]')
    exit(1)

num_samples, input_file, output_file = sys.argv[1:]
num_samples = int(num_samples)

data = np.fromfile(input_file, dtype=np.float32)
data = data.reshape(-1, 7)
sample_idx = np.random.choice(data.shape[0], num_samples, replace=False)
sampled = data[sample_idx]
sampled.tofile(output_file)
