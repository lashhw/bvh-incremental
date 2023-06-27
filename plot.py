import numpy as np
import matplotlib.pyplot as plt

full_trv = np.fromfile('full_trv.bin', dtype=np.uint)
full_ist = np.fromfile('full_ist.bin', dtype=np.uint)
reduced_trv = np.fromfile('reduced_trv.bin', dtype=np.uint)
reduced_ist = np.fromfile('reduced_ist.bin', dtype=np.uint)
reduced_trv_full = np.fromfile('reduced_trv_full.bin', dtype=np.uint)
reduced_trv_reduced = np.fromfile('reduced_trv_reduced.bin', dtype=np.uint)

print(f'full_trv avg: {np.mean(full_trv)}')
print(f'full_ist avg: {np.mean(full_ist)}')
print(f'reduced_trv avg: {np.mean(reduced_trv)}')
print(f'reduced_ist avg: {np.mean(reduced_ist)}')
print(f'reduced_trv_full avg: {np.mean(reduced_trv_full)}')
print(f'reduced_trv_reduced avg: {np.mean(reduced_trv_reduced)}')

# plt.hist(full_trv, 100, density=True, alpha=0.5)
# plt.hist(reduced_trv, 100, density=True, alpha=0.5)
# plt.show()
