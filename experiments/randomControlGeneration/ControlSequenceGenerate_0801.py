import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.io import savemat
from datetime import datetime

# 时间向量
t = np.linspace(0, 240, 24000)  # 60秒，采样率100 Hz

# 安全频率区间（rad/s）
frequencies_1 = [0.5, 1.5, 3.0]      # for u1(t) < ω_safe
frequencies_2 = [0.8, 2.2, 3.5]      # for u2(t) < ω_safe
amplitudes = [0.05, 0.03, 0.02]      # 控制幅度小，避免非线性激发
phases = np.random.uniform(0, 2*np.pi, size=3)

# 构造控制输入
u1_raw = sum(a * np.sin(w * t + p) for a, w, p in zip(amplitudes, frequencies_1, phases))
u2_raw = sum(a * np.sin(w * t + p) for a, w, p in zip(amplitudes, frequencies_2, phases))

# First normalize to [-1, 1] range
u1 = u1_raw / np.max(np.abs(u1_raw))  # Normalize to [-1, 1]
u2 = u2_raw / np.max(np.abs(u2_raw))  # Normalize to [-1, 1]

# Apply tapering to ensure smooth start and end at origin
taper_length = int(0.05 * len(t))  # 5% of total length for tapering
taper_window = np.ones_like(t)

# Create smooth tapering window (cosine taper)
for i in range(taper_length):
    taper_window[i] = 0.5 * (1 - np.cos(np.pi * i / taper_length))
    taper_window[-(i+1)] = 0.5 * (1 - np.cos(np.pi * i / taper_length))

# Apply tapering to force start and end to zero
u1 = u1 * taper_window
u2 = u2 * taper_window

# Scale to ensure the trajectory stays within [-1, 1] range
max_abs_u1 = np.max(np.abs(u1))
max_abs_u2 = np.max(np.abs(u2))
max_abs_overall = max(max_abs_u1, max_abs_u2)

# Scale both signals by the same factor to maintain trajectory shape
if max_abs_overall > 1.0:
    scale_factor = 0.95 / max_abs_overall  # Use 0.95 to ensure we stay within [-1,1]
    u1 *= scale_factor
    u2 *= scale_factor

# Verify that start and end points are at origin
print(f"Initial point: u1[0]={u1[0]:.6f}, u2[0]={u2[0]:.6f}")
print(f"End point: u1[-1]={u1[-1]:.6f}, u2[-1]={u2[-1]:.6f}")
print(f"u1 range: [{np.min(u1):.6f}, {np.max(u1):.6f}]")
print(f"u2 range: [{np.min(u2):.6f}, {np.max(u2):.6f}]")
print(f"Maximum absolute value in trajectory: {max(np.max(np.abs(u1)), np.max(np.abs(u2))):.6f}")

u = np.stack([u1, u2], axis=1)  # shape: [len(t), 2]

N = len(t)
dt = t[1] - t[0]
freqs = fftfreq(N, dt)

U1 = np.abs(fft(u1)) / N
U2 = np.abs(fft(u2)) / N

plt.figure(figsize=(10, 4))
plt.plot(freqs[:N//2], U1[:N//2], label='u1 spectrum')
plt.plot(freqs[:N//2], U2[:N//2], label='u2 spectrum')
plt.axvline(3.0, color='r', linestyle='--', label='cutoff ~3Hz')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)
plt.title("Control Input Spectrum")

# Create a figure with 3 subplots
plt.figure(figsize=(15, 5))

# Plot u1(t) vs time
plt.subplot(131)
plt.plot(t, u1, 'b-', label='u1(t)')
plt.xlabel('Time (s)')
plt.ylabel('u1')
plt.grid(True)
plt.legend()
plt.title('Control Input u1(t)')

# Plot u2(t) vs time
plt.subplot(132)
plt.plot(t, u2, 'r-', label='u2(t)')
plt.xlabel('Time (s)')
plt.ylabel('u2')
plt.grid(True)
plt.legend()
plt.title('Control Input u2(t)')

# Plot u1 vs u2 (control input space)
plt.subplot(133)
plt.plot(u1, u2, 'k-', alpha=0.6)
plt.plot(u1[0], u2[0], 'go', label='Start')  # Mark start point
plt.plot(u1[-1], u2[-1], 'ro', label='End')  # Mark end point
plt.xlabel('u1')
plt.ylabel('u2')
plt.grid(True)
plt.legend()
plt.title('Control Input Space (u1 vs u2)')
plt.axis('equal')  # Make axes equal scale

plt.tight_layout()
plt.show()

# Save the final trajectory data
trajectory_data = {
    'x': u1,
    'y': u2,
    'time': t
}

# Save the final trajectory data
current_date = datetime.now().strftime("%Y%m%d")  # Format: YYYYMMDD
filename = f"trajectory_data_{current_date}_v1.mat"  # Generate filename with current date
savemat(filename, trajectory_data)

print(f"Data saved to {filename}")