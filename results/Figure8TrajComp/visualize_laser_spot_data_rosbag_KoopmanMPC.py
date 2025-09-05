from mcap.reader import make_reader
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from datetime import datetime
import os
from scipy.io import savemat, loadmat

def parse_mcap_file(mcap_file, topic):
    """
    Parse an MCAP file and extract messages from a specific topic.

    Args:
        mcap_file (str): Path to the MCAP file.
        topic (str): Topic to extract messages from.

    Returns:
        list: Extracted messages data
    """
    extracted_data = []  # Store extracted messages
    try:
        with open(mcap_file, "rb") as f:
            reader = make_reader(f)
            print(f"Parsing MCAP file: {mcap_file}")
            print(f"Extracting messages from topic: {topic}")
            
            for schema, channel, message in reader.iter_messages():
                if channel.topic == topic:
                    extracted_data.append({
                        "time": message.log_time,
                        "publish_time": message.publish_time,
                        "data": message.data.hex()  # Convert bytes to hex string
                    })
        
        print(f"Extracted {len(extracted_data)} messages from topic '{topic}'")
        return extracted_data
        
    except Exception as e:
        print(f"Error parsing MCAP file: {e}")
        return []

def parse_actuator_motion_data(data_source, source_type="json"):
    """
    Parse actuator motion data and extract pixel and actual coordinates.
    
    Args:
        data_source (str): Path to the data file or raw data list.
        source_type (str): Type of data source ("json" or "mcap" or "raw")
    
    Returns:
        tuple: (timestamps, pixel_coords, actual_coords)
    """
    try:
        # Load data based on source type
        if source_type == "json":
            with open(data_source, 'r') as f:
                data = json.load(f)
        elif source_type == "mcap":
            # Extract topic and file path if needed
            mcap_file, topic = data_source, "/ActuatorMotion"
            data = parse_mcap_file(mcap_file, topic)
        elif source_type == "raw":
            data = data_source
        else:
            raise ValueError(f"Unsupported source type: {source_type}")
        
        if not data:
            print("No data found in the source.")
            return None, None, None
        
        timestamps = []
        pixel_coords = []
        actual_coords = []
        
        # 正则表达式模式
        pattern = r"Pixel: \(([0-9.-]+), ([0-9.-]+)\), Actual: \(([0-9.-]+), ([0-9.-]+)\)"
        
        print(f"Processing {len(data)} messages...")
        valid_messages = 0
        
        for entry in data:
            # 将十六进制数据转换为字符串
            hex_data = entry['data']
            try:
                # 将十六进制转换为字节，然后解码为字符串
                byte_data = bytes.fromhex(hex_data)
                msg_string = byte_data.decode('utf-8', errors='ignore')
                
                # 使用正则表达式提取坐标
                match = re.search(pattern, msg_string)
                
                if match:
                    pixel_x = float(match.group(1))
                    pixel_y = float(match.group(2))
                    actual_x = float(match.group(3))
                    actual_y = float(match.group(4))
                    
                    timestamps.append(entry['time'])
                    pixel_coords.append((pixel_x, pixel_y))
                    actual_coords.append((actual_x, actual_y))
                    valid_messages += 1
                    
            except Exception as e:
                continue
        
        print(f"成功解析 {valid_messages} 个有效消息")
        
        # 转换为numpy数组便于处理
        timestamps = np.array(timestamps)
        if len(pixel_coords) > 0:
            pixel_coords = np.array(pixel_coords)
            actual_coords = np.array(actual_coords)
        
        return timestamps, pixel_coords, actual_coords
        
    except Exception as e:
        print(f"Error parsing data: {e}")
        return None, None, None


def load_data_from_mat(mat_file):
    """
    Load coordinate data from MATLAB .mat file.
    
    Args:
        mat_file (str): Path to the .mat file
    
    Returns:
        tuple: (timestamps, pixel_coords, actual_coords)
    """
    try:
        data = loadmat(mat_file)
        
        timestamps = data['timestamps'].flatten()
        pixel_coords = data['pixel_coords']
        actual_coords = data['actual_coords']
        
        print(f"Loaded data from {mat_file}")
        print(f"  - Timestamps: {len(timestamps)} points")
        print(f"  - Pixel coordinates: {pixel_coords.shape}")
        print(f"  - Actual coordinates: {actual_coords.shape}")
        
        return timestamps, pixel_coords, actual_coords
        
    except Exception as e:
        print(f"Error loading data from .mat file: {e}")
        return None, None, None
    
def save_extracted_data_to_mat(timestamps, pixel_coords, actual_coords, output_file):
    """
    Save extracted coordinate data to MATLAB .mat file.
    
    Args:
        timestamps: Array of timestamps
        pixel_coords: Array of pixel coordinates
        actual_coords: Array of actual coordinates
        output_file: Output .mat file path
    """
    try:
        # Prepare data structure for MATLAB
        data_dict = {
            'timestamps': timestamps,
            'pixel_coords': pixel_coords,
            'actual_coords': actual_coords,
            'info': {
                'description': 'Laser spot tracking data',
                'timestamp_unit': 'nanoseconds',
                'coordinate_unit': 'mm',
                'creation_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        }
        
        # Save to .mat file
        savemat(output_file, data_dict)
        print(f"Data saved to {output_file}")
        print(f"  - Timestamps: {len(timestamps)} points")
        print(f"  - Pixel coordinates: {pixel_coords.shape}")
        print(f"  - Actual coordinates: {actual_coords.shape}")
        
    except Exception as e:
        print(f"Error saving data to .mat file: {e}")

def slice_motion_data(timestamps, actual_coords, start_idx, end_idx):
    """
    截取timestamps, pixel_coords, actual_coords的某一段。

    Args:
        timestamps (np.ndarray): 时间戳数组
        pixel_coords (np.ndarray): 像素坐标数组
        actual_coords (np.ndarray): 实际坐标数组
        start_idx (int): 起始索引（包含）
        end_idx (int): 结束索引（不包含）

    Returns:
        tuple: (timestamps_slice, pixel_coords_slice, actual_coords_slice)
    """
    return (
        timestamps[start_idx:end_idx],
        actual_coords[start_idx:end_idx]
    )

def load_reference_trajectory(mat_file="figure8_extracted_data.mat"):
    try:
        if not os.path.exists(mat_file):
            print(f"Reference trajectory file not found: {mat_file}")
            return None, None
            
        data = loadmat(mat_file)
        
        # 打印所有可用的键，帮助调试
        print(f"Available keys in {mat_file}: {list(data.keys())}")
        
        # 检查期望的键是否存在
        if 'figure8_x' not in data or 'figure8_y' not in data:
            print(f"Expected keys 'figure8_x' and 'figure8_y' not found")
            print(f"Available keys: {[k for k in data.keys() if not k.startswith('__')]}")
            return None, None
        
        figure8_x = data['figure8_x'].flatten() - 4
        figure8_y = data['figure8_y'].flatten()
        
        print(f"Loaded reference trajectory from {mat_file}")
        print(f"  - Figure8 X: {len(figure8_x)} points")
        print(f"  - Figure8 Y: {len(figure8_y)} points")
        
        return figure8_x, figure8_y
        
    except Exception as e:
        print(f"Error loading reference trajectory from .mat file: {e}")
        return None, None

def process_multiple_files(file_configs):
    """
    Process multiple MCAP files and extract coordinate data.
    
    Args:
        file_configs (list): List of dictionaries with 'mcap_file', 'mat_file', 'label' keys
    
    Returns:
        dict: Dictionary containing data for each configuration
    """
    results = {}
    topic = "/ActuatorMotion"
    
    for config in file_configs:
        mcap_file = config['mcap_file']
        mat_file = config['mat_file']
        label = config['label']
        
        print(f"\n{'='*50}")
        print(f"Processing {label}")
        print(f"{'='*50}")
        
        if os.path.exists(mcap_file):
            print(f"Found MCAP file: {mcap_file}")
            
            # Parse data from MCAP file
            timestamps, pixel_coords, actual_coords = parse_actuator_motion_data(mcap_file, "mcap")
            
            # Save data to .mat file
            if timestamps is not None and len(pixel_coords) > 0:  
                print(f"Saving data to .mat file: {mat_file}")
                save_extracted_data_to_mat(timestamps, pixel_coords, actual_coords, mat_file)  

                results[label] = {
                    'timestamps': timestamps,
                    'pixel_coords': pixel_coords,
                    'actual_coords': actual_coords,
                    'color': config.get('color', 'blue'),
                    'marker': config.get('marker', 'o')
                }
                print(f"Successfully processed {len(pixel_coords)} data points for {label}")
            else:
                print(f"No valid data found in {mcap_file}")
        else:
            print(f"MCAP file not found: {mcap_file}")

    return results

def calculate_tracking_errors(actual_coords, time_seconds, ref_x, ref_y, ref_time):
    """
    Calculate ISE and RMSE between actual trajectory and reference trajectory.
    Also calculates errors for three time segments: 0-8s, 8-12s, 12-14s.
    
    Args:
        actual_coords: Actual coordinates array (N x 2)
        ref_x: Reference X coordinates
        ref_y: Reference Y coordinates  
        time_seconds: Time array for integration
        ref_time: Reference time array
        
    Returns:
        dict: Dictionary containing ISE and RMSE values
    """
    try:
        # Interpolate reference trajectory to match actual trajectory time points
        from scipy.interpolate import interp1d
        
        # Interpolate reference to match actual time points
        interp_ref_x = interp1d(ref_time, ref_x, kind='linear', fill_value='extrapolate')
        interp_ref_y = interp1d(ref_time, ref_y, kind='linear', fill_value='extrapolate')
        
        # Get reference coordinates at actual time points
        ref_x_interp = interp_ref_x(time_seconds)
        ref_y_interp = interp_ref_y(time_seconds)
        
        # Calculate errors
        error_x = actual_coords[:, 0] - ref_x_interp
        error_y = actual_coords[:, 1] - ref_y_interp
        
        # Calculate squared errors
        squared_error = error_x**2 + error_y**2
        
        # Overall statistics
        if len(time_seconds) > 1:
            dt = np.mean(np.diff(time_seconds))
            ise_total = np.trapz(squared_error, dx=dt)
        else:
            ise_total = 0
            
        rmse_total = np.sqrt(np.mean(squared_error))
        rmse_x_total = np.sqrt(np.mean(error_x**2))
        rmse_y_total = np.sqrt(np.mean(error_y**2))
        
        # Calculate segment-wise errors
        segments = {
            'Segment_1 (0-8s)': (0, 8),
            'Segment_2 (8-12s)': (8, 12),
            'Segment_3 (12-14s)': (12, 14)
        }
        
        segment_results = {}
        
        for seg_name, (start_time, end_time) in segments.items():
            # Find indices for this time segment
            mask = (time_seconds >= start_time) & (time_seconds < end_time)
            
            if np.any(mask):
                seg_time = time_seconds[mask]
                seg_error_x = error_x[mask]
                seg_error_y = error_y[mask]
                seg_squared_error = squared_error[mask]
                
                # Calculate ISE for this segment
                if len(seg_time) > 1:
                    seg_dt = np.mean(np.diff(seg_time))
                    seg_ise = np.trapz(seg_squared_error, dx=seg_dt)
                else:
                    seg_ise = 0
                
                # Calculate RMSE for this segment
                seg_rmse = np.sqrt(np.mean(seg_squared_error))
                seg_rmse_x = np.sqrt(np.mean(seg_error_x**2))
                seg_rmse_y = np.sqrt(np.mean(seg_error_y**2))
                
                segment_results[seg_name] = {
                    'ISE': seg_ise,
                    'RMSE': seg_rmse,
                    'RMSE_X': seg_rmse_x,
                    'RMSE_Y': seg_rmse_y,
                    'Max_Error': np.sqrt(np.max(seg_squared_error)),
                    'Mean_Error': np.sqrt(np.mean(seg_squared_error)),
                    'Data_Points': np.sum(mask)
                }
            else:
                segment_results[seg_name] = {
                    'ISE': 0,
                    'RMSE': 0,
                    'RMSE_X': 0,
                    'RMSE_Y': 0,
                    'Max_Error': 0,
                    'Mean_Error': 0,
                    'Data_Points': 0
                }
        
        # Combine all results
        results = {
            'Total': {
                'ISE': ise_total,
                'RMSE': rmse_total,
                'RMSE_X': rmse_x_total,
                'RMSE_Y': rmse_y_total,
                'Max_Error': np.sqrt(np.max(squared_error)),
                'Mean_Error': np.sqrt(np.mean(squared_error))
            },
            'Segments': segment_results
        }
        
        # Print results
        print(f"\n{'='*60}")
        print(f"TRACKING ERROR ANALYSIS")
        print(f"{'='*60}")
        
        # Print total results
        print(f"\nOVERALL RESULTS:")
        print(f"ISE (Integral of Squared Error): {results['Total']['ISE']:.4f} mm²·s")
        print(f"RMSE (Root Mean Squared Error): {results['Total']['RMSE']:.4f} mm")
        print(f"RMSE X: {results['Total']['RMSE_X']:.4f} mm")
        print(f"RMSE Y: {results['Total']['RMSE_Y']:.4f} mm")
        print(f"Maximum Error: {results['Total']['Max_Error']:.4f} mm")
        print(f"Mean Error: {results['Total']['Mean_Error']:.4f} mm")
        
        # Print segment results
        print(f"\nSEGMENT-WISE RESULTS:")
        for seg_name, seg_data in segment_results.items():
            print(f"\n{seg_name}:")
            print(f"  Data Points: {seg_data['Data_Points']}")
            print(f"  ISE: {seg_data['ISE']:.4f} mm²·s")
            print(f"  RMSE: {seg_data['RMSE']:.4f} mm")
            print(f"  RMSE X: {seg_data['RMSE_X']:.4f} mm")
            print(f"  RMSE Y: {seg_data['RMSE_Y']:.4f} mm")
            print(f"  Max Error: {seg_data['Max_Error']:.4f} mm")
        
        print(f"{'='*60}")
        
        return results
    
    except Exception as e:
        print(f"Error calculating tracking errors: {e}")
        return None

def visualize_multiple_trajectories(results_dict):
    """
    Visualize a single coordinate trajectory with time-based color mapping.
    (Assumes only one trajectory in results_dict)
    """
    if not results_dict:
        print("No valid coordinate data to visualize.")
        return

    # 只取第一条轨迹
    label, data = next(iter(results_dict.items()))
    timestamps = data['timestamps']
    actual_coords = data['actual_coords']

    if timestamps is None or len(actual_coords) == 0:
        print("No valid coordinate data to visualize.")
        return

    timestamps, actual_coords = slice_motion_data(timestamps, actual_coords, 6363, 6824)

    # Load reference trajectory
    ref_x, ref_y = load_reference_trajectory()

    # 设置字体为 Times New Roman
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 10

    # 转换时间戳为相对时间（秒）
    start_time = timestamps[0]
    relative_times = (timestamps - start_time) / 1e9  # 纳秒转秒

    # 设置字体大小
    title_size = 16
    label_size = 14
    tick_size = 12

    # 创建可视化
    fig, ax = plt.subplots(figsize=(4, 5))

    # Plot reference trajectory first (if available)
    if ref_x is not None and ref_y is not None:
        ax.plot(ref_x, ref_y, 'k--', linewidth=3, alpha=0.8, 
                label='Reference Trajectory', zorder=5)
        print(f"Added reference trajectory with {len(ref_x)} points")

    # 实际坐标轨迹
    ax.plot(actual_coords[:, 0], actual_coords[:, 1], 'g-', alpha=0.6, linewidth=1)
    scatter = ax.scatter(actual_coords[:, 0], actual_coords[:, 1], c=relative_times, cmap='plasma', s=12)
    cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.2, shrink=0.8)
    cbar.set_label('Time [s]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax.set_xlabel('X [mm]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax.set_ylabel('Y [mm]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    # ax.set_title('Laser Spot Trajectory', fontsize=title_size, fontfamily='Times New Roman', fontweight='bold', pad=15)
    ax.grid(True)
    # ax.axis('equal')

    # 设置xlim和ylim
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)

    # 设置坐标轴刻度标签的字体大小和字体
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')

    # 设置colorbar刻度标签字体
    cbar.ax.tick_params(labelsize=tick_size)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')

    plt.tight_layout()

    # 先保存图像，再显示
    current_date = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_filename = f"Koopman_actual_coordinate_trajectory_{current_date}.png"
    plt.savefig(plot_filename, dpi=600, bbox_inches='tight')
    print(f"可视化图像已保存为: {plot_filename}")

    # 然后显示图像
    plt.show()

def visualize_time_series(results_dict):
    """
    Visualize Y and Z coordinates over time for multiple trajectories.
    
    Args:
        results_dict (dict): Dictionary containing data for each frequency
    """
    if not results_dict:
        print("No valid coordinate data to visualize.")
        return
    
    # Set font to Times New Roman
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 10
    
    # Font sizes
    title_size = 14
    label_size = 12
    tick_size = 10
    legend_size = 9

    # Load reference trajectory
    ref_x, ref_y = load_reference_trajectory()
    if ref_x is not None and ref_y is not None:
        n = len(ref_x)
        # 三段时间长度
        t1 = np.linspace(0, 8, n)
        t2 = np.linspace(8, 12, n)
        t3 = np.linspace(12, 14, n)
        # 拼接
        ref_x_all = np.concatenate([ref_x, ref_x, ref_x])
        ref_y_all = np.concatenate([ref_y, ref_y, ref_y])
        ref_t_all = np.concatenate([t1, t2, t3])
        # 组合为trajectory
        ref_traj = np.stack([ref_x_all, ref_y_all, ref_t_all], axis=1)
        print(f"拼接后的参考轨迹 shape: {ref_traj.shape}")

    # Create subplots for Y and Z coordinates
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4))
    
    # Plot each trajectory
    label, data = next(iter(results_dict.items()))
    timestamps = data['timestamps']
    actual_coords = data['actual_coords']
    
    timestamps, actual_coords = slice_motion_data(timestamps, actual_coords, 6363, 6824)
    
    # Convert timestamps to relative time (seconds from start)
    if len(timestamps) > 0:
        # Convert nanoseconds to seconds
        time_seconds = (timestamps - timestamps[0]) / 1e9
        
        color = data['color']
        
        # Plot Y coordinate over time
        ax1.plot(time_seconds, actual_coords[:, 1], 
                color=color, linewidth=2, label=f'{label}', alpha=0.8)
        
        # Plot ref_y over time
        if ref_y is not None:
            ax1.plot(ref_t_all, ref_y_all, color='gray', linestyle='--', linewidth=2, label='Reference Y', alpha=0.8)

        # Plot Z coordinate over time (assuming Z is the third column if available)
        # If your data only has X,Y coordinates, we'll plot X instead
        if actual_coords.shape[1] > 2:
            z_coords = actual_coords[:, 2]
        else:
            # Use X coordinates as substitute for Z
            z_coords = actual_coords[:, 0]
            
        ax2.plot(time_seconds, z_coords, 
                color=color, linewidth=2, label=f'{label}', alpha=0.8)
        
        # Plot ref_x over time
        if ref_x is not None:
            ax2.plot(ref_t_all, ref_x_all, color='gray', linestyle='--', linewidth=2, label='Reference X', alpha=0.8)

    # Calculate tracking errors
    if ref_x is not None and ref_y is not None:
        error_results = calculate_tracking_errors(actual_coords, time_seconds, ref_x_all, ref_y_all, ref_t_all)

    # Configure Y coordinate plot
    ax1.set_ylabel('Y [mm]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    # ax1.legend(fontsize=legend_size, loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax1.set_ylim(-12, 12)

    # Configure Z coordinate plot (or X if Z not available)
    coord_label = 'Z [mm]' if actual_coords.shape[1] > 2 else 'X [mm]'
    ax2.set_xlabel('Time [s]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax2.set_ylabel(coord_label, fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=legend_size, loc='upper center', bbox_to_anchor=(0.5, -0.4), ncol=2)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.set_ylim(-12, 12)
    
    # Set font family for all tick labels
    for ax in [ax1, ax2]:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontfamily('Times New Roman')
    
    plt.tight_layout()
    
    # Save the plot
    current_date = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_filename = f"Koopman_time_series_coordinates_{current_date}.png"
    plt.savefig(plot_filename, dpi=600, bbox_inches='tight')
    print(f"Time series visualization saved as: {plot_filename}")
    
    # Show the plot
    plt.show()

def main():
    """
    Main function to parse and visualize multiple actuator motion data files.
    """
    # Configuration for multiple files
    # 选择plasma色系中的一个颜色
    plasma_color = plt.get_cmap('plasma')(0.6)[:3]  # 取RGB，不要alpha
    file_configs = [
        {
            'mcap_file': "rosbag_20250817_174352_0.mcap",
            'mat_file': "Koopman.mat",  # Changed from json_file to mat_file
            'label': "Koopman",
            'color': plasma_color,
            'marker': 'o'
        }
    ]
    
    print("Starting multi-file processing...")
    print("="*60)
    
    # Process all files
    results = process_multiple_files(file_configs)
    
    # Visualize combined trajectories
    if results:
        print(f"\nSuccessfully processed {len(results)} frequency configurations")
        print("Creating combined visualization...")
        visualize_multiple_trajectories(results)
        
        # Add time series visualization
        print("Creating time series visualization...")
        visualize_time_series(results)
    else:
        print("No valid data found in any of the files.")

if __name__ == "__main__":
    main()