from mcap.reader import make_reader
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from datetime import datetime
import os

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
                    
                    if valid_messages <= 5:  # 只显示前5个坐标
                        print(f'提取的像素坐标: ({pixel_x}, {pixel_y})')
                        print(f'提取的实际坐标: ({actual_x}, {actual_y})')
                else:
                    # 只显示前几个无法解析的消息
                    if len(timestamps) < 3:
                        print(f'无法解析消息: {msg_string[:50]}...')
                    
            except Exception as e:
                if len(timestamps) < 3:
                    print(f"Error decoding hex data: {e}")
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

def visualize_coordinates(timestamps, pixel_coords, actual_coords):
    """
    Visualize actual coordinates trajectory only.
    
    Args:
        timestamps: Array of timestamps
        pixel_coords: Array of pixel coordinates (x, y)
        actual_coords: Array of actual coordinates (x, y)
    """
    if timestamps is None or len(actual_coords) == 0:
        print("No valid coordinate data to visualize.")
        return
    
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
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # 实际坐标轨迹
    ax.plot(actual_coords[:, 0], actual_coords[:, 1], 'g-', alpha=0.6, linewidth=1)
    scatter = ax.scatter(actual_coords[:, 0], actual_coords[:, 1], c=relative_times, cmap='plasma', s=10)
    cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.2, shrink=0.8)
    cbar.set_label('Time [s]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax.set_xlabel('X [mm]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax.set_ylabel('Y [mm]', fontsize=label_size, fontfamily='Times New Roman', fontweight='bold')
    ax.set_title('Laser Spot Trajectory', fontsize=title_size, fontfamily='Times New Roman', fontweight='bold', pad=15)
    ax.grid(True)
    ax.axis('equal')
    
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
    plot_filename = f"actual_coordinate_trajectory_{current_date}.png"
    plt.savefig(plot_filename, dpi=600, bbox_inches='tight')
    print(f"可视化图像已保存为: {plot_filename}")
    
    # 然后显示图像
    plt.show()
    
def save_extracted_data(data, output_file):
    """
    Save extracted data to JSON file.
    
    Args:
        data: Extracted data list
        output_file: Output file path
    """
    try:
        with open(output_file, "w") as out_file:
            json.dump(data, out_file, indent=4)
        print(f"Data saved to {output_file}")
    except Exception as e:
        print(f"Error saving data: {e}")

def main():
    """
    Main function to parse and visualize actuator motion data.
    支持从MCAP文件或JSON文件读取数据
    """
    # 配置文件路径
    mcap_file = "rosbag_20250709_201227-20250709T122130Z-1-001/rosbag_20250709_201227/rosbag_20250709_201227_0.mcap"
    json_file = "actuator_motion_data.json"
    topic = "/ActuatorMotion"
    
    # 检查文件存在性并选择数据源
    if os.path.exists(mcap_file):
        print(f"发现MCAP文件: {mcap_file}")
        print(f"开始从MCAP文件解析数据...")
        
        # 从MCAP文件解析数据
        timestamps, pixel_coords, actual_coords = parse_actuator_motion_data(mcap_file, "mcap")
        
        # 可选：保存提取的原始数据到JSON
        if timestamps is not None and len(pixel_coords) > 0:
            print("同时保存原始数据到JSON文件...")
            raw_data = parse_mcap_file(mcap_file, topic)
            save_extracted_data(raw_data, json_file)
        
    elif os.path.exists(json_file):
        print(f"发现JSON文件: {json_file}")
        print(f"开始从JSON文件解析数据...")
        
        # 从JSON文件解析数据
        timestamps, pixel_coords, actual_coords = parse_actuator_motion_data(json_file, "json")
        
    else:
        print(f"未找到数据文件:")
        print(f"  MCAP文件: {mcap_file}")
        print(f"  JSON文件: {json_file}")
        print("请确保至少有一个文件存在。")
        return
    
    # 可视化数据
    if timestamps is not None and len(pixel_coords) > 0:
        print(f"成功提取 {len(pixel_coords)} 个坐标数据点")
        visualize_coordinates(timestamps, pixel_coords, actual_coords)
    else:
        print("未找到有效的坐标数据")

if __name__ == "__main__":
    main()