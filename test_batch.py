import subprocess
import re
import argparse
import sys
import os

def run_spiral_command():
    """运行 ./spiral 命令并捕获输出"""
    cmd = ["./spiral", "12", "8", "1024", "--batch"]
    print(f"Running command: {' '.join(cmd)} ...")
    
    try:
        # capture_output=True 将标准输出捕获到 result.stdout
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print("Standard Error Output:", e.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: Executable './spiral' not found in the current directory.")
        sys.exit(1)

def parse_log_data(log_content):
    """解析日志内容，对重复项求和，提取单项和大小"""
    
    # 定义数据结构
    metrics = {
        "single_times": {}, # 存储只出现一次的时间
        "sum_times": {      # 存储需要累加的时间
            "time_expand_query": 0.0,
            "time_reorient_ct": 0.0,
            "time_mul_db_with_query": 0.0,
            "time_fold_further_dim": 0.0,
            "time_convert_mod": 0.0
        },
        "sizes": {}         # 存储大小信息
    }

    # 正则表达式
    # 匹配: time_xxx uses 123.45 ms
    pattern_uses = re.compile(r'(time_\w+)\s+uses\s+([\d\.]+)\s+ms')
    # 匹配: time_xxx is 123.45 (例如 time_batch_query)
    pattern_is = re.compile(r'(time_\w+)\s+is\s+([\d\.]+)')
    # 匹配: num_bytes_B: 12345
    pattern_size_colon = re.compile(r'(num_bytes_B):\s+(\d+)')
    # 匹配: size_of_xxx is 12345
    pattern_size_is = re.compile(r'(size_of_\w+)\s+is\s+(\d+)')

    lines = log_content.splitlines()
    
    for line in lines:
        line = line.strip()
        
        # 1. 解析 "uses X ms"
        match = pattern_uses.search(line)
        if match:
            key, val = match.group(1), float(match.group(2))
            if key in metrics["sum_times"]:
                metrics["sum_times"][key] += val
            else:
                metrics["single_times"][key] = val
            continue

        # 2. 解析 "is X" (时间)
        match = pattern_is.search(line)
        if match:
            key, val = match.group(1), float(match.group(2))
            metrics["single_times"][key] = val
            continue

        # 3. 解析大小 (num_bytes_B)
        match = pattern_size_colon.search(line)
        if match:
            metrics["sizes"][match.group(1)] = int(match.group(2))
            continue
            
        # 4. 解析大小 (size_of_...)
        match = pattern_size_is.search(line)
        if match:
            metrics["sizes"][match.group(1)] = int(match.group(2))
            continue

    return metrics

def save_report(metrics, output_file):
    """格式化数据并写入文件"""
    s_times = metrics["single_times"]
    sum_times = metrics["sum_times"]
    sizes = metrics["sizes"]

    # 准备输出内容
    output_lines = []
    output_lines.append("=" * 50)
    output_lines.append(f"Spiral Batch Experiment Report")
    output_lines.append("=" * 50)
    
    output_lines.append("\n[1] Key Time Metrics (Summed for Batch)")
    output_lines.append("-" * 40)
    # 打印累加的时间
    for k, v in sum_times.items():
        output_lines.append(f"{k:<30} : {v:.4f} ms")
    
    output_lines.append("\n[2] Single Step Times")
    output_lines.append("-" * 40)
    # 打印一些关键的单次时间 (你可以根据需要添加更多)
    keys_to_show = [
        "time_set_constants", "time_encode_db", 
        "time_batch_query", "time_decode_respond"
    ]
    for k in keys_to_show:
        if k in s_times:
            output_lines.append(f"{k:<30} : {s_times[k]:.4f} ms")
        else:
            output_lines.append(f"{k:<30} : N/A")

    output_lines.append("\n[3] Data Sizes")
    output_lines.append("-" * 40)
    for k, v in sizes.items():
        # 自动转换可读性更强的单位 (MB/KB)
        readable = f"{v} bytes"
        if v > 1024*1024:
            readable += f" (~{v/1024/1024:.2f} MB)"
        elif v > 1024:
            readable += f" (~{v/1024:.2f} KB)"
        output_lines.append(f"{k:<30} : {readable}")

    # 写入文件
    try:
        with open(output_file, 'w') as f:
            f.write("\n".join(output_lines))
        print(f"\nSuccess! Report saved to: {output_file}")
    except IOError as e:
        print(f"Error writing to file: {e}")

if __name__ == "__main__":
    # 命令行参数处理
    parser = argparse.ArgumentParser(description="Run Spiral experiment and parse logs.")
    parser.add_argument("-o", "--output", type=str, default="experiment_result.txt", 
                        help="Path to save the result report (default: experiment_result.txt)")
    
    args = parser.parse_args()

    # 1. 运行命令
    log_output = run_spiral_command()
    
    # 2. 解析数据
    parsed_metrics = parse_log_data(log_output)
    
    # 3. 保存报告
    save_report(parsed_metrics, args.output)