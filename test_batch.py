import subprocess
import re
import argparse
import sys
import os

def run_spiral_command():
    """运行 ./spiral 命令，实时打印输出并捕获用于解析"""
    cmd = ["./spiral", "12", "8", "1024", "--batch"]
    print(f"Running command: {' '.join(cmd)} ...")
    print("-" * 60) # 分割线，方便查看开始

    captured_output = []
    
    try:
        # 使用 Popen 建立管道，bufsize=1 表示行缓冲，text=True 表示处理文本
        # stdout=subprocess.PIPE 用于捕获输出
        # stderr=subprocess.STDOUT 将错误输出也合并到标准输出中
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            bufsize=1,
            encoding='utf-8' # 显式指定编码，防止解码错误
        )

        # 实时逐行读取
        for line in process.stdout:
            # 1. 实时打印到终端 (这就让你能看到进度了)
            sys.stdout.write(line)
            sys.stdout.flush() 
            
            # 2. 同时存入列表供后续解析
            captured_output.append(line)

        # 等待进程结束
        process.wait()

        if process.returncode != 0:
            print(f"\n[Error] Command failed with return code {process.returncode}")
            sys.exit(process.returncode)

        return "".join(captured_output)

    except FileNotFoundError:
        print("Error: Executable './spiral' not found in the current directory.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def parse_log_data(log_content):
    """解析日志内容，对重复项求和，提取单项和大小"""
    
    metrics = {
        "single_times": {},
        "sum_times": {
            "time_expand_query": 0.0,
            "time_reorient_ct": 0.0,
            "time_mul_db_with_query": 0.0,
            "time_fold_further_dim": 0.0,
            "time_convert_mod": 0.0
        },
        "sizes": {}
    }

    # 正则表达式
    pattern_uses = re.compile(r'(time_\w+)\s+uses\s+([\d\.]+)\s+ms')
    pattern_is = re.compile(r'(time_\w+)\s+is\s+([\d\.]+)')
    pattern_size_colon = re.compile(r'(num_bytes_B):\s+(\d+)')
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
    for k, v in sum_times.items():
        output_lines.append(f"{k:<30} : {v:.4f} ms")
    
    output_lines.append("\n[2] Single Step Times")
    output_lines.append("-" * 40)
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
        readable = f"{v} bytes"
        if v > 1024*1024:
            readable += f" (~{v/1024/1024:.2f} MB)"
        elif v > 1024:
            readable += f" (~{v/1024:.2f} KB)"
        output_lines.append(f"{k:<30} : {readable}")

    # 自动创建父目录 (防止 benchmark_results 文件夹不存在报错)
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"\n[Info] Created directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}")
            return

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
    parser.add_argument("-o", "--output", type=str, default="./benchmark_results/batch_log.txt", 
                        help="Path to save the result report (default: ./benchmark_results/batch_log.txt)")
    
    args = parser.parse_args()

    # 1. 运行命令 (现在会实时打印了)
    log_output = run_spiral_command()
    
    # 2. 解析数据
    parsed_metrics = parse_log_data(log_output)
    
    # 3. 保存报告
    save_report(parsed_metrics, args.output)