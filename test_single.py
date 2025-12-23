import subprocess
import re
import csv
import time
from collections import defaultdict
import os

# ================= 配置区域 =================

# 1. 你要执行的命令
COMMAND = ["./spiral", "12", "8", "45","--single"]

# 2. 运行次数
RUN_TIMES = 10 

# 3. 【新功能】文件存放目录 (可以是相对路径 "./results"，也可以是绝对路径 "/home/user/data")
OUTPUT_DIR = "/app/benchmark_results" 

# 4. 文件名
LOG_FILENAME = "benchmark_raw_log.txt"   # 日志文件名
CSV_FILENAME = "benchmark_data.csv"      # 表格文件名

# ===========================================

def parse_output(output_text):
    """ 解析输出文本，提取 metrics """
    metrics = {}
    
    # 1. 匹配 time_xxxx uses 123.456 ms
    time_uses_pattern = re.findall(r'(time_\w+)\s+uses\s+([\d\.]+)\s+ms', output_text)
    for key, val in time_uses_pattern:
        metrics[key] = float(val)

    # 2. 匹配 time_xxxx is 123.456
    time_is_pattern = re.findall(r'(time_\w+)\s+is\s+([\d\.]+)', output_text)
    for key, val in time_is_pattern:
        metrics[key] = float(val)

    # 3. 匹配 size_of_xxxx is 12345
    size_pattern = re.findall(r'(size_of_\w+)\s+is\s+(\d+)', output_text)
    for key, val in size_pattern:
        metrics[key] = int(val)

    return metrics

def main():
    # --- 路径处理逻辑 ---
    # 自动创建目录（如果不存在）
    if not os.path.exists(OUTPUT_DIR):
        try:
            os.makedirs(OUTPUT_DIR)
            print(f"已创建输出目录: {OUTPUT_DIR}")
        except OSError as e:
            print(f"创建目录失败: {e}")
            return

    # 拼接完整路径
    full_log_path = os.path.join(OUTPUT_DIR, LOG_FILENAME)
    full_csv_path = os.path.join(OUTPUT_DIR, CSV_FILENAME)
    # ------------------

    # 1. 初始化日志文件
    with open(full_log_path, "w", encoding="utf-8") as f:
        f.write(f"Benchmark Started at {time.ctime()}\n")
        f.write(f"Command: {' '.join(COMMAND)}\n")
        f.write("="*50 + "\n\n")

    print(f"开始测试。")
    print(f"1. 原始日志将保存在: {full_log_path}")
    print(f"2. 数据表格将保存在: {full_csv_path}")
    print(f"总共运行次数: {RUN_TIMES}\n")

    all_runs_data = [] 
    all_keys = set()   

    # 2. 循环运行
    for i in range(1, RUN_TIMES + 1):
        print(f"正在运行第 {i}/{RUN_TIMES} 次...", end="", flush=True)
        
        try:
            # 执行命令
            result = subprocess.run(COMMAND, capture_output=True, text=True, check=True)
            
            # 写入原始日志
            with open(full_log_path, "a", encoding="utf-8") as f:
                f.write(f"--- Run #{i} ---\n")
                f.write(result.stdout)
                f.write("\n" + "-"*30 + "\n\n")

            # 解析数据
            current_metrics = parse_output(result.stdout)
            
            if current_metrics:
                current_metrics["run_id"] = i 
                all_runs_data.append(current_metrics)
                all_keys.update(current_metrics.keys())
                print(" 完成")
            else:
                print(" [警告: 无数据]")

        except subprocess.CalledProcessError as e:
            print(f" 失败! (查看日志)")
            with open(full_log_path, "a", encoding="utf-8") as f:
                f.write(f"--- Run #{i} FAILED ---\n")
                f.write(e.stderr)
        except Exception as e:
            print(f" 未知错误: {e}")
            return

    # 3. 计算平均值并构建 Average 行
    if all_runs_data:
        data_keys = [k for k in all_keys if k != "run_id"]
        
        avg_row = {"run_id": "Average"}
        
        for key in data_keys:
            values = [d[key] for d in all_runs_data if key in d]
            if values:
                avg_val = sum(values) / len(values)
                avg_row[key] = round(avg_val, 4)
            else:
                avg_row[key] = 0

        # 4. 导出 CSV
        sorted_keys = sorted(data_keys)
        time_keys = [k for k in sorted_keys if "time" in k]
        size_keys = [k for k in sorted_keys if "size" in k]
        other_keys = [k for k in sorted_keys if k not in time_keys and k not in size_keys]
        
        fieldnames = ["run_id"] + other_keys + time_keys + size_keys
        
        with open(full_csv_path, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_runs_data)
            writer.writerow(avg_row)
        
        print(f"\n成功! CSV 数据已保存至: {full_csv_path}")
        
        print("\n" + "="*50)
        print(f"{'METRIC (Average)':<30} | {'VALUE':<15}")
        print("-" * 50)
        for key in fieldnames:
            if key != "run_id" and key in avg_row:
                print(f"{key:<30} | {avg_row[key]}")
        print("="*50)

if __name__ == "__main__":
    main()