import requests
import re

# 1. 设置 URL
base_url = "http://172.17.0.15:14909"  # 请替换为实际 URL
file_url = base_url + "/0"

print(f"[*] 正在下载文件: {file_url}")

try:
    r = requests.get(file_url)

    if r.status_code == 200:
        content = r.content  # 获取二进制内容
        print(f"[+] 下载成功，文件大小: {len(content)} 字节")
        print("[*] 文件头 (Hex):", content[:10].hex())

        # 2. 尝试在文件中搜索 flag 关键字
        # Flag 通常以 flag{ 或 <?php 开头
        print("\n[*] 正在搜索 Flag ...")

        # 尝试转换成文本查找
        try:
            text_content = content.decode('utf-8', errors='ignore')

            # 搜索 flag{...}
            flags = re.findall(r'flag\{.*?\}', text_content)
            if flags:
                print(f"\n[SUCCESS] 找到 Flag 了！！！")
                for f in flags:
                    print(f">>> {f}")
            else:
                print("[-] 未找到 'flag{' 格式的字符串。")

            # 搜索 PHP 源码
            if "<?php" in text_content:
                print("\n[+] 发现 PHP 源码片段:")
                start_index = text_content.find("<?php")
                print(text_content[start_index:start_index + 200])

            # 如果上面都没找到，直接打印最后 500 个字符
            print("\n[*] 文件末尾 500 字符 (Flag 可能在这里):")
            print("-" * 50)
            print(text_content[-500:])
            print("-" * 50)

        except Exception as e:
            print(f"解析错误: {e}")

    else:
        print(f"[-] 下载失败，状态码: {r.status_code}")
        print("请先运行之前的攻击脚本生成文件 '0'")

except Exception as e:
    print(f"Error: {e}")