def solve():
    # [cite_start]从文件中提取的加密字符串片段 [cite: 21]
    # 注意：片段3中的反斜杠和引号需要转义
    part1 = "dnpiz#g\\"
    part2 = "X|3DMa'a"
    part3 = "2PM{\\`\"n"

    # 拼接完整的密文
    ciphertext = part1 + part2 + part3

    # 已知的 Flag 前缀
    known_prefix = "vmc"

    # 1. 计算密钥 (Key)
    key = []
    print("正在计算密钥...")
    for i in range(len(known_prefix)):
        # Key = 密文 ^ 明文
        k = ord(ciphertext[i]) ^ ord(known_prefix[i])
        key.append(k)
        print(f"Key[{i}] = {k}")

    print(f"推导出的完整密钥序列: {key}")
    print("-" * 30)

    # 2. 解密全文字符串
    flag = ""
    for i in range(len(ciphertext)):
        # 循环使用密钥 [18, 3, 19]
        k = key[i % len(key)]
        decoded_char = chr(ord(ciphertext[i]) ^ k)
        flag += decoded_char

    print(f"解密结果 (Flag): {flag}")


if __name__ == "__main__":
    solve()