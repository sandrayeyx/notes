# 1. 从 .data 段提取的加密数据 (dest 数组)
dest_raw = [
    0x42, 0x3F, 0x53, 0x24, 0x2C, 0x66, 0x5D, 0x5F, 0x7F, 0x0A, 0x27,
    0x39, 0x5F, 0x1D, 0x0B, 0x0F, 0x2E, 0x00, 0x6B, 0x2B, 0x5B, 0x19,
    0x5C, 0x41
]

# 2. 虚拟机初始化的参数 (从 .rdata 的 MOV 指令提取)
seed = 0x47  # R[3]
multiplier = 23333  # R[4]
increment = 19260817  # R[5]
modulus = 127  # R[6]

flag = ""

# 3. 模拟 VM 生成逻辑并异或还原 Flag
print("正在计算 Flag...")

for i in range(len(dest_raw)):
    # 对应 VM 中的 MUL, ADD, MOD 指令
    generated_val = (seed * multiplier + increment) % modulus

    # 对应 VM 中的 MOV R[3], Result (更新种子)
    seed = generated_val

    # Flag[i] = dest[i] ^ generated_val
    # 因为 main 函数里是 dest[i] ^= Str[i]，如果输入正确，
    # 结果应该等于 generated_val。
    char_code = dest_raw[i] ^ generated_val
    flag += chr(char_code)

print(f"Flag 结果: {flag}")