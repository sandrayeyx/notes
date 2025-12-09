def solve():
    # 从提供的汇编数据中提取的 byte_403060 数组 (0x403060 - 0x403082)
    enc_bytes = [
        0x3E, 0x0D, 0xED, 0xBE, 0x4A, 0x8A, 0x7D, 0xBC, 0x7C, 0xFC,
        0x2E, 0x2A, 0x79, 0x9D, 0x6A, 0x1A, 0xCC, 0x3D, 0x4A, 0xF8,
        0x3C, 0x79, 0x69, 0x39, 0xD9, 0xDD, 0x9D, 0xA9, 0x69, 0x4C,
        0x8C, 0xDD, 0x59, 0xE9, 0xD7
    ]

    # 复制一份以进行操作
    flag = list(enc_bytes)

    # 步骤 1: 逆向 ROR 4 (循环右移 4 位)
    # 在 8 位字节中，右移 4 位相当于交换高 4 位和低 4 位。
    # 逆操作就是再做一次 ROR 4。
    for i in range(len(flag)):
        flag[i] = ((flag[i] >> 4) | (flag[i] << 4)) & 0xFF

    # 步骤 2: 逆向累加 (Forward Accumulation)
    # 原逻辑: v10[i] += v10[i+1]
    # 逆逻辑: v10[i] -= v10[i+1]
    # 注意：必须从倒数第二个数开始，往前遍历
    for i in range(len(flag) - 2, -1, -1):
        flag[i] = (flag[i] - flag[i+1]) & 0xFF

    # 转换为字符串
    result = "".join([chr(x) for x in flag])
    print(f"Flag: {result}")

if __name__ == "__main__":
    solve()