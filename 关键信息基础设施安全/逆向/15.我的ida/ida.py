# 逆向脚本
# 目标格式: vmc{...}

def solve():
    # 原始题目(V&N 2020)的加密数据 (byte_403060)
    # 如果你的题目是修改版，请在 IDA 中跳转到 0x403060 复制 35 个字节替换这里
    # 注意：这里使用常见的 Poison 题目数据进行演示
    enc = [
        0x3E, 0x0D, 0x0ED, 0x0BE, 0x4A, 0x8A, 0x7D, 0x0BC, 0x7C, 0x0FC, 0x2E, 0x2A, 0x79, 0x9D, 0x6A, 0x1A,
        0x0CC, 0x3D, 0x4A, 0x0F8, 0x3C, 0x79, 0x69, 0x39, 0x0D9, 0x0DD, 0x9D, 0x0A9, 0x69, 0x4C, 0x8C, 0x0DD,
        0x59, 0x0E9, 0x0D7
    ]

    # 确保数据是列表格式以便修改
    flag = list(enc)

    # --- 第一步逆向：逆向 ROR 4 ---
    # 在 8bit 中，ROR 4 和 ROL 4 效果一样，都是交换高低 4 位
    # 例如: 0xC5 (1100 0101) -> 0x5C (0101 1100)
    for i in range(len(flag)):
        flag[i] = ((flag[i] >> 4) | (flag[i] << 4)) & 0xFF

    # --- 第二步逆向：逆向累加 ---
    # 正向: v10[i] += v10[i+1] (从 0 到 33)
    # 逆向: v10[i] -= v10[i+1] (从 33 到 0)
    # 必须从后往前减
    for i in range(len(flag) - 2, -1, -1):
        flag[i] = (flag[i] - flag[i+1]) & 0xFF

    # --- 输出结果 ---
    result = "".join([chr(x) for x in flag])
    print(f"Flag: {result}")

if __name__ == "__main__":
    solve()