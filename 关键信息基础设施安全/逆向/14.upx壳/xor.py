def solve_from_code():
    # 1. 获取 C 代码中的硬编码密文 (Str2)
    # 注意：反汇编代码中的字符串包含转义字符
    # "dnpiz#g\\X|3DMa'a2PM{\\`\"n"
    # 在 Python 中表示如下：
    encrypted_str = "dnpiz#g\\X|3DMa'a2PM{\\`\"n"

    # 将字符串转换为可变的字节列表 (ASCII 值)
    flag_chars = list(map(ord, encrypted_str))
    length = len(flag_chars)  # 应该是 24

    print(f"密文长度: {length}")
    print(f"密文内容: {encrypted_str}")

    # 2. 逆向操作 - 第二轮循环 (Loop 2)
    # C代码: Str[j] ^= v5[j % 3]; 其中 v5 = {3, 22, 7}
    v5 = [3, 22, 7]
    for j in range(length):
        flag_chars[j] ^= v5[j % 3]

    # 3. 逆向操作 - 第一轮循环 (Loop 1)
    # C代码根据 i % 3 的值异或不同的数
    for i in range(length):
        mod = i % 3
        if mod == 0:
            # else: Str[i] ^= 0x11u;
            flag_chars[i] ^= 0x11
        elif mod == 1:
            # if ( i % 3 == 1 ) Str[i] ^= 0x15u;
            flag_chars[i] ^= 0x15
        else:
            # elseStr[i] ^= 0x14u;
            flag_chars[i] ^= 0x14

    # 4. 组合结果
    flag = "".join(map(chr, flag_chars))
    print("-" * 30)
    print(f"解密出的 Flag: {flag}")


if __name__ == "__main__":
    solve_from_code()