from sympy import Matrix
import sys

# 密文数据
ciphertext = [
    62, 117, 16, 
    108, 57, 10, 
    112, 73, 44, 
    48, 4, 94, 
    74, 0, 105, 
    98, 3, 12, 
    21, 56, 100, 
    31, 8, 73, 
    120, 107, 10, 
    70, 25, 102, 
    122, 20, 80, 
    84, 4, 3, 
    62, 82, 126
]

MOD = 127

# 辅助函数：判断字符是否合法
def is_valid(c):
    if c == 0x20: return True
    try:
        ch = chr(c)
        return ch.isalnum() or ch in '{}_'
    except:
        return False

# 辅助函数：解密
def try_decrypt(key_matrix_inv, cipher_vecs):
    decrypted_str = ""
    try:
        for i in range(len(cipher_vecs)):
            # 向量乘矩阵: [1x3] * [3x3]
            v = Matrix(1, 3, cipher_vecs[i])
            res = (v * key_matrix_inv) % MOD
            
            # 检查每个字符是否合法
            for val in res:
                if not is_valid(int(val)):
                    return None
            
            decrypted_str += "".join([chr(int(val)) for val in res])
        return decrypted_str
    except:
        return None

def main():
    # 将密文预处理为分块列表
    cipher_vecs = []
    for i in range(len(ciphertext) // 3):
        cipher_vecs.append(ciphertext[3*i : 3*i+3])

    print("[*] 正在运行解密脚本 (Python版)...")

    # 构造已知的部分矩阵 (密文部分)
    # Block 0, Block 1, Block Last
    me_rows = [
        ciphertext[0:3],       # Block 0
        ciphertext[3:6],       # Block 1
        ciphertext[-3:]        # Block Last
    ]
    me = Matrix(me_rows)

    # ---------------------------------------------------------
    # 策略 3: 填充 2 个空格 (最常见情况) -> ['}', 0x20, 0x20]
    # ---------------------------------------------------------
    print("[-] 尝试策略 3 (填充2个空格)...")
    for guess1 in range(0x20, 127):
        if not is_valid(guess1): continue
        for guess2 in range(0x20, 127):
            if not is_valid(guess2): continue
            
            # 构造明文矩阵 P (mc)
            # Row 1: vmc
            # Row 2: { + guess1 + guess2
            # Row 3: } + space + space
            mc = Matrix([
                [ord('v'), ord('m'), ord('c')],
                [ord('{'), guess1, guess2],
                [ord('}'), 0x20, 0x20]
            ])

            # 计算 Key
            # P * K = C  => K = P^-1 * C
            if mc.det() % MOD == 0: continue # 不可逆则跳过
            
            try:
                p_inv = mc.inv_mod(MOD)
                key = (p_inv * me) % MOD
                
                # 计算解密矩阵 (Key^-1)
                if key.det() % MOD == 0: continue
                key_inv = key.inv_mod(MOD)

                # 尝试解密
                flag = try_decrypt(key_inv, cipher_vecs)
                if flag:
                    print(f"\n[+] 成功解密 (策略3): {flag}")
                    return
            except:
                continue

    # ---------------------------------------------------------
    # 策略 1: 填充 3 个空格 -> [0x20, 0x20, 0x20]
    # ---------------------------------------------------------
    print("[-] 尝试策略 1 (填充3个空格)...")
    for guess1 in range(0x20, 127):
        for guess2 in range(0x20, 127):
            mc = Matrix([
                [ord('v'), ord('m'), ord('c')],
                [ord('{'), guess1, guess2],
                [0x20, 0x20, 0x20]
            ])
            if mc.det() % MOD == 0: continue
            try:
                p_inv = mc.inv_mod(MOD)
                key = (p_inv * me) % MOD
                if key.det() % MOD == 0: continue
                key_inv = key.inv_mod(MOD)
                flag = try_decrypt(key_inv, cipher_vecs)
                if flag:
                    print(f"\n[+] 成功解密 (策略1): {flag}")
                    return
            except:
                continue

    # ---------------------------------------------------------
    # 策略 2: 填充 1 个空格 -> [guess3, '}', 0x20]
    # ---------------------------------------------------------
    print("[-] 尝试策略 2 (填充1个空格)...")
    # 这个计算量较大，如果前面没出结果再跑这个
    for guess1 in range(0x20, 127):
        if not is_valid(guess1): continue
        for guess2 in range(0x20, 127):
            if not is_valid(guess2): continue
            for guess3 in range(0x20, 127):
                if not is_valid(guess3): continue
                
                mc = Matrix([
                    [ord('v'), ord('m'), ord('c')],
                    [ord('{'), guess1, guess2],
                    [guess3, ord('}'), 0x20]
                ])
                if mc.det() % MOD == 0: continue
                try:
                    p_inv = mc.inv_mod(MOD)
                    key = (p_inv * me) % MOD
                    if key.det() % MOD == 0: continue
                    key_inv = key.inv_mod(MOD)
                    flag = try_decrypt(key_inv, cipher_vecs)
                    if flag:
                        print(f"\n[+] 成功解密 (策略2): {flag}")
                        return
                except:
                    continue

    print("[-] 未找到 Flag")

if __name__ == "__main__":
    main()