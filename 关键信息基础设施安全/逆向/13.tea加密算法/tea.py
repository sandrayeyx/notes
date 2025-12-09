import struct


def decrypt(v, k):
    v0, v1 = v[0], v[1]
    sum = 0xC6EF3720  # delta * 32
    delta = 0x9e3779b9
    k0, k1, k2, k3 = k[0], k[1], k[2], k[3]

    for _ in range(32):
        v1 -= ((v0 << 4) + k2) ^ (v0 + sum) ^ ((v0 >> 5) + k3)
        v1 &= 0xFFFFFFFF
        v0 -= ((v1 << 4) + k0) ^ (v1 + sum) ^ ((v1 >> 5) + k1)
        v0 &= 0xFFFFFFFF
        sum -= delta
        sum &= 0xFFFFFFFF

    return [v0, v1]


if __name__ == "__main__":
    # 提取的密钥
    key = [0x11223344, 0x55667788, 0x99aabbcc, 0xddeeff00]

    # 修正后的密文 (3组，每组2个32位整数)
    cipher_blocks = [
        # 修复了不可见字符: 06 D6 1D 82 -> 0x821DD606
        [0x821DD606, 0x56EF64F8],
        # 修复了末尾字符: 64 28 78 81 -> 0x81782864
        [0x4B3F4472, 0x81782864],
        # 第三组无误
        [0x7C9337E4, 0xE7B6B67A]
    ]

    flag = b""
    for block in cipher_blocks:
        decrypted = decrypt(block, key)
        # 将解密后的整数转换为字节串 (小端序)
        flag += struct.pack("<II", decrypted[0], decrypted[1])

    print("Flag: " + flag.decode('utf-8'))