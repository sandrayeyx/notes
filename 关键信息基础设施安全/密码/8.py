from pwn import *
from Crypto.Util.number import isPrime, long_to_bytes, getPrime
from sympy.ntheory import discrete_log
import sys

# 题目提供的 IP 和端口
HOST = '172.17.0.15'
PORT = 16147

def get_smooth_prime(bits=1024):
    """
    构造一个弱素数 P，使得 P-1 是光滑的（由小素数相乘得到）。
    这使得离散对数问题 (DLP) 变得很容易解决。
    """
    print(f"[*] 正在生成 {bits} 位的弱素数 (Smooth Prime)...")
    
    # 1. 累乘小素数直到达到目标位数附近
    res = 2
    cur = 3
    while res.bit_length() < bits - 16: # 预留一些空间给 k
        res *= cur
        # 简单的寻找下一个素数的方法
        cur += 2
        while not isPrime(cur):
            cur += 2
            
    # 2. 寻找倍数 k，使得 p = res * k + 1 是素数
    # 这样 p-1 = res * k，大部分因子都是我们刚刚累乘的小素数
    k = 1
    while True:
        p = res * k + 1
        if p.bit_length() >= bits and isPrime(p):
            print(f"[+] 找到弱素数 P! (k={k})")
            return p
        k += 1

def solve():
    # 1. 生成弱素数
    p = get_smooth_prime(1024)
    
    # 2. 连接服务器
    try:
        r = remote(HOST, PORT)
    except:
        print(f"[-] 无法连接到 {HOST}:{PORT}")
        return

    # 3. 发送 P
    r.recvuntil(b'P = ')
    r.sendline(str(p).encode())
    
    # 4. 接收参数
    # 格式:
    # Alice公钥: ...
    # Bob公钥: ...
    # 密文: ...
    
    try:
        r.recvuntil(b': ')
        alice_pub = int(r.recvline().strip())
        print(f"[*] 收到 Alice 公钥: {alice_pub}")
        
        r.recvuntil(b': ')
        bob_pub = int(r.recvline().strip())
        print(f"[*] 收到 Bob 公钥: {bob_pub}")
        
        r.recvuntil(b': ')
        cipher_txt = int(r.recvline().strip())
        print(f"[*] 收到密文: {cipher_txt}")
    except ValueError:
        print("[-] 接收数据解析失败，可能 P 不符合服务器要求")
        r.close()
        return

    # 5. 攻击：使用 Pohlig-Hellman 算法求离散对数
    # 我们要求解 7^a = alice_pub (mod p) 中的 a
    print("[*] 正在利用 Pohlig-Hellman 算法计算私钥 a ...")
    
    # sympy.ntheory.discrete_log 会自动检测 p-1 的因子并使用 Pohlig-Hellman
    # 因为我们构造的 p-1 非常光滑，这步会非常快
    pri_key_a = discrete_log(p, alice_pub, 7)
    
    print(f"[+] 破解成功! 私钥 a = {pri_key_a}")

    # 6. 计算共享密钥 Superkey
    # superkey = B^a % p
    shared_secret = pow(bob_pub, pri_key_a, p)
    print(f"[*] 共享密钥 S = {shared_secret}")

    # 7. 解密
    # c = m * S % p  ==>  m = c * inverse(S, p) % p
    # pow(x, -1, p) 是求模逆元
    m = (cipher_txt * pow(shared_secret, -1, p)) % p
    
    flag = long_to_bytes(m)
    print(f"\n[+] FLAG: {flag.decode(errors='ignore')}")
    
    r.close()

if __name__ == '__main__':
    solve()