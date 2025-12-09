import math
from sympy import sieve
from sympy.ntheory import discrete_log
from Crypto.Util.number import long_to_bytes, inverse

# --- Data from smooth_rsa.py ---
n = 66122138029123637485622099063384580083505557745377597177203540940107575483455058441534745742726077100259626684486447094393469907413324927330740820933519482530159725896029132352481400945391730343310634496954420537559010455362501044047450178627331989654988592965545106987007529060454334215987903766593
c1 = 21643050659076843907390847868727807681135941718467556563674681275224858864609049163572978045117239167231777317871629250999940530618738378951567389401980792651444944928056052429547365695100287961220002659806950315830891983573749530855181159697406975406101106163304152386749257224249365054344832956654
c2 = 48845665578329487665088087274045255799985193039351789200728024209343318696909901579635279695976739506024148740932384412362259951474389120780787194715296649493703067927280758279464021100626229299819246164493699108809798819622321771099243521972536365215150318715191817229298659462419417494815097421651
e = 65537


def solve():
    print("[-] Attempting to factor n using Pollard's p-1...")

    # 1. Factor N
    # The smoothness bound B is roughly 2^20 (size of blocks used in generation)
    B = 2 ** 20
    a = 2

    # Iterate through primes. We compute a = a^(LCM(1..B)) mod n
    # We approximate the LCM exponent by ensuring every prime p <= B
    # contributes p^k where p^k is large enough to cover the order.
    # Since p-1 is constructed from 20-bit primes and small bases, this works perfectly.
    primes = sieve.primerange(2, B + 1000)

    for i, p_val in enumerate(primes):
        # Calculate exponent k to cover high multiplicities of small factors
        # For small bases like 2 or 3, k needs to be high (e.g. ~512 bits)
        # For large bases near 2^20, k=1 is sufficient.
        k = int(math.log(2 ** 512, p_val))

        # Modular exponentiation update
        a = pow(a, pow(p_val, k), n)

        # Check GCD periodically to speed up
        if i % 1000 == 0:
            g = math.gcd(a - 1, n)
            if 1 < g < n:
                p = g
                break
    else:
        # Final check if loop finishes
        g = math.gcd(a - 1, n)
        if 1 < g < n:
            p = g
        else:
            print("[!] Failed to factor n.")
            return

    q = n // p
    print(f"[+] Factored n!\n    p = {p}\n    q = {q}")

    # 2. Decrypt Part 1 (RSA)
    print("[-] Decrypting Part 1 (Standard RSA)...")
    phi = (p - 1) * (q - 1)
    d = inverse(e, phi)
    m1 = pow(c1, d, n)
    flag_part1 = long_to_bytes(m1)
    print(f"[+] Part 1: {flag_part1}")

    # 3. Decrypt Part 2 (DLP)
    print("[-] Decrypting Part 2 (Discrete Log)...")
    # Equation: c2 = e^m2 mod n
    # Reduce to mod p: c2 = e^m2 mod p
    # Solves m2 = log_e(c2) in Z_p
    # SymPy's discrete_log uses Pohlig-Hellman automatically for smooth orders
    try:
        m2 = discrete_log(p, c2 % p, e)
        flag_part2 = long_to_bytes(m2)
        print(f"[+] Part 2: {flag_part2}")
    except Exception as err:
        print(f"[!] DLP Failed: {err}")
        return

    # 4. Combine
    full_flag = flag_part1 + flag_part2
    print(f"\n[+] Full Flag: {full_flag.decode()}")


if __name__ == "__main__":
    solve()