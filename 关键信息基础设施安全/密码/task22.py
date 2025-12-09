from Crypto.Util.number import *
from gmpy2 import *
from pwn import *
from time import sleep
from ecdsa import ecdsa as ec
import hashlib

n = 6277101735386680763835789423176059013767194773182842284081
r = []
s = []
e = []

host = '172.17.0.15'
port = '16033'

def get_data(stop=True):
    io = remote(host, port)
    io.recvuntil(">")
    io.sendline("1")
    data = io.recvline(False).strip().decode().split(",")
    io.close()
    r.append(int(data[0]))
    s.append(int(data[1]))
    e.append(int(data[2]))
    if stop:
        sleep(1)

# Getting data twice
for i in range(2):
    get_data()

def get_secret(s, e):
    return (s[0] * e[1] - s[1] * e[0]) * inverse(s[1] - s[0], n) * inverse(r[1], n) % n

secret = get_secret(s, e)
print(secret)

g = ec.generator_192
PUBKEY = ec.Public_key(g, g * secret)
PRIVKEY = ec.Private_key(PUBKEY, secret)

io = remote(host, port)
# context.log_level = 'debug'

io.recvuntil(">")
io.sendline("2")

io.recvuntil('Get signature for md5("')
to_check = io.recv(17)
print(to_check)

to_check = int(hashlib.md5(to_check).hexdigest(), 16)
signature = PRIVKEY.sign(to_check, 114514)

io.recvline(False)
io.sendline(str(signature.r) + "," + str(signature.s))

io.interactive()
