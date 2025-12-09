def bxor(b1, b2):	# 对二进制数据进行异或
    result = b""
    for b_1, b_2 in zip(b1, b2):
        result += bytes([b_1 ^ b_2])
    return result

authcode = "fbc78b91a1755467467391d0162d1ad27d712502e9c88f7dd3aebc5ce37bf75b"	# 这里是得到的数据
#f7dadea72293748d1aba42e01bdfb35c7144ce8236c3193f5af3506d19d6182a
user_key = bytes.fromhex(authcode)[:16]	# 拿出iv
code = bytes.fromhex(authcode)[16:]		# 结果密文
user_key = bxor(user_key, b'AdminAdmin!_____')	# 异或目标明文
user_key = bxor(user_key, b'HUSTCTFer!______')	# 异或原明文
authcode = user_key.hex() + code.hex()	# 最终结果:经过构造的authcode
print(authcode)