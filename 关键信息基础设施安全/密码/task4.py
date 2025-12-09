import socket
from Crypto.Cipher import AES

#可以用，用调试，直接运行不行

def bxor(b1, b2):	# 对二进制数据进行异或
    result = b""
    for b_1, b_2 in zip(b1, b2):
        result += bytes([b_1 ^ b_2])
    return result



# Server details
SERVER_IP = '172.17.0.15'
SERVER_PORT = 16095

def read_line(sock):
    body = b""
    while True:
        ch = sock.recv(1)
        if ch == b"\n":
            break
        body += ch
    return body


def get_token(sock):
    # Read initial server message
    print(sock.recv(1024).decode())
    
    while True:
        # Read the menu
        menu = sock.recv(1024).decode()
        print(menu)
        
        if "Enter your choice:" in menu:
            break
    
    # Choose to create a HUSTCTFer account
    sock.send(b"1\n")
    
    # Read the token message
    response = sock.recv(1024).decode()
    print(response)
    
    if ":" in response:
        tokens = response.split(":")
        if len(tokens) >= 2:
            user_token = tokens[1].strip()
            print(f"User token: {user_token}")
        else:
            print("Invalid token format received.")
            return None
    else:
        print("Invalid token line format received.")
        return None
    
    return user_token

def main():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.connect((SERVER_IP, SERVER_PORT))
        
        # Get the token from the server
        user_token = get_token(sock)
        
        if user_token is None:
            print("Failed to retrieve valid user token. Exiting.")
            return
        

        authcode = user_token	# 这里是得到的数据
        user_key = bytes.fromhex(authcode)[:16]	# 拿出iv
        code = bytes.fromhex(authcode)[16:]		# 结果密文
        user_key = bxor(user_key, b'AdminAdmin!_____')	# 异或目标明文
        user_key = bxor(user_key, b'HUSTCTFer!______')	# 异或原明文
        authcode = user_key.hex() + code.hex()	# 最终结果:经过构造的authcode
        print(authcode)
        # Known plaintext/ciphertext pair
        
        # Choose to login
        sock.send(b"3\n")
        
        # Read prompt for token
        print(sock.recv(1024).decode())

        sock.send(authcode.encode() + b"\n")
        
        # Read response
        response = sock.recv(1024).decode()
        print(response)
        
        if "Hello Admin! Here is your FLAG:" in response:
            flag = sock.recv(1024).decode()
            print(f"FLAG: {flag}")
        else:
            print("Failed to get the flag.")

if __name__ == "__main__":
    main()

