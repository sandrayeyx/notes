import base64

print("="*60)
print("欢迎来到技能五子棋学校!")
print("="*60)
print("传统的五子棋😯就是把五个子连成一条线👓")
print("好无趣😫好无聊😩")
print("而技能五子棋🤓☝就是在传统的五子棋🧐加入技能🤭")
print("好好玩🥰💥💥💥要爆了💥💥💥")
print()
print("校长技能五和学生子棋正在守护着一个password...")
print("只有掌握了真正的技能五子棋精髓,才能获得这个password!")
print("="*60)
print()

password = input('请输入你的password: ')

if len(password) != 42:
    print("❌ 技能五: 飞沙走石! 你的password长度不对,被我扔进什刹海了!")
    print(f"   (当前长度: {len(password)}, 要求: 42)")
    exit(1)

if password.count("_") != 6:
    print("❌ 技能五: 力拔山兮! 我把棋盘都摔了,下划线数量不对!")
    print(f"   (当前: {password.count('_')}个, 要求: 6个)")
    exit(1)

if password[0:6] != "skill_":
    print("❌ 子棋: 两极反转! 翻转棋盘,发现你不懂技能!")
    print(f"   (password[0:6] = '{password[0:6]}', 要求: 'skill_')")
    exit(1)

parts = password.split("_")
if list(map(len, parts)) != [5, 4, 2, 4, 4, 5, 12]:
    print("❌ 技能五: 换个顺序! 分段长度不对!")
    print(f"   (当前分段长度: {list(map(len, parts))}, 要求: [5,4,2,4,4,5,12])")
    exit(1)

if base64.b64encode(password[6:10].encode('utf-8')).decode()[:4] != 'Zml2':
    print("❌ 王教练: See You Again! 手刀劈下,编码不对!")
    print(f"   (password[6:10] base64编码前4位 = '{base64.b64encode(password[6:10].encode()).decode()[:4]}', 要求: 'Zml2')")
    exit(1)

if not password[11:13].isdigit():
    print("❌ 王教练: 胜天半子! 这里应该是数字!")
    exit(1)
    
number = int(password[11:13])
if number <= 10 or number >= 30 or pow(number, 2, 7) != 4:
    print("❌ 王教练: 你的数字修为不够!")
    print(f"   (password[11:13] = {number}, 要求: 10<n<30 且 n²%7=4)")
    exit(1)

if "sing" not in password or "danc" not in password:
    print("❌ 技能五: 你有多久没有在下五子棋的时候又唱又跳了?")
    print("   (必须包含 'sing' 和 'danc')")
    exit(1)

if password[19:23].replace("i", "1").replace("n", "0").replace("g", "9") != "s109":
    print("❌ 王教练: 不要把我徒弟当臭狗一样玩耍!")
    print(f"   (password[19:23] 替换后 = '{password[19:23].replace('i','1').replace('n','0').replace('g','9')}', 要求: 's109')")
    exit(1)

if password[14:18].swapcase() != "DANC":
    print("❌ 子棋: 大小写有问题!")
    print(f"   (password[14:18].swapcase() = '{password[14:18].swapcase()}', 要求: 'DANC')")
    exit(1)

if ord(password[28]) + ord(password[30]) != 202:
    print("❌ 子棋: ASCII值相加不对!")
    print(f"   (ord(password[28]) + ord(password[30]) = {ord(password[28])+ord(password[30])}, 要求: 202)")
    exit(1)

num_str = password[24] + password[31] + password[39]
if not num_str.isdigit() or int(num_str) * 37 != 18500:
    print("❌ 技能五: 计算有误!")
    print(f"   (password[24]+password[31]+password[39] = '{num_str}', {num_str}×37 应该 = 18500)")
    exit(1)

new = ""
for i in password[35:39]:
    new += chr(ord(i) + 1)
if new != "gvob":
    print("❌ 王教练: 字符位移不对!")
    print(f"   (password[35:39] 每个字符+1 = '{new}', 要求: 'gvob')")
    exit(1)

if password[24:28] != "5app":
    print("❌ 技能五: 五子棋的真谛是快乐! 但你缺少关键元素!")
    print(f"   (password[24:28] = '{password[24:28]}', 要求: '5app')")
    exit(1)

if "5app" in password and "fun" in password:
    print()
    print("🎉"*20)
    print("✨ 恭喜你! ✨")
    print("技能五: 太好了! 你领悟了五子棋的真谛!")
    print("子棋: 老师! 老师! 他做到了!")
    print("王教练: 不错! 你已经掌握了技能五子棋的精髓!")
    print()
    print("🎵 现在,让我们一起唱起来，跳起来! 🎵")
    print("技~能~五~子~棋~💃🕺")
    print("传统的CTF就是猜密码，好无趣好无聊~")
    print("而技能CTF就是在传统的CTF加入技能，好好玩~")
    print("💥💥💥 要爆了! 💥💥💥")
    print()
    print(f"你的password是: {password}")
    print("🎉"*20)
else:
    print("❌ 还差一点点! 五子棋的真谛是什么?")
    exit(1)