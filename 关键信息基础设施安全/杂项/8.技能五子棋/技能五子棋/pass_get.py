import itertools
import zipfile
import string
import sys

# ================= 题目配置 =================
ZIP_FILENAME = "flag.zip"  # 你的压缩包文件名


# ===========================================

def brute_force_with_crc_check():
    try:
        zf = zipfile.ZipFile(ZIP_FILENAME)
        # 自动获取 ZIP 内的第一个文件名（哪怕是乱码）
        target_file = zf.namelist()[0]
        print(f"[*] 目标文件: {target_file}")
    except Exception as e:
        print(f"[!] ZIP 加载失败: {e}")
        return

    # === 1. 设置高概率的固定片段 (基于CTF常识) ===
    # fiv? -> 99% 是 'e' (five)
    char_fiv = ['e']

    # NN -> 99% 是 '19' (标准五子棋路数), 但为了保险我们保留所有符合数学逻辑的
    valid_numbers = ['19', '12', '16', '23', '26']

    # 5app? -> 99% 是 'e' (apple), 对应桥接字符也是 'e' (101+101=202)
    # 如果这一步跑不出来，可以将下面改为 string.ascii_lowercase 来跑全量
    bridge_pairs = [('e', 'e')]

    # === 2. 只有这些是真正的“未知数” ===
    # 密码结构: ..._e0{DDD}funa0{EE}
    # 我们需要爆破 DDD (3位) 和 EE (2位)
    # 范围: 小写字母 (CTF flag通常是小写英语句子)
    charset = string.ascii_lowercase

    total_combinations = len(valid_numbers) * (26 ** 3) * (26 ** 2)
    print(f"[*] 开始智能爆破...")
    print(f"[*] 尝试组合数: {total_combinations} (预计耗时: 几分钟)")

    count = 0

    # 开始循环
    for num in valid_numbers:
        for char_a in char_fiv:  # five
            for bridge_a, bridge_b in bridge_pairs:  # 5appe -> e0

                # 爆破中间 3 位 (DDD)
                for ddd in itertools.product(charset, repeat=3):
                    str_ddd = "".join(ddd)

                    # 爆破结尾 2 位 (EE)
                    for ee in itertools.product(charset, repeat=2):
                        str_ee = "".join(ee)

                        count += 1
                        if count % 10000 == 0:
                            print(f"\r[*] 已扫描: {count} 个密码...", end="")

                        # 拼凑完整密码
                        # 模板: skill_five_{num}_danc_sing_5appe_e0{DDD}funa0{EE}
                        password = (
                            f"skill_fiv{char_a}_"
                            f"{num}_"
                            f"danc_"
                            f"sing_"
                            f"5app{bridge_a}_"
                            f"{bridge_b}0{str_ddd}funa0{str_ee}"
                        )

                        # === 关键步骤：利用 CRC 校验 ===
                        try:
                            # 必须使用 .read() 读取整个文件流
                            # 只有密码完全正确，CRC 才会通过，否则抛出 BadZipFile
                            with zf.open(target_file, pwd=password.encode('utf-8')) as f:
                                f.read()

                            # 如果代码能走到这里，说明没有任何报错！
                            print(f"\n\n[★] 找到唯一真密码！CRC 校验通过！")
                            print(f"[★] Password: {password}")
                            print(f"[+] 语义分析: ... e0 {str_ddd} fun a0 {str_ee} ...")
                            return

                        except (zipfile.BadZipFile, RuntimeError):
                            # BadZipFile: 密码凑巧过了头校验，但 CRC 挂了 (碰撞)
                            # RuntimeError: 密码头校验就挂了 (错误密码)
                            continue
                        except Exception as e:
                            print(f"\n[!] 异常: {e}")


if __name__ == "__main__":
    brute_force_with_crc_check()