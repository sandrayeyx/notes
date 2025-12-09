import itertools
import zipfile
import string
import concurrent.futures
import sys

# ================= 配置 =================
ZIP_FILENAME = "flag.zip"


# =======================================

def check_password_chunk(passwords):
    """
    工作线程函数：负责验证一批密码
    """
    try:
        # 每个线程必须有自己独立的 ZipFile 句柄，否则会冲突
        local_zf = zipfile.ZipFile(ZIP_FILENAME)
        target_file = local_zf.namelist()[0]

        for pwd in passwords:
            try:
                # 优化：只读1个字节做快速筛选
                # 这比读取全文快几千倍，足以筛掉 99.6% 的错误密码
                with local_zf.open(target_file, pwd=pwd.encode('utf-8')) as f:
                    f.read(1)

                    # 如果没报错，再读全文校验 CRC（防止碰撞）
                with local_zf.open(target_file, pwd=pwd.encode('utf-8')) as f:
                    f.read()

                return pwd  # 找到了！
            except (RuntimeError, zipfile.BadZipFile):
                continue
            except Exception:
                continue
        return None
    except:
        return None


def main():
    print(f"[*] 正在生成智能字典...")

    # === 1. 核心常量 ===
    prefix = "skill_five_19_danc_sing_5appe_e0"
    mid_fix = "funa0"

    # === 2. 智能候选集 (优先跑这些) ===
    # 既然是 fun a... 很有可能是 have fun and / good fun and
    likely_ddd = ['hav', 'has', 'had', 'joy', 'god', 'bad', 'the', 'and', 'for', 'get']
    likely_ee = ['nd', 'th', 'st', 'er', 'in', 'on', 'at']

    priority_passwords = []

    # 生成高优先级的组合 (几百个)
    for d in likely_ddd:
        for e in likely_ee:
            priority_passwords.append(f"{prefix}{d}{mid_fix}{e}")

    # === 3. 备用全量集 (如果优先级没跑出来) ===
    # 为了节省内存，我们使用生成器，或者分批生成
    # 这里为了演示，我们只生成剩余的小写字母组合

    # 这里的技巧是：我们假设 DDD 还是比较像单词的，EE 也是
    # 如果要跑 26^5 全量，列表会太大。
    # 我们直接把生成的任务分发给线程。

    # 4. 执行高优先级检查
    print(f"[*] 阶段1: 检查 {len(priority_passwords)} 个高频语义组合...")
    if result := check_password_chunk(priority_passwords):
        print(f"\n[★] 密码破解成功: {result}")
        return

    # 5. 执行全量暴力破解 (多线程)
    print(f"[*] 阶段2: 启动多线程暴力破解 (请耐心等待)...")

    charset = string.ascii_lowercase

    # 这里的策略：只爆破 EE (2位)，遍历 DDD (3位)
    # 每次把 DDD 变动作为一批任务分发

    # 创建线程池 (根据你的CPU核心数自动调整)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []

        # 遍历 DDD (26*26*26 = 17576 种情况)
        for ddd_tuple in itertools.product(charset, repeat=3):
            ddd = "".join(ddd_tuple)

            # 为每个 DDD 生成 26*26 = 676 个 EE 的组合，作为一个任务包
            batch_passwords = []
            for ee_tuple in itertools.product(charset, repeat=2):
                ee = "".join(ee_tuple)
                batch_passwords.append(f"{prefix}{ddd}{mid_fix}{ee}")

            # 提交给线程池
            future = executor.submit(check_password_chunk, batch_passwords)
            futures.append(future)

            # 内存保护：如果堆积太多任务，稍微等一下（防止内存爆炸）
            if len(futures) > 1000:
                done, not_done = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                # 检查已完成的任务是否有结果
                for f in done:
                    if res := f.result():
                        print(f"\n\n[★] 密码破解成功: {res}")
                        # 强制终止所有线程
                        executor.shutdown(wait=False)
                        sys.exit(0)
                futures = list(not_done)
                print(f"\r[*] 正在推进... 当前前缀: {ddd}", end="")

        # 处理剩余结果
        for f in concurrent.futures.as_completed(futures):
            if res := f.result():
                print(f"\n\n[★] 密码破解成功: {res}")
                return


if __name__ == "__main__":
    main()