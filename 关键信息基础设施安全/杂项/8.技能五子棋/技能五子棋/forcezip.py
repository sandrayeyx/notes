import zipfile


def try_zip(dict_file, zip_path):
    print(f"开始尝试解压 {zip_path} ...")
    with open(dict_file, 'r') as f:
        passwords = f.readlines()

    success = False
    with zipfile.ZipFile(zip_path) as zf:
        for p in passwords:
            pwd = p.strip()
            try:
                # 尝试解压第一个文件来验证密码
                zf.extractall(pwd=pwd.encode('utf-8'))
                print(f"\n✨ 成功找到密码: {pwd}")
                success = True
                break
            except (RuntimeError, zipfile.BadZipFile):
                # 密码错误
                pass
            except Exception as e:
                pass

    if not success:
        print("\n❌ 字典中的密码都不对。")

# 使用方法
# try_zip("pass_dict.txt", "secret.zip")