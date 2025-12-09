import zlib
import struct

# 题目中的关键数据
# IHDR chunk type (4 bytes)
chunk_type = b'IHDR'
# Width: 0cd8 (3288)
width = struct.pack('>I', 0x0cd8)
# Params (Bit depth, Color type, etc.): 08 02 00 00 00
params = b'\x08\x02\x00\x00\x00'
# 目标 CRC (题目截图里的 39 ca 35 4d)
expected_crc = 0x39ca354d

# 爆破高度 (从 1000 开始尝试到 4000)
for height_int in range(1000, 4000):
    height = struct.pack('>I', height_int)
    # 计算当前组合的 CRC
    data = chunk_type + width + height + params
    calculated_crc = zlib.crc32(data) & 0xffffffff

    if calculated_crc == expected_crc:
        print(f"找到真实高度了! Hex: {hex(height_int)}")
        print(f"请将 '00 00 03 e8' 修改为 -> '{height.hex()}'")
        break
else:
    print("范围内未找到，请尝试扩大爆破范围。")