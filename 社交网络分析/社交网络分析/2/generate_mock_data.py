import json
import random
import math
import numpy as np
from collections import defaultdict

# --- 配置参数 ---
N_USERS = 500  # 模拟用户数量 (论文为 ~500k)
N_CITIES = 10  # 模拟城市数量
AVG_NON_GEO_FRIENDS = 2.5  # 论文估算的非地理朋友 [cite: 175]
AVG_GEO_FRIENDS = 5.5  # 论文估算的地理朋友 [cite: 176, 184]

print(f"开始生成 {N_USERS} 个用户和 {N_CITIES} 个城市的模拟数据...")

# --- 1. 生成 us_cities.json ---

cities = {}
city_ids = [f"city_{i}" for i in range(N_CITIES)]

for city_id in city_ids:
    cities[city_id] = {
        "latitude": random.uniform(24.0, 49.0),  # 模拟美国纬度
        "longitude": random.uniform(-125.0, -66.0),  # 模拟美国经度
        "population": random.randint(1000, 5000000)
    }

# 将城市数据写入文件
try:
    with open('us_cities.json', 'w') as f:
        json.dump(cities, f, indent=2)
    print("成功创建 'us_cities.json'")
except Exception as e:
    print(f"创建 'us_cities.json' 时出错: {e}")

# --- 2. 生成 livejournal_data.json ---

users = {}
friendships = defaultdict(set)
user_locations = {}  # 缓存: user_id -> (lon, lat)
user_ids = [f"user_{i}" for i in range(N_USERS)]

# (A) 分配用户到城市
for user_id in user_ids:
    # 随机分配一个城市
    assigned_city_id = random.choice(city_ids)
    users[user_id] = {
        "location": assigned_city_id,
        "username": f"user_{user_id}"
    }
    user_locations[user_id] = (cities[assigned_city_id]['longitude'], cities[assigned_city_id]['latitude'])


def haversine(lon1, lat1, lon2, lat2):
    """ 简单的球面距离估算 (用于模拟) """
    R = 6371  # 地球半径
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


# (B) 为友谊链接预计算距离和排名 (这模拟了论文的核心思想 [cite: 243])
print("正在预计算用户排名 (模拟论文的排名机制)...")
all_ranks = {}
for u_id in user_ids:
    u_lon, u_lat = user_locations[u_id]
    distances = []
    for v_id in user_ids:
        if u_id == v_id:
            continue
        v_lon, v_lat = user_locations[v_id]
        dist = haversine(u_lon, u_lat, v_lon, v_lat)
        distances.append((v_id, dist))

    # 按距离排序
    distances.sort(key=lambda x: x[1])

    # 计算排名 [cite: 243]
    ranks_for_u = {}
    current_rank = 0
    last_dist = -1
    for v_id, dist in distances:
        if dist > last_dist:
            current_rank = len(ranks_for_u)  # 排名 = 比我近的人数
        ranks_for_u[v_id] = current_rank
        last_dist = dist
    all_ranks[u_id] = ranks_for_u

# (C) 生成友谊链接
print("正在生成友谊链接...")
for u_id in user_ids:

    # 1. 添加非地理朋友 (随机) [cite: 168]
    num_non_geo = np.random.poisson(AVG_NON_GEO_FRIENDS)
    for _ in range(num_non_geo):
        v_id = random.choice(user_ids)
        if v_id != u_id:
            friendships[u_id].add(v_id)

    # 2. 添加地理朋友 (基于排名) [cite: 244, 246]
    num_geo = np.random.poisson(AVG_GEO_FRIENDS)
    ranks = all_ranks[u_id]

    valid_targets = []
    weights = []

    for v_id, rank in ranks.items():
        if rank < 1:  # 忽略 rank 0
            continue
        valid_targets.append(v_id)
        weights.append(1.0 / rank)  # 论文的核心: Pr ∝ 1/rank [cite: 246]

    if not valid_targets or not weights:
        continue

    # 归一化概率
    total_weight = sum(weights)
    probabilities = [w / total_weight for w in weights]

    # 抽样
    if num_geo > 0 and len(valid_targets) > num_geo:
        geo_friends = np.random.choice(
            valid_targets,
            size=num_geo,
            replace=False,
            p=probabilities
        )
        for v_id in geo_friends:
            friendships[u_id].add(v_id)

# (D) 格式化友谊为列表 (JSON需要)
final_friendships = {}
for u_id, friend_set in friendships.items():
    final_friendships[u_id] = list(friend_set)

# (E) 组装并写入文件
livejournal_data = {
    "users": users,
    "friendships": final_friendships
}

try:
    with open('livejournal_data.json', 'w') as f:
        json.dump(livejournal_data, f)  # 为节省空间，不使用 indent
    print("成功创建 'livejournal_data.json'")
except Exception as e:
    print(f"创建 'livejournal_data.json' 时出错: {e}")

print("\n--- 模拟数据生成完毕 ---")
print("您现在可以运行 'replicate_livejournal_analysis.py' 脚本了。")