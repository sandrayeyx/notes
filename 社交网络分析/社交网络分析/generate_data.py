import pandas as pd

import random

# 生成500个用户，1000条关注关系

users = list(range(1, 501))

edges = []

for _ in range(1000):

    source = random.choice(users)

    target = random.choice(users)

    if source != target:   # 排除自关注

          edges.append((source, target))

# 确保存在3-5个密集社区（手动增加部分连接）

for community in [range(1, 101), range(101, 201), range(201, 301)]:

    for _ in range(100):  # 每个社区内增加100条内部连接

        source = random.choice(community)

        target = random.choice(community)

        if source != target:

            edges.append((source, target))

# 去重并保存

df = pd.DataFrame(edges, columns=['source', 'target']).drop_duplicates()

df = df.sample(1000).reset_index(drop=True)  # 确保最终1000条

df.to_csv('social_network.csv', index=False)

