import networkx as nx
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- 1. 数据集生成 (保持原题逻辑) ---
random.seed(42)
np.random.seed(42)

nodes = list(range(50))
edges = []
while len(edges) < 150:
    u = random.choice(nodes)
    v = random.choice(nodes)
    if u != v and (u, v) not in edges and (v, u) not in edges:
        weight = random.randint(1, 10)
        edges.append((u, v, weight))

G = nx.Graph()
G.add_nodes_from(nodes)
G.add_weighted_edges_from(edges)

# --- 2. 预处理：为路径类指标计算“距离” ---
# 互动频次越高，距离越近。因此定义 distance = 1 / weight
for u, v, data in G.edges(data=True):
    data['distance'] = 1.0 / data['weight']

# --- 3. 计算中心性指标 ---

# (1) 度中心性 (标准化)
degree_cen = nx.degree_centrality(G)

# (2) 加权度 (通常指节点强度 Strength，不自动标准化，这里手动归一化以便比较)
weighted_degree = dict(G.degree(weight='weight'))
max_wd = max(weighted_degree.values())
weighted_degree_norm = {k: v / max_wd for k, v in weighted_degree.items()}

# (3) 介数中心性 (无权)
betweenness_cen = nx.betweenness_centrality(G, weight=None)

# (4) 加权介数中心性 (使用倒数权重作为距离)
weighted_betweenness_cen = nx.betweenness_centrality(G, weight='distance')

# (5) 接近中心性 (使用倒数权重作为距离)
closeness_cen = nx.closeness_centrality(G, distance='distance')

# (6) 特征向量中心性 (使用原始权重作为强度)
eigenvector_cen = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)

# --- 4. 整理数据与筛选 Top 5 ---
metrics_data = {
    "Degree": degree_cen,
    "Weighted_Degree": weighted_degree_norm,
    "Betweenness": betweenness_cen,
    "Weighted_Betweenness": weighted_betweenness_cen,
    "Closeness_Weighted": closeness_cen,
    "Eigenvector_Weighted": eigenvector_cen
}

# 创建 DataFrame 并保留4位小数
df = pd.DataFrame(metrics_data)
df = df.round(4)

# 打印各指标 Top 5
print("=== 各中心性指标 Top 5 节点 ===")
top_nodes_set = set() # 用于后续可视化标记
for col in df.columns:
    top5 = df.nlargest(5, col)
    print(f"\n[{col}] Top 5:")
    print(top5)
    # 记录排名前3的节点用于可视化
    top3 = df.nlargest(3, col).index.tolist()
    top_nodes_set.update(top3)

# --- 5. 可视化核心节点 ---
plt.figure(figsize=(12, 10))

# 布局算法
pos = nx.spring_layout(G, seed=42, k=0.3)

# 基础节点绘制
nx.draw_networkx_nodes(G, pos, node_color='lightgray', node_size=100, alpha=0.6)

# 突出显示核心节点 (至少在两种指标中排名前3的节点，或者简单地突出所有指标的前3)
# 这里我们突出显示“加权介数”和“加权特征向量”的前3名
wb_top3 = df.nlargest(3, 'Weighted_Betweenness').index.tolist()
ev_top3 = df.nlargest(3, 'Eigenvector_Weighted').index.tolist()

# 绘制不同类型的核心节点
# 红色：信息枢纽 (加权介数高)
nx.draw_networkx_nodes(G, pos, nodelist=wb_top3, node_color='red', node_size=300, label='High Weighted Betweenness')
# 蓝色：意见领袖 (加权特征向量高)
nx.draw_networkx_nodes(G, pos, nodelist=ev_top3, node_color='blue', node_size=300, label='High Weighted Eigenvector')
# 紫色：重叠节点 (如果是同一个节点)
overlap = list(set(wb_top3) & set(ev_top3))
if overlap:
    nx.draw_networkx_nodes(G, pos, nodelist=overlap, node_color='purple', node_size=400, label='Overlapped Core')

# 绘制边 (根据权重调整粗细)
weights = [G[u][v]['weight'] * 0.3 for u, v in G.edges()]
nx.draw_networkx_edges(G, pos, width=weights, alpha=0.4)

# 标签
nx.draw_networkx_labels(G, pos, font_size=8)

plt.title("Campus Club Interaction Network: Core Nodes Analysis")
plt.legend()
plt.axis('off')
plt.show()