import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# --- 配置（确保中文字体显示正常）---
#  plt.rcParams['font.sans-serif'] = ['SimHei']  # 'SimHei' 是黑体
# plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# （在某些环境中，可能需要指定字体路径）
# my_font = fm.FontProperties(fname="/path/to/your/font.ttf") # 替换为你系统中的中文字体路径
# 如果找不到字体，请取消注释上面两行，并提供你系统上的中文字体文件路径（如 .ttf 文件）
# 在 Colab 或标准环境中，可以尝试以下：
try:
    plt.rcParams['font.sans-serif'] = ['Source Han Sans CN', 'SimHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
except Exception as e:
    print(f"中文字体设置警告: {e}。图表中的中文可能显示为方框。")
    print("请尝试安装 'Source Han Sans CN' 或 'SimHei' 字体，或在代码中指定本地字体路径。")

# --------------------------
# 任务 2.1: 数据清洗
# --------------------------
print("--- 任务 2.1: 数据清洗开始 ---")

# (1) 加载数据
try:
    members_df = pd.read_csv("社团成员信息.csv")
    online_df = pd.read_csv("线上互动数据.csv")
    offline_df = pd.read_csv("线下活动数据.csv")
except FileNotFoundError:
    print("错误：未找到 CSV 数据文件。请先运行数据生成代码。")
    # 停止执行，因为后续步骤依赖这些文件
    exit()

print(f"原始数据量：线上 {len(online_df)}, 线下 {len(offline_df)}")

# (2) 清洗线上数据 (online_df)
# 标准 1：处理缺失值（receiver_id 缺失的数据无法构成关系）
online_cleaned = online_df.dropna(subset=['receiver_id']).copy()
print(f"清洗后 (去除缺失 receiver_id): {len(online_cleaned)} 条")

# 标准 2：处理重复数据（完全一致的记录）
# 模拟数据中包含时间戳，完全重复的概率较低，但按题意执行
online_cleaned = online_cleaned.drop_duplicates()
print(f"清洗后 (去除重复记录): {len(online_cleaned)} 条")

# (3) 清洗线下数据 (offline_df)
# 标准 3：处理异常值（duration 过长）
# 我们将时长超过 60 分钟的视为异常值（根据模拟代码，异常值为 120-180）
# 处理方式：采用“封顶”（Capping）的方式，将其设置为 60 分钟，保留该次互动
anomaly_threshold = 60
offline_cleaned = offline_df.copy()
offline_cleaned['duration_clipped'] = offline_cleaned['duration'].clip(upper=anomaly_threshold)
print(
    f"线下数据异常值处理：{len(offline_df[offline_df['duration'] > anomaly_threshold])} 条记录的时长被封顶至 {anomaly_threshold} 分钟。")

# --------------------------
# 任务 2.1 (续): 数据融合与网络构建
# --------------------------
print("\n--- 数据融合与网络构建 ---")

# (4) 统一权重标准
# 目标：将线上 (weight 1-5) 和线下 (duration_clipped) 转换为统一的互动强度 (strength)
# 策略：
# 1. 线上：strength = weight
# 2. 线下：将 duration (0-60分钟) 映射到 1-5 的区间
#    我们使用 pd.cut 进行分箱：(0, 10] -> 1, (10, 20] -> 2, ..., (50, 60] -> 6
#    为了匹配 1-5，我们调整分箱
bins = [0, 10, 20, 30, 45, 61]  # 0-10, 11-20, 21-30, 31-45, 46-60+
labels = [1, 2, 3, 4, 5]
offline_cleaned['strength'] = pd.cut(offline_cleaned['duration_clipped'], bins=bins, labels=labels, right=True,
                                     include_lowest=True).astype(int)

# (5) 构建边列表 (Edge List)
# 线上数据（有向）
# 修正：原始列名为 'weight'，在此处选取并重命名为 'strength' 以便统一
online_edges = online_cleaned[['sender_id', 'receiver_id', 'weight']].rename(columns={'weight': 'strength'})

# 线下数据（无向，需转换为有向）
# 一次线下互动 (A, B) 视为 A->B 和 B->A 两次强度相同的互动
offline_edges_1 = offline_cleaned[['person_a', 'person_b', 'strength']].rename(
    columns={'person_a': 'sender_id', 'person_b': 'receiver_id'})
offline_edges_2 = offline_cleaned[['person_b', 'person_a', 'strength']].rename(
    columns={'person_b': 'sender_id', 'person_a': 'receiver_id'})
offline_edges = pd.concat([offline_edges_1, offline_edges_2])

# (6) 合并并聚合
all_edges = pd.concat([online_edges, offline_edges])

# 聚合：将所有 A->B 的互动强度相加
# 这将作为网络中边的最终权重 (final_weight)
edge_list_df = all_edges.groupby(['sender_id', 'receiver_id'])['strength'].sum().reset_index()
edge_list_df = edge_list_df.rename(columns={'strength': 'final_weight'})

print(f"网络构建完成：共 {len(members_df)} 个节点, {len(edge_list_df)} 条聚合后的有向边。")

# (7) 构建 NetworkX 图 (DiGraph)
G = nx.DiGraph()

# 添加节点，并赋予 'club' 属性
for _, row in members_df.iterrows():
    G.add_node(row['member_id'], club=row['club'], name=row['name'])

# 添加带权重的边
for _, row in edge_list_df.iterrows():
    G.add_edge(row['sender_id'], row['receiver_id'], weight=row['final_weight'])

# --------------------------
# 任务 2.2: 网络指标计算
# --------------------------
print("\n--- 任务 2.2: 网络指标计算 ---")

# (1) 指标 1: 度中心性 (Degree Centrality)
# 在有向图中，分为入度 (In-Degree) 和出度 (Out-Degree)
in_degree = G.in_degree(weight='weight')  # 按权重计算（信息接收强度）
out_degree = G.out_degree(weight='weight')  # 按权重计算（信息发送强度）

# 将指标存回节点属性
nx.set_node_attributes(G, dict(in_degree), 'in_degree_weighted')
nx.set_node_attributes(G, dict(out_degree), 'out_degree_weighted')

# (2) 指标 2: 介数中心性 (Betweenness Centrality)
# 计算时考虑权重（权重代表强度，路径长度=1/强度）
# 注意：在大型网络中计算介数非常耗时
print("正在计算介数中心性 (可能需要一点时间)...")
# weight='weight' 在介数中通常指“成本”或“距离”。
# 我们的 'final_weight' 是强度，所以距离应为 1/weight。
# 为避免除零和简化，我们先计算无权重的介数，这在SNA中也很常见
betweenness = nx.betweenness_centrality(G, normalized=True, weight=None)
nx.set_node_attributes(G, betweenness, 'betweenness')
print("介数中心性计算完成。")

# (3) 指标 3: 网络密度 (Network Density)
# 对于有向图 D = E / (N * (N-1))
density = nx.density(G)
print(f"全网络密度: {density:.4f}")

# (4) 分社团计算指标（用于对比）
club_graphs = {}
club_densities = {}
club_avg_degrees = {}
clubs = members_df['club'].unique()

for club in clubs:
    # 提取子图
    club_members = members_df[members_df['club'] == club]['member_id'].tolist()
    sub_G = G.subgraph(club_members)
    club_graphs[club] = sub_G

    # 计算子图密度
    club_density = nx.density(sub_G)
    club_densities[club] = club_density

    # 计算子图平均（加权）度
    if sub_G.number_of_nodes() > 0:
        avg_degree = sum(dict(sub_G.degree(weight='weight')).values()) / sub_G.number_of_nodes()
        club_avg_degrees[club] = avg_degree
    else:
        club_avg_degrees[club] = 0

print("\n分社团指标：")
print("社团密度:", club_densities)
print("社团平均加权度:", club_avg_degrees)

# --------------------------
# 任务 2.3: 可视化图表
# --------------------------
print("\n--- 任务 2.3: 生成可视化图表 ---")

# (1) 图表一：网络社群结构图 (Network Graph)
plt.figure(figsize=(16, 12))

# 节点颜色映射
color_map = {'人工智能社': '#FF6B6B', '环境治理社': '#6BFFB8', '文化传播社': '#6B8EFF'}
node_colors = [color_map[G.nodes[node]['club']] for node in G.nodes()]

# 节点大小映射（使用出度中心性，代表影响力）
# 我们使用无权重的出度，以便更清晰地展示结构
node_sizes = [G.out_degree(node) * 50 + 10 for node in G.nodes()]  # 基础大小10，乘以系数

# 边权重映射（透明度）
edge_weights = nx.get_edge_attributes(G, 'weight').values()
max_weight = max(edge_weights) if edge_weights else 1.0
edge_alphas = [weight / max_weight * 0.5 for weight in edge_weights]  # 用透明度表示权重

# 布局
print("正在计算网络布局 (Fruchterman-Reingold)...")
# pos = nx.spring_layout(G, k=0.3, iterations=50) # spring_layout 较常用
pos = nx.fruchterman_reingold_layout(G, seed=42)  # F-R 布局
print("布局计算完成。")

# 绘制
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
nx.draw_networkx_edges(G, pos, edge_color='#AAAAAA', alpha=0.3, width=1.0)  # 使用统一的细线，避免混乱
# nx.draw_networkx_edges(G, pos, edge_color='#AAAAAA', alpha=edge_alphas) # 使用透明度（如果图不乱）

# 添加图例
legend_handles = [plt.Line2D([0], [0], marker='o', color='w', label=club,
                             markerfacecolor=color, markersize=10) for club, color in color_map.items()]
plt.legend(handles=legend_handles, title="社团 (节点颜色)", loc="upper right")

plt.title("学术兴趣社群整体网络图 (节点大小: 出度中心性)", fontsize=20,
          fontproperties=fm.FontProperties(size=20) if 'my_font' in locals() else None)
plt.axis('off')
plt.savefig("network_graph.png")
print("图表一 (network_graph.png) 已保存。")
# plt.show()


# (2) 图表二：社团网络密度对比柱状图
plt.figure(figsize=(10, 6))

club_names = list(club_densities.keys())
densities = list(club_densities.values())
colors = [color_map[name] for name in club_names]

bars = plt.bar(club_names, densities, color=colors, edgecolor='black')
plt.ylabel("网络密度 (Network Density)", fontsize=12,
           fontproperties=fm.FontProperties(size=12) if 'my_font' in locals() else None)
plt.title("三个社团内部网络密度对比", fontsize=16,
          fontproperties=fm.FontProperties(size=16) if 'my_font' in locals() else None)
plt.xticks(fontproperties=fm.FontProperties(size=12) if 'my_font' in locals() else None)

# 添加数据标签
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.001, f'{yval:.4f}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig("density_comparison.png")
print("图表二 (density_comparison.png) 已保存。")
# plt.show()

print("\n--- 任务 2 全部完成 ---")