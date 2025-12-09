import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math

# ==========================================
# 1. 数据录入 (虚拟数据集)
# ==========================================

# 3.1 学者基础信息 & 3.3 引用指标
# 将两个表合并录入
data_scholars = {
    'ID': ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10'],
    'Area': ['AI', 'AI', 'DM', 'DM', 'CN', 'CN', 'AI', 'DM', 'CN', 'AI'],
    'Total_Papers': [22, 18, 25, 15, 20, 16, 19, 21, 17, 23],
    'Top_Papers': [10, 8, 9, 6, 7, 5, 8, 9, 6, 10],
    'Total_Citations': [920, 750, 1050, 480, 880, 620, 720, 950, 550, 900],
    'High_Cited_Papers': [6, 4, 8, 2, 7, 3, 4, 7, 3, 6],
    'Cross_Cit_Ratio': [0.35, 0.18, 0.52, 0.12, 0.48, 0.20, 0.22, 0.45, 0.15, 0.30], # 百分比转小数
    'Max_Single_Cit': [128, 95, 156, 68, 142, 85, 92, 135, 72, 118]
}

# 3.2 学者合作关系
# 格式：(学者A, 学者B, 合作总数, 跨方向合作数)
data_edges = [
    ('S1', 'S2', 5, 0),
    ('S1', 'S3', 3, 3),
    ('S2', 'S7', 4, 0),
    ('S3', 'S8', 6, 0),
    ('S3', 'S5', 2, 2),
    ('S5', 'S6', 5, 0),
    ('S5', 'S8', 3, 3),
    ('S7', 'S10', 4, 0),
    ('S8', 'S10', 2, 2),
    ('S9', 'S5', 3, 0),
    ('S4', 'S3', 2, 0)
]

df = pd.DataFrame(data_scholars)

# ==========================================
# 2. 数据预处理与指标计算
# ==========================================

# 2.1 h指数估算 (任务1.2)
# 逻辑：h指数 <= sqrt(总被引) 且 <= 发表总数，同时 >= 已知高被引论文数
def estimate_h_index(row):
    h_theoretical = int(math.sqrt(row['Total_Citations']))
    h_cap = min(row['Total_Papers'], h_theoretical)
    return max(row['High_Cited_Papers'], h_cap)

df['h_index_est'] = df.apply(estimate_h_index, axis=1)

# 2.2 跨学科影响力基础数据 (任务1.3)
# 计算每个学者的"加权成果数"：同方向 + 2*跨方向
# 先建立一个字典存储每个人的加权分
weighted_output_map = {sid: 0 for sid in df['ID']}

for u, v, total, cross in data_edges:
    same = total - cross
    score = same + (cross * 2)
    weighted_output_map[u] += score
    weighted_output_map[v] += score

df['Weighted_Output'] = df['ID'].map(weighted_output_map)

# ==========================================
# 3. 网络构建与中心性计算 (任务2.1)
# ==========================================

G = nx.Graph()

# 添加节点和属性
for idx, row in df.iterrows():
    G.add_node(row['ID'], area=row['Area'])

# 添加边和权重
# 权重设定：使用合作总数作为连接强度 (weight)，用于计算加权度
for u, v, total, cross in data_edges:
    G.add_edge(u, v, weight=total, is_cross=(cross > 0))

# 计算中心性指标
# 1. 度中心性 (Degree): 连接的广泛程度
degree_dict = nx.degree_centrality(G)
# 2. 中介中心性 (Betweenness): 桥梁作用 (使用无权图计算拓扑位置)
betweenness_dict = nx.betweenness_centrality(G)
# 3. 特征向量中心性 (Eigenvector): 核心圈层地位
eigenvector_dict = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)

df['Degree_Cent'] = df['ID'].map(degree_dict)
df['Betweenness_Cent'] = df['ID'].map(betweenness_dict)
df['Eigenvector_Cent'] = df['ID'].map(eigenvector_dict)

# ==========================================
# 4. 综合影响力指标计算 (任务2.2 & 2.3)
# ==========================================

# 4.1 综合引用影响力得分 (归一化后加权)
# 权重假设：总被引(40%) + 高被引(40%) + h指数(20%)
def normalize(series):
    return (series - series.min()) / (series.max() - series.min())

df['Norm_Cit'] = normalize(df['Total_Citations'])
df['Norm_High'] = normalize(df['High_Cited_Papers'])
df['Norm_h'] = normalize(df['h_index_est'])

df['Citation_Score'] = 0.4*df['Norm_Cit'] + 0.4*df['Norm_High'] + 0.2*df['Norm_h']

# 4.2 跨学科影响力得分
# 公式：加权成果数 * 跨学科被引占比
df['Cross_Impact_Score'] = df['Weighted_Output'] * df['Cross_Cit_Ratio']

# ==========================================
# 5. 结果排名与展示
# ==========================================

# 生成排名
df['Rank_Citation'] = df['Citation_Score'].rank(ascending=False)
df['Rank_Cross'] = df['Cross_Impact_Score'].rank(ascending=False)

print("--- 最终分析结果表 (部分列) ---")
cols_to_show = ['ID', 'Area', 'h_index_est', 'Degree_Cent', 'Betweenness_Cent',
                'Citation_Score', 'Rank_Citation', 'Cross_Impact_Score', 'Rank_Cross']
print(df[cols_to_show].sort_values('Rank_Citation').to_string(index=False))

# ==========================================
# 6. 网络可视化绘图 (重点)
# ==========================================

plt.figure(figsize=(12, 10))

# 布局算法
pos = nx.spring_layout(G, k=0.8, seed=42)  # k值调节节点间距

# 6.1 节点绘制
# 颜色映射：不同领域不同颜色
color_map = {'AI': '#FF6B6B', 'DM': '#4ECDC4', 'CN': '#45B7D1'}
node_colors = [color_map[G.nodes[n]['area']] for n in G.nodes()]

# 大小映射：根据"综合引用影响力"调整节点大小
# 将分数放大以便观察
node_sizes = [df[df['ID']==n]['Citation_Score'].values[0] * 3000 + 500 for n in G.nodes()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, edgecolors='gray')

# 6.2 标签绘制
nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', font_color='white')

# 6.3 边绘制
# 区分跨学科合作 (虚线) 和 同学科合作 (实线)
cross_edges = [(u, v) for u, v, d in G.edges(data=True) if d['is_cross']]
same_edges = [(u, v) for u, v, d in G.edges(data=True) if not d['is_cross']]

# 边的粗细对应权重 (合作数量)
weights = nx.get_edge_attributes(G, 'weight')

nx.draw_networkx_edges(G, pos, edgelist=same_edges, width=[weights[e] for e in same_edges],
                       edge_color='gray', style='solid', alpha=0.6)
nx.draw_networkx_edges(G, pos, edgelist=cross_edges, width=[weights[e] for e in cross_edges],
                       edge_color='#FFD93D', style='dashed', alpha=0.8) # 跨学科用黄色虚线高亮

# 图例制作
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='AI Scholar', markerfacecolor='#FF6B6B', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='DM Scholar', markerfacecolor='#4ECDC4', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='CN Scholar', markerfacecolor='#45B7D1', markersize=10),
    Line2D([0], [0], color='gray', lw=2, label='Intra-disciplinary (Solid)'),
    Line2D([0], [0], color='#FFD93D', lw=2, linestyle='--', label='Cross-disciplinary (Dashed)')
]
plt.legend(handles=legend_elements, loc='upper right')

plt.title("Scholar Collaboration & Influence Network Analysis", fontsize=15)
plt.axis('off')

# 保存图片
plt.savefig("scholar_network_analysis.png", dpi=300, bbox_inches='tight')
plt.show()

print("\n图表已生成：scholar_network_analysis.png")
print("节点大小代表引用影响力，虚线代表跨学科合作。")