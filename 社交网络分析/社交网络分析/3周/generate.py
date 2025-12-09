import pandas as pd
import numpy as np
import random
from datetime import datetime, timedelta
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# --------------------------
# 零、全局配置（字体）
# --------------------------
# 确保中文字体显示正常
try:
    plt.rcParams['font.sans-serif'] = ['Source Han Sans CN', 'SimHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
except Exception as e:
    print(f"中文字体设置警告: {e}。图表中的中文可能显示为方框。")
    print("请尝试安装 'Source Han Sans CN' 或 'SimHei' 字体，或在代码中指定本地字体路径。")
    # 例如: my_font = fm.FontProperties(fname="/path/to/your/font.ttf")
    # 然后在绘图时使用 fontproperties=my_font

# ======================================================================
#
# 第一部分：数据模拟与生成（来自您的 Prompt）
#
# ======================================================================

# --------------------------
# 1. 基础数据：社团与成员信息
# --------------------------
print("--- [第1部分] 开始生成模拟数据 ---")
# 社团信息（3个学术社团）
clubs = ["人工智能社", "环境治理社", "文化传播社"]
# 生成成员数据（每个社团100人，含唯一ID、姓名、社团）
members = []
for club in clubs:
    for i in range(100):  # 每个社团100名成员
        member_id = f"{club[0]}{i + 1:03d}"  # 生成如"人001"的唯一ID
        name = f"{club[2:]}{i + 1}"  # 生成如"智能1"的简化姓名
        members.append([member_id, name, club])

# 转为DataFrame（成员信息表）
members_df = pd.DataFrame(members, columns=["member_id", "name", "club"])


# --------------------------
# 2. 模拟线上互动数据（微信/QQ群）
# --------------------------
def generate_online_data(members_df, start_date, end_date, n_records=5000):
    """生成线上互动数据（含群聊、私聊）"""
    # 时间范围：随机生成[start_date, end_date]内的时间
    date_range = pd.date_range(start=start_date, end=end_date, freq="T")  # 按分钟级
    online_data = []

    for _ in range(n_records):
        # 随机选择互动双方（可同社团或跨社团，但同社团概率更高）
        sender = random.choice(members_df["member_id"].tolist())
        receiver = random.choice(members_df["member_id"].tolist())
        while sender == receiver:  # 排除自互动
            receiver = random.choice(members_df["member_id"].tolist())

        # 互动类型：群聊发言、@提及、文件分享、私聊
        interact_type = random.choices(
            ["群聊发言", "@提及", "文件分享", "私聊"],
            weights=[0.4, 0.3, 0.1, 0.2],  # 群聊和@更频繁
            k=1
        )[0]

        # 互动时间
        interact_time = random.choice(date_range)

        # 互动权重（1-5，代表强度）
        weight = random.randint(1, 5)

        # 加入噪声：10%概率生成重复数据，5%概率缺失receiver
        if random.random() < 0.1:  # 10%重复数据
            online_data.append([sender, receiver, interact_type, interact_time, weight])
        if random.random() < 0.05:  # 5%缺失值
            receiver = np.nan

        online_data.append([sender, receiver, interact_type, interact_time, weight])

    # 转为DataFrame
    online_df = pd.DataFrame(
        online_data,
        columns=["sender_id", "receiver_id", "interact_type", "time", "weight"]
    )
    return online_df


# --------------------------
# 3. 模拟线下活动数据
# --------------------------
def generate_offline_data(members_df, n_events=15):
    """生成线下活动互动数据（学术沙龙、小组研讨等）"""
    offline_data = []
    events = [f"{club}学术沙龙第{i + 1}期" for club in clubs for i in range(5)]  # 每个社团5场活动

    for event in events:
        # 活动日期：近3个月内的周末
        event_date = (datetime.now() - timedelta(days=random.randint(1, 90))).replace(
            hour=random.choice([14, 15, 16]),  # 下午2-4点
            minute=0, second=0
        )

        # 参与人数：每场活动30-60人（随机从对应社团抽取）
        club = event.split("学术")[0]  # 从活动名提取所属社团
        participants = members_df[members_df["club"] == club]["member_id"].sample(
            random.randint(30, 60)
        ).tolist()

        # 生成互动记录：小组讨论、跨组请教、合作演示
        n_interacts = len(participants) * random.randint(3, 5)  # 每人3-5次互动
        for _ in range(n_interacts):
            a = random.choice(participants)
            b = random.choice(participants)
            while a == b:
                b = random.choice(participants)

            interact_type = random.choices(
                ["小组讨论", "跨组请教", "合作演示"],
                weights=[0.6, 0.3, 0.1],
                k=1
            )[0]

            # 互动时长（分钟）：加入异常值（1%概率出现极长时长）
            duration = random.randint(5, 30)
            if random.random() < 0.01:
                duration = random.randint(120, 180)  # 异常值（2-3小时）

            offline_data.append([event, event_date, a, b, interact_type, duration])

    # 转为DataFrame
    offline_df = pd.DataFrame(
        offline_data,
        columns=["event_name", "event_date", "person_a", "person_b", "interact_type", "duration"]
    )
    return offline_df


# --------------------------
# 4. 生成并保存数据
# --------------------------
# (注意：这里不再使用 if __name__ == "__main__": 以确保在脚本中顺序执行)

# 生成线上数据（近1个月）
online_df_gen = generate_online_data(
    members_df,
    start_date=datetime.now() - timedelta(days=30),
    end_date=datetime.now(),
    n_records=5000
)

# 生成线下数据
offline_df_gen = generate_offline_data(members_df)

# 保存为CSV（供学生后续处理）
members_df.to_csv("社团成员信息.csv", index=False, encoding="utf-8-sig")
online_df_gen.to_csv("线上互动数据.csv", index=False, encoding="utf-8-sig")
offline_df_gen.to_csv("线下活动数据.csv", index=False, encoding="utf-8-sig")

print("\n数据生成完成：")
print(f"成员信息：{len(members_df)}条")
print(f"线上互动：{len(online_df_gen)}条")
print(f"线下活动：{len(offline_df_gen)}条")
print("CSV 文件已保存到本地。")

# ======================================================================
#
# 第二部分：数据预处理、表达与分析（来自我的回答）
#
# ======================================================================

# --------------------------
# 任务 2.1: 数据清洗
# --------------------------
print("\n--- [第2部分] 任务 2.1: 数据清洗开始 ---")

# (1) 加载数据 (现在这些文件必定存在)
members_df = pd.read_csv("社团成员信息.csv")
online_df = pd.read_csv("线上互动数据.csv")
offline_df = pd.read_csv("线下活动数据.csv")

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
in_degree = G.in_degree(weight='weight')  # 按权重计算（信息接收强度）
out_degree = G.out_degree(weight='weight')  # 按权重计算（信息发送强度）
nx.set_node_attributes(G, dict(in_degree), 'in_degree_weighted')
nx.set_node_attributes(G, dict(out_degree), 'out_degree_weighted')

# (2) 指标 2: 介数中心性 (Betweenness Centrality)
print("正在计算介数中心性 (可能需要一点时间)...")
betweenness = nx.betweenness_centrality(G, normalized=True, weight=None)
nx.set_node_attributes(G, betweenness, 'betweenness')
print("介数中心性计算完成。")

# (3) 指标 3: 网络密度 (Network Density)
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
node_sizes = [G.out_degree(node) * 50 + 10 for node in G.nodes()]  # 基础大小10，乘以系数

# 布局
print("正在计算网络布局 (Fruchterman-Reingold)...")
pos = nx.fruchterman_reingold_layout(G, seed=42)
print("布局计算完成。")

# 绘制
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
nx.draw_networkx_edges(G, pos, edge_color='#AAAAAA', alpha=0.3, width=1.0)

# 添加图例
legend_handles = [plt.Line2D([0], [0], marker='o', color='w', label=club,
                             markerfacecolor=color, markersize=10) for club, color in color_map.items()]
plt.legend(handles=legend_handles, title="社团 (节点颜色)", loc="upper right")

plt.title("学术兴趣社群整体网络图 (节点大小: 出度中心性)", fontsize=20)
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
plt.ylabel("网络密度 (Network Density)", fontsize=12)
plt.title("三个社团内部网络密度对比", fontsize=16)

# 添加数据标签
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.001, f'{yval:.4f}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig("density_comparison.png")
print("图表二 (density_comparison.png) 已保存。")
# plt.show()

print("\n--- 脚本执行完毕 ---")