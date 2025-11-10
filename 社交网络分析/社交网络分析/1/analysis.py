#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
社交网络分析：关键用户识别与社区发现

本脚本将执行以下任务：
1. 加载社交网络数据（如果 'social_network.csv' 不存在，则创建一个模拟数据）
2. 构建有向图并进行基本统计
3. 计算中心性指标，识别关键用户
4. 使用 Louvain 算法进行社区发现
5. 可视化网络图（节点大小=入度，节点颜色=社区）
6. 分析特定问题

所需库 (请确保已安装):
pip install pandas networkx python-louvain matplotlib
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from collections import Counter
import community  # 这是 'python-louvain' 库
import matplotlib
import subprocess

# --- 1. 数据预处理 ---

def load_and_preprocess_data(filename='social_network.csv'):
    """
    读取数据集，构建有向图，并统计基本特征。
    """
    print("\n--- 1. 数据预处理 ---")

    # 检查数据文件是否存在
    if not os.path.exists(filename):
        other_script_name = "generate_data.py"
        print(f"\n--- 正在调用外部脚本: {other_script_name} ---")
        try:
            # 'check=True' 会检查脚本是否成功运行 (返回码 0)
            # 如果 'other_script.py' 运行失败, 'analysis.py' 也会在这里停止
            subprocess.run(["python", other_script_name], check=True)

            print(f"--- {other_script_name} 运行完成，返回 'analysis.py' ---")

        except FileNotFoundError:
            print(f"错误: 找不到脚本 '{other_script_name}'。")
            print("请确保它和 'analysis.py' 在同一个文件夹中。")
            return  # 如果找不到文件，可以选择停止主脚本
        except subprocess.CalledProcessError as e:
            print(f"错误: '{other_script_name}' 运行时失败: {e}")
            return  # 如果子脚本出错，可以选择停止主脚本

    # 读取数据集
    try:
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"错误: 无法读取 '{filename}'。")
        return None

    # 构建有向图
    G = nx.from_pandas_edgelist(df, 'source', 'target', create_using=nx.DiGraph())

    # 确保图中包含所有 1-500 的用户 (即使某些用户没有边)
    all_nodes = set(range(1, 501))
    G.add_nodes_from(all_nodes)

    # 统计网络基本特征
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # 平均入度/出度
    # G.in_degree() 返回 (node, degree) 对的迭代器
    avg_in_degree = sum(d for n, d in G.in_degree()) / num_nodes
    avg_out_degree = sum(d for n, d in G.out_degree()) / num_nodes

    print(f"总用户数 (节点): {num_nodes}")
    print(f"总关注关系数 (边): {num_edges}")
    print(f"平均入度: {avg_in_degree:.4f}")
    print(f"平均出度: {avg_out_degree:.4f}")

    return G


# --- 2. 关键用户识别 ---

def identify_key_users(G):
    """
    计算入度中心性和中介中心性，并识别 Top 10 用户。
    """
    print("\n--- 2. 关键用户识别 ---")

    # a. 入度中心性 (反映被关注程度)
    # 字典: {node: centrality}
    in_degree_centrality = nx.in_degree_centrality(G)

    # b. 中介中心性 (反映作为“桥梁”的能力)
    # 注意：计算中介中心性可能需要一些时间
    print("正在计算中介中心性 (可能需要几秒钟)...")
    betweenness_centrality = nx.betweenness_centrality(G)
    print("计算完成。")

    # c. 辅助函数：排序并获取 Top N
    def get_top_n(centrality_dict, n=10):
        # 将字典按值（value）降序排序
        sorted_items = sorted(centrality_dict.items(), key=lambda item: item[1], reverse=True)
        return sorted_items[:n]

    # d. 获取并打印 Top 10
    top_10_in_degree = get_top_n(in_degree_centrality, 10)
    top_10_betweenness = get_top_n(betweenness_centrality, 10)

    print("\n** 入度中心性 Top 10 (意见领袖/受关注者) **")
    for node, score in top_10_in_degree:
        print(f"用户 {node}: {score:.4f} (原始入度: {G.in_degree(node)})")

    print("\n** 中介中心性 Top 10 (桥梁/连接者) **")
    for node, score in top_10_betweenness:
        print(f"用户 {node}: {score:.4f}")

    print("\n** 角色分析 **")
    print("入度中心性高的用户是网络中的 '名人'，被大量用户直接关注。")
    print("中介中心性高的用户是 '桥梁'，他们连接了不同的群体，在信息跨群体传播中至关重要。")

    # 返回中心性结果以供后续使用
    return in_degree_centrality, betweenness_centrality


# --- 3. 社区发现 ---

def find_communities(G):
    """
    使用 Louvain 算法（在无向图上）进行社区划分。
    """
    print("\n--- 3. 社区发现 (Louvain 算法) ---")

    # 任务要求：忽略边的方向，按无向图处理
    # 我们创建一个无向图副本 G_undirected
    G_undirected = G.to_undirected()

    # 使用 Louvain 算法 (community.best_partition)
    # 它返回一个字典 {node: community_id}
    print("正在运行 Louvain 算法进行社区划分...")

    partition = community.best_partition(G_undirected)

    print("社区划分完成。")

    # 统计每个社区的用户数量
    community_counts = Counter(partition.values())

    print(f"总共发现 {len(community_counts)} 个社区。")

    # 识别最大的 3 个社区
    top_3_communities = community_counts.most_common(3)

    print("\n** 最大的 3 个社区 **")
    for i, (comm_id, count) in enumerate(top_3_communities):
        print(f"排名 {i + 1}: 社区 {comm_id}, 用户数: {count}")

    # 返回分区结果和 Top 3 列表
    return partition, top_3_communities


# --- 4. 可视化与分析 ---

def visualize_and_analyze(G, in_degree_centrality, betweenness_centrality, partition, top_3_communities):
    """
    绘制网络结构图并回答分析问题。
    """
    print("\n--- 4. 可视化与分析 ---")

    # --- 4a. 可视化 ---
    print("准备网络可视化 (500 个节点可能较慢且密集)...")

    # 1. 布局 (使用无向图进行布局，效果通常更好)
    # spring_layout 是常用的布局算法
    G_undirected = G.to_undirected()
    pos = nx.spring_layout(G_undirected, k=0.15, iterations=50, seed=42)

    # 2. 节点大小 (按入度中心性)
    min_size = 20
    max_size_scale = 100

    node_sizes = [G.in_degree(node) * max_size_scale + min_size for node in G.nodes()]

    # 3. 节点颜色 (按社区)
    # 获取社区 ID 的数量
    num_communities = len(set(partition.values()))

    # -------------------
    #  *** 修改点 1: 获取 colormap ***
    # -------------------
    # 新 API 不再接受第二个参数
    cmap = matplotlib.colormaps.get_cmap('tab20')

    # -------------------
    #  *** 修改点 2: 应用 colormap ***
    # -------------------
    # 我们使用模运算 (%) 来循环使用 colormap 中的颜色 ('tab20' 有 20 种颜色, cmap.N == 20)
    node_colors = [cmap(partition[node] % cmap.N) for node in G.nodes()]

    # 4. 绘制
    plt.figure(figsize=(18, 18))

    nx.draw_networkx(
        G,
        pos,
        node_size=node_sizes,
        node_color=node_colors,
        with_labels=False,  # 500 个节点太多，不显示标签
        width=0.1,  # 边线宽度
        alpha=0.7,  # 节点透明度
        arrowsize=5  # 箭头大小 (因为是有向图)
    )

    plt.title("社交网络分析：社区（颜色）与影响力（节点大小=入度）", fontsize=24)
    plt.axis('off')  # 关闭坐标轴

    # 5. 保存图像
    output_filename = "social_network_visualization.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"可视化结果已保存为 '{output_filename}'")
    # plt.show()

    # --- 4b. 分析问题 ---
    print("\n** 分析问题 **")
    print("问题：最大社区中的核心用户（入度最高）是否同时是跨社区的“桥梁”（中介中心性较高）？")

    # 1. 找到最大社区 ID
    largest_comm_id = top_3_communities[0][0]

    # 2. 筛选出该社区的所有用户
    nodes_in_largest_comm = [node for node, comm in partition.items() if comm == largest_comm_id]

    # 3. 在该社区中找到入度最高的用户
    largest_comm_in_degrees = {node: in_degree_centrality[node] for node in nodes_in_largest_comm}
    core_user = max(largest_comm_in_degrees, key=largest_comm_in_degrees.get)
    core_user_in_degree = largest_comm_in_degrees[core_user]

    # 4. 查找该核心用户的中介中心性
    core_user_betweenness = betweenness_centrality[core_user]

    # 5. 分析该中介中心性是否“较高”
    top_50_betweenness_nodes = [node for node, score in
                                sorted(betweenness_centrality.items(), key=lambda item: item[1], reverse=True)[:50]]

    print(f"\n分析结果:")
    print(f"最大社区 ID: {largest_comm_id} (共 {len(nodes_in_largest_comm)} 个用户)")
    print(f"该社区的核心用户 (入度最高): 用户 {core_user}")
    print(f"  - 其入度中心性: {core_user_in_degree:.4f} (在社区内最高)")
    print(f"  - 其中介中心性: {core_user_betweenness:.4f}")

    if core_user in top_50_betweenness_nodes:
        print(f"结论：**是/偏是**。该核心用户 (用户 {core_user}) 的中介中心性排在全网前 50 名 (Top 10%)。")
        print("这表明他/她不仅是社区内的核心，同时也扮演着跨社区的“桥梁”角色。")
    else:
        print(f"结论：**否**。该核心用户 (用户 {core_user}) 的中介中心性未排进全网前 50 名 (Top 10%)。")
        print("这表明他/她很可能是该社区的“内部意见领袖”，但并非连接不同社区的关键“桥梁”。")

# --- 主函数 ---

def main():
    # 任务 1
    G = load_and_preprocess_data('social_network.csv')
    if G is None:
        return

    # 任务 2
    in_degree, betweenness = identify_key_users(G)

    # 任务 3
    partition, top_3_communities = find_communities(G)

    # 任务 4
    visualize_and_analyze(G, in_degree, betweenness, partition, top_3_communities)

    print("\n--- 所有任务已完成 ---")


if __name__ == "__main__":
    # 设置 matplotlib 解决中文显示问题
    try:
        # 优先尝试 'Microsoft YaHei' (微软雅黑), 适用于 Windows
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
        plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
        print("已设置中文字体为 'Microsoft YaHei'。")
    except Exception as e:
        try:
            # 如果失败, 尝试 'SimHei' (黑体)
            plt.rcParams['font.sans-serif'] = ['SimHei']
            plt.rcParams['axes.unicode_minus'] = False
            print("设置 'Microsoft YaHei' 失败, 尝试 'SimHei'。")
        except Exception as e2:
            # 如果都失败, 打印警告, 图像中的中文将无法显示
            print(f"设置中文字体失败 (SimHei / Microsoft YaHei 均未找到): {e2}")
            print("警告：图像中的中文标题可能显示为方块。")

    main()