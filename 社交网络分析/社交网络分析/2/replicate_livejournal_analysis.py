import json
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from haversine import haversine
from collections import defaultdict, Counter
import random
import math

# --- 0. 中文字体设置 (来自我们之前的修正) ---
# 解决 matplotlib 中文显示问题
try:
    plt.rcParams['font.sans-serif'] = ['SimHei']  # Windows (黑体)
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
except Exception:
    try:
        plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']  # Mac
        plt.rcParams['axes.unicode_minus'] = False
    except Exception:
        print("警告：未找到中文字体 'SimHei' 或 'Arial Unicode MS'。")
        print("图表中的中文可能显示为方块。")


# --- 1. 导入和加载数据 (来自 image_348263) ---
def load_data(user_file='livejournal_data.json', city_file='us_cities.json'):
    """
    加载 LiveJournal 用户/朋友数据和美国城市数据
    """
    print("正在加载数据...")
    try:
        with open(user_file, 'r') as f:
            data = json.load(f)
            users = data['users']  # 假设 {user_id: {location: "city_id", ...}}
            friendships = data['friendships']  # 假设 {user_id: [friend_id, ...]}

        with open(city_file, 'r') as f:
            cities = json.load(f)  # 假设 {city_id: {latitude: ..., longitude: ..., population: ...}}

        print(f"加载了 {len(users)} 个用户, {len(cities)} 个城市。")
        return users, friendships, cities
    except FileNotFoundError as e:
        print(f"错误: 必需的数据文件未找到: {e.filename}")
        print("请确保 'livejournal_data.json' 和 'us_cities.json' 在同一目录下。")
        return None, None, None
    except Exception as e:
        print(f"加载数据时出错: {e}")
        return None, None, None


# --- 2. 构建网络图 (来自 image_348263) ---
def build_graph(users, friendships, cities):
    """
    使用 NetworkX 构建图
    """
    print("正在构建网络图...")
    G = nx.DiGraph()  # 论文提到链接是有向的

    for user_id, user_info in users.items():
        city_id = user_info.get('location')
        if city_id and city_id in cities:
            # *** BUG 修复 ***
            # haversine 库期望 (lat, lon)，而不是 (lon, lat)
            # 修正前: pos=(cities[city_id]['longitude'], cities[city_id]['latitude'])
            G.add_node(user_id,
                       city_id=city_id,
                       pos=(cities[city_id]['latitude'], cities[city_id]['longitude']))

    # 仅添加图中已存在的节点之间的边
    for u, friends in friendships.items():
        if u in G:
            # 确保 friend 是字符串列表 (来自模拟数据)
            if isinstance(friends, list):
                for v in friends:
                    if v in G:
                        G.add_edge(u, v)

    # 移除没有位置信息或没有朋友的孤立节点
    isolated = [n for n in G.nodes() if G.degree(n) == 0]
    G.remove_nodes_from(isolated)

    print(f"网络图构建完成: {G.number_of_nodes()} 个节点, {G.number_of_edges()} 条边。")
    return G


# --- 3. 路由算法 (来自 image_348263) ---
def geogreedy_routing(G, source, target_user, max_steps=100):
    """
    实现 GEOGREEDY 算法 (论文核心)
    'target_user' 是目标用户的 ID
    """
    # 目标位置是目标用户的 (lat, lon)
    target_pos = G.nodes[target_user]['pos']
    # 目标城市 ID 是目标用户的城市 ID
    target_city_id = G.nodes[target_user]['city_id']

    path = [source]
    current = source

    for _ in range(max_steps):
        current_pos = G.nodes[current]['pos']

        # haversine 库函数现在接收正确的 (lat, lon) 元组
        current_dist = haversine(current_pos, target_pos)

        # 成功: 到达目标城市
        if G.nodes[current]['city_id'] == target_city_id:
            return path, "success"

        friends = list(G.successors(current))
        if not friends:
            return path, "dead_end"  # 没有朋友

        best_friend = None
        min_dist = current_dist

        # 寻找最近的朋友
        for friend in friends:
            friend_pos = G.nodes[friend]['pos']
            dist = haversine(friend_pos, target_pos)
            if dist < min_dist:
                min_dist = dist
                best_friend = friend

        # 失败: 陷入局部最小值
        if best_friend is None:
            return path, "local_minimum"

        current = best_friend
        path.append(current)

    return path, "max_steps"  # 超过最大步数


def enhanced_geogreedy_routing(G, source, target_user, max_steps=100):
    """
    论文中的“修改版”实验：如果GEOGREEDY失败，随机跳到同城
    """
    path, reason = geogreedy_routing(G, source, target_user, max_steps)

    if reason == "local_minimum":
        # 路由失败，尝试跳转到同城随机用户
        current = path[-1]
        current_city_id = G.nodes[current]['city_id']

        # 找到所有在同一个城市的人
        city_residents = [n for n in G.nodes() if G.nodes[n]['city_id'] == current_city_id and n != current]

        if city_residents:
            # 随机选择一个同城的人
            next_hop = random.choice(city_residents)
            # 从这个人开始，继续尝试路由
            new_path, new_reason = geogreedy_routing(G, next_hop, target_user, max_steps - len(path))
            return path + new_path, new_reason

    return path, reason


def run_routing_experiment(G, routing_func, n_trials=5000):
    """
    运行 N 次路由实验并报告统计数据
    """
    print(f"--- 正在运行路由实验 ({routing_func.__name__}) ---")
    successful_trials = 0
    successful_lengths = []

    # 随机选择源和目标（论文是随机选源和目标用户）
    all_users = list(G.nodes())
    if len(all_users) < 2:
        print("错误：图中用户不足2个，无法运行实验。")
        return 0, 0, []

    for _ in range(n_trials):
        source, target = random.sample(all_users, 2)
        # 传递目标 *用户ID*
        path, reason = routing_func(G, source, target, max_steps=50)

        if reason == "success":
            successful_trials += 1
            successful_lengths.append(len(path) - 1)  # 步数

    success_rate = successful_trials / n_trials if n_trials > 0 else 0
    avg_length = np.mean(successful_lengths) if successful_lengths else 0

    return success_rate, avg_length, successful_lengths


# --- 4. 基于排名的友谊 (来自 image_348263) ---
def analyze_rank_based_friendship(G, sample_size=100):
    """
    验证 P(r) ∝ 1/r (论文核心贡献 2)
    """
    print(f"--- 正在分析基于排名的友谊 (采样 {sample_size} 个节点) ---")

    all_users = list(G.nodes())
    if not all_users:
        print("错误：图中没有用户。")
        return [], []

    sample_users = random.sample(all_users, min(sample_size, len(all_users)))

    rank_probs = defaultdict(lambda: [0, 0])  # rank -> [friend_count, total_count]

    for i, u in enumerate(sample_users):
        if (i + 1) % 10 == 0:
            print(f"  正在处理采样节点 {i + 1}/{sample_size}...")

        u_pos = G.nodes[u]['pos']
        friends_of_u = set(G.successors(u))

        # 计算 u 到所有其他 v 的距离
        distances = []
        for v in all_users:
            if u == v:
                continue
            v_pos = G.nodes[v]['pos']
            distances.append((v, haversine(u_pos, v_pos)))  # haversine 正常工作

        # 按距离排序
        distances.sort(key=lambda x: x[1])

        # 计算排名 (rank = 排名比我近的人数)
        ranks = {}
        current_rank = 0
        last_dist = -1
        for v, dist in distances:
            if dist > last_dist:
                current_rank = len(ranks)
            ranks[v] = current_rank
            last_dist = dist

        # 统计友谊
        for v, rank in ranks.items():
            if rank < 1:  # 忽略 rank 0
                continue
            rank_probs[rank][1] += 1  # total_count
            if v in friends_of_u:
                rank_probs[rank][0] += 1  # friend_count

    # 计算最终概率
    ranks = []
    probabilities = []
    for rank, (friend_count, total_count) in rank_probs.items():
        if total_count > 0:
            prob = friend_count / total_count
            if prob > 0:  # 避免 log(0)
                ranks.append(rank)
                probabilities.append(prob)

    print("分析完成，正在生成图表...")
    return ranks, probabilities


# --- 5. 网络结构分析 (来自 image_348263 和 348581) ---
def analyze_network_structure(G):
    """
    分析网络的基本属性：度分布，聚类系数等
    """
    print("--- 正在分析网络结构 ---")

    if G.number_of_nodes() == 0:
        print("  网络为空，跳过结构分析。")
        return

    # 度分布
    in_degrees = [d for n, d in G.in_degree()]
    out_degrees = [d for n, d in G.out_degree()]

    print(f"  平均入度: {np.mean(in_degrees):.2f}")
    print(f"  平均出度: {np.mean(out_degrees):.2f}")

    # 聚类系数 (论文为 0.2)
    undirected_G = G.to_undirected()
    avg_clustering = nx.average_clustering(undirected_G)
    print(f"  平均聚类系数 (无向图估算): {avg_clustering:.3f}")

    # 绘制度分布 (log-log)
    plt.figure(figsize=(12, 6))

    def plot_degree_dist(degrees, title):
        if not degrees:
            return
        counts = Counter(degrees)
        deg, cnt = zip(*counts.items())
        plt.scatter(deg, cnt, alpha=0.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.title(title)
        plt.xlabel("度 (k)")
        plt.ylabel("N(k)")
        plt.grid(True, which="both", ls="--", alpha=0.4)

    plt.subplot(1, 2, 1)
    plot_degree_dist(in_degrees, "入度分布 (Log-Log)")

    plt.subplot(1, 2, 2)
    plot_degree_dist(out_degrees, "出度分布 (Log-Log)")

    plt.suptitle("网络度分布 (复现图1)")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show(block=False)

    # 路径长度 (论文图2)
    print("  正在采样计算最短路径长度...")
    try:
        largest_cc_nodes = max(nx.connected_components(undirected_G), key=len)
        subgraph_nodes = random.sample(list(largest_cc_nodes), min(100, len(largest_cc_nodes)))
        subgraph = G.subgraph(subgraph_nodes)

        path_lengths = []
        for source in subgraph.nodes():
            for target in subgraph.nodes():
                if source != target:
                    try:
                        l = nx.shortest_path_length(subgraph, source=source, target=target)
                        path_lengths.append(l)
                    except nx.NetworkXNoPath:
                        pass

        if path_lengths:
            print(f"  平均最短路径长度 (采样): {np.mean(path_lengths):.2f}")
            plt.figure(figsize=(8, 5))
            plt.hist(path_lengths, bins=range(1, max(path_lengths) + 2), density=True, alpha=0.7)
            plt.title("最短路径长度分布 (复现图2 - 采样)")
            plt.xlabel("路径长度 k")
            plt.ylabel("f(k)")
            plt.show(block=False)

    except Exception as e:
        print(f"  计算最短路径时出错 (可能是图太小): {e}")


# --- 7. 进阶分析：多维度路由 (来自 image_348581) ---
def add_interests_to_nodes(G, num_interests=100, avg_interests=5):
    """
    为节点添加模拟的“兴趣”属性，用于多维路由
    """
    print("--- 正在为节点添加模拟兴趣 ---")
    for node in G.nodes():
        num = np.random.poisson(avg_interests)
        interests = set(random.sample(range(num_interests), min(num, num_interests)))
        G.nodes[node]['interests'] = interests


def multidimensional_routing(G, source, target_user, lambda_s=0.5, max_steps=100):
    """
    多维度路由算法
    """
    target_pos = G.nodes[target_user]['pos']
    target_interests = G.nodes[target_user].get('interests', set())
    target_city_id = G.nodes[target_user]['city_id']

    path = [source]
    current = source

    for _ in range(max_steps):
        current_pos = G.nodes[current]['pos']

        if G.nodes[current]['city_id'] == target_city_id:
            return path, "success"

        friends = list(G.successors(current))
        if not friends:
            return path, "dead_end"

        best_friend = None
        min_score = float('inf')

        current_geodist = haversine(current_pos, target_pos)

        for friend in friends:
            friend_pos = G.nodes[friend]['pos']
            geodist = haversine(friend_pos, target_pos)

            friend_interests = G.nodes[friend].get('interests', set())
            union_size = len(target_interests.union(friend_interests))
            intersect_size = len(target_interests.intersection(friend_interests))

            social_dist = 1.0 - (intersect_size / union_size) if union_size > 0 else 1.0

            # 组合分数 (非常简化的模型)
            # 我们需要一个方法来平衡地理距离和社会距离
            # 假设社会距离 1.0 约等于 1000km 的地理距离
            try:
                score = (lambda_s * geodist) + ((1 - lambda_s) * social_dist * 1000)
            except Exception:
                score = float('inf')

            if score < min_score:
                min_score = score
                best_friend = friend

        if best_friend is None or min_score >= (lambda_s * current_geodist):
            return path, "local_minimum"

        current = best_friend
        path.append(current)

    return path, "max_steps"


# --- 6. 模拟验证与论文对比 (来自 image_348581) ---
def main():
    """
    主执行函数
    """
    # 1. 加载和构建
    users, friendships, cities = load_data()
    if users is None:
        return

    G = build_graph(users, friendships, cities)
    if G.number_of_nodes() == 0:
        print("图为空，无法继续。请检查数据文件。")
        return

    # 2. 分析网络结构
    # analyze_network_structure(G) # 耗时较长，可以取消注释以运行

    # 3. 运行 GEOGREEDY 路由实验
    paper_success_rate = 0.13
    paper_avg_length = 4.12
    n_trials = 2000  # 论文为 500,000，我们用 2000

    success_rate, avg_length, lengths = run_routing_experiment(G, geogreedy_routing, n_trials)

    print("\n--- GEOGREEDY 路由结果对比 ---")
    print(f"  模拟成功率: {success_rate * 100:.2f}% (论文: {paper_success_rate * 100:.2f}%)")
    print(f"  模拟平均长度: {avg_length:.2f} (论文: {paper_avg_length:.2f}%)")

    # 绘制路径长度分布 (图2)
    if lengths:
        plt.figure(figsize=(8, 5))
        plt.hist(lengths, bins=range(1, max(lengths, default=20) + 2), density=True, alpha=0.7)
        plt.title("GEOGREEDY 成功路径长度分布 (复现图2)")
        plt.xlabel("路径长度 k")
        plt.ylabel("f(k)")
        plt.show(block=False)
    else:
        print("  没有成功的 GEOGREEDY 路径可供绘制。")

    # 4. 验证基于排名的友谊
    ranks, probs = analyze_rank_based_friendship(G, sample_size=100)  # 论文为 10,000，非常耗时

    if ranks and probs:
        log_r = np.log10(ranks)
        log_p = np.log10(probs)

        plt.figure(figsize=(10, 6))
        plt.scatter(log_r, log_p, alpha=0.3, label="数据点 (binned)")

        try:
            # 拟合
            m, c = np.polyfit(log_r, log_p, 1)
            fit_line = m * log_r + c
            plt.plot(log_r, fit_line, color='red', linestyle='--',
                     label=f"拟合线 (斜率 m = {m:.2f})")
            print(f"--- 基于排名的友谊 结果 ---")
            print(f"  log-log 图拟合斜率: {m:.2f} (论文发现 ~ -1.25)")
        except Exception as e:
            print(f"  拟合P(r) vs r时出错: {e}")

        plt.title("友谊概率 vs 排名 (复现图5)")
        plt.xlabel("Log10(Rank r) (排名)")
        plt.ylabel("Log10(Friendship Probability P(r)) (友谊概率)")
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.show(block=False)  # 改为非阻塞
    else:
        print("  没有排名数据可供绘制。")

    # 5. 运行多维度路由
    add_interests_to_nodes(G)

    multi_dim_func = lambda g, s, t, max_steps: multidimensional_routing(g, s, t, lambda_s=0.5, max_steps=max_steps)
    multi_dim_func.__name__ = "multidimensional_routing (lambda=0.5)"

    success_rate_multi, avg_length_multi, _ = run_routing_experiment(G, multi_dim_func, n_trials)

    print("\n--- 多维度路由 结果对比 ---")
    print(f"  GEOGREEDY 成功率: {success_rate * 100:.2f}%")
    print(f"  多维度路由 成功率: {success_rate_multi * 100:.2f}%")
    print(f"  多维度路由 平均长度: {avg_length_multi:.2f}")

    # --- 显示所有图表 ---
    print("\n--- 分析完成 ---")
    print("正在显示所有图表... (如果窗口未响应，请关闭图表窗口以结束程序)")
    plt.show()  # 阻塞，直到所有图表窗口被关闭


if __name__ == "__main__":
    main()