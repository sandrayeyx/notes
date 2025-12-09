import pandas as pd
import numpy as np
import jieba
import re
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import Counter
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

# =================配置区域=================
DATA_FILE = 'product_reviews_dataset.csv'
# 图片保存文件夹名称
IMAGE_SAVE_DIR = 'result_images'

# 词典文件映射
DICT_FILES = {
    'pos_eval': '正面评价词语（中文）.txt',
    'pos_sent': '正面情感词语（中文）.txt',
    'neg_eval': '负面评价词语（中文）.txt',
    'neg_sent': '负面情感词语（中文）.txt',
    'degree': '程度级别词语（中文）.txt'
}

# 绘图设置 (解决中文乱码)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False
# 修复 seaborn 风格导致的字体问题
sns.set_style("whitegrid", {"font.sans-serif": ['SimHei', 'Arial Unicode MS', 'Microsoft YaHei']})


# ================= 工具函数 =================

def load_dict_file(filepath):
    words = set()
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('中文') or line[0].isdigit():
                    continue
                words.add(line)
    except Exception as e:
        print(f"警告: 无法读取词典文件 {filepath} ({e})")
    return words


# ================= 任务1：数据预处理 =================

def load_data():
    print(">>> 正在加载数据...")
    try:
        df = pd.read_csv(DATA_FILE)
    except FileNotFoundError:
        print(f"错误：找不到文件 {DATA_FILE}")
        return None, None, None, None, None

    if 'comment_time' in df.columns:
        df['comment_time'] = pd.to_datetime(df['comment_time'])
        df['month'] = df['comment_time'].dt.to_period('M')

    def clean(text):
        if not isinstance(text, str): return ""
        text = re.sub(r'[0-9\W_]+', ' ', text)
        words = jieba.lcut(text)
        stopwords = {'的', '了', '是', '我', '在', '和', '都', '就', '也', '很', '太', '这', '那', '有', '用'}
        return " ".join([w for w in words if w not in stopwords and len(w.strip()) > 0])

    print(">>> 正在清洗文本...")
    df['clean_content'] = df['comment_content'].apply(clean)

    df['stratify'] = df['product_category'] + "_" + df['sentiment_label']
    valid_counts = df['stratify'].value_counts()
    df_valid = df[df['stratify'].isin(valid_counts[valid_counts > 1].index)].copy()

    X_train, X_test, y_train, y_test = train_test_split(
        df_valid['clean_content'], df_valid['sentiment_label'],
        test_size=0.3, random_state=42, stratify=df_valid['stratify']
    )

    return df, X_train, X_test, y_train, y_test


# ================= 任务2：情感分析技术 =================

class DictionaryAnalyzer:
    def __init__(self):
        self.pos_words = set()
        self.neg_words = set()

    def load_resources(self):
        print(">>> [方案A] 加载词典...")
        self.pos_words.update(load_dict_file(DICT_FILES['pos_eval']))
        self.pos_words.update(load_dict_file(DICT_FILES['pos_sent']))
        self.neg_words.update(load_dict_file(DICT_FILES['neg_eval']))
        self.neg_words.update(load_dict_file(DICT_FILES['neg_sent']))
        print(f"    正面词库: {len(self.pos_words)} 个, 负面词库: {len(self.neg_words)} 个")

    def predict(self, text):
        words = text.split()
        score = 0
        for w in words:
            if w in self.pos_words: score += 1
            if w in self.neg_words: score -= 1
        return '正面' if score > 0 else ('负面' if score < 0 else '中性')


def run_ml_models(X_train, y_train, X_test):
    print(">>> [方案B] 训练机器学习模型...")
    tfidf = TfidfVectorizer(max_features=3000)
    X_train_tf = tfidf.fit_transform(X_train)
    X_test_tf = tfidf.transform(X_test)

    model = LogisticRegression(max_iter=1000)
    model.fit(X_train_tf, y_train)
    y_pred = model.predict(X_test_tf)

    return model, tfidf, y_pred


# ================= 任务3：可视化分析 (含保存功能) =================

def visualize_and_save(df, best_model_results, metrics_summary):
    """
    生成可视化图表并保存到本地文件夹
    """
    print(f"\n>>> [任务3] 开始生成可视化图表并保存至 '{IMAGE_SAVE_DIR}'...")

    # 确保保存目录存在
    if not os.path.exists(IMAGE_SAVE_DIR):
        os.makedirs(IMAGE_SAVE_DIR)
        print(f"    已创建文件夹: {IMAGE_SAVE_DIR}")

    # --- 图表 1: 情感分布可视化 ---
    fig1, axes1 = plt.subplots(1, 2, figsize=(16, 6))

    # (1) 各产品品类的情感标签分布
    ct = pd.crosstab(df['product_category'], df['sentiment_label'], normalize='index')
    ct.plot(kind='bar', stacked=True, ax=axes1[0], colormap='viridis', alpha=0.8)
    axes1[0].set_title('各品类情感分布占比')
    axes1[0].set_ylabel('占比')

    # (2) 时间维度情感趋势
    if 'month' in df.columns:
        trend = df.groupby(['month', 'sentiment_label']).size().unstack().fillna(0)
        trend.plot(ax=axes1[1], marker='o')
        axes1[1].set_title('时间维度情感趋势')

    plt.tight_layout()
    # 保存图片
    save_path1 = os.path.join(IMAGE_SAVE_DIR, '1_sentiment_distribution.png')
    plt.savefig(save_path1, dpi=300)
    print(f"    已保存: {save_path1}")
    plt.show()

    # --- 图表 2: 用户属性与情感关联 ---
    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

    sns.countplot(data=df, x='user_age_range', hue='sentiment_label', ax=axes2[0], palette='Set2')
    axes2[0].set_title('不同年龄段的情感分布')

    sns.countplot(data=df, x='user_gender', hue='sentiment_label', ax=axes2[1], palette='Set1')
    axes2[1].set_title('性别与情感关联')

    top_provinces = df['user_province'].value_counts().head(10).index
    prov_data = df[df['user_province'].isin(top_provinces)]
    prov_matrix = pd.crosstab(prov_data['user_province'], prov_data['sentiment_label'])
    sns.heatmap(prov_matrix, annot=True, fmt='d', cmap='YlGnBu', ax=axes2[2])
    axes2[2].set_title('Top10 地域情感热力图')

    plt.tight_layout()
    # 保存图片
    save_path2 = os.path.join(IMAGE_SAVE_DIR, '2_user_demographics.png')
    plt.savefig(save_path2, dpi=300)
    print(f"    已保存: {save_path2}")
    plt.show()

    # --- 图表 3: 关键词可视化 ---
    print("    正在计算关键词频率...")
    fig3, axes3 = plt.subplots(1, 2, figsize=(16, 6))

    def get_top_words(texts, n=15):
        all_words = []
        for t in texts:
            all_words.extend(t.split())
        return Counter(all_words).most_common(n)

    pos_texts = df[df['sentiment_label'] == '正面']['clean_content']
    pos_top = get_top_words(pos_texts)

    neg_texts = df[df['sentiment_label'] == '负面']['clean_content']
    neg_top = get_top_words(neg_texts)

    # 使用修复后的 seaborn 写法
    if pos_top:
        words, counts = zip(*pos_top)
        sns.barplot(x=list(counts), y=list(words), hue=list(words), legend=False, ax=axes3[0], palette='Greens_d')
        axes3[0].set_title('正面评论高频词 (Top 15)')

    if neg_top:
        words, counts = zip(*neg_top)
        sns.barplot(x=list(counts), y=list(words), hue=list(words), legend=False, ax=axes3[1], palette='Reds_d')
        axes3[1].set_title('负面评论高频词 (Top 15)')

    plt.tight_layout()
    # 保存图片
    save_path3 = os.path.join(IMAGE_SAVE_DIR, '3_keywords_analysis.png')
    plt.savefig(save_path3, dpi=300)
    print(f"    已保存: {save_path3}")
    plt.show()

    # --- 图表 4: 模型效果可视化 ---
    fig4, axes4 = plt.subplots(1, 2, figsize=(14, 6))

    metrics_df = pd.DataFrame(metrics_summary).T
    if not metrics_df.empty:
        sns.barplot(x=metrics_df.index, y=metrics_df['Accuracy'], hue=metrics_df.index, legend=False, ax=axes4[0],
                    palette='viridis')
        axes4[0].set_title('不同情感分析方案准确率对比')
        axes4[0].set_ylim(0, 1.1)
        for i, v in enumerate(metrics_df['Accuracy']):
            axes4[0].text(i, v + 0.01, f'{v:.2f}', ha='center')

    y_true = best_model_results['y_true']
    y_pred = best_model_results['y_pred']
    cm = confusion_matrix(y_true, y_pred, labels=['正面', '中性', '负面'])
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['正面', '中性', '负面'],
                yticklabels=['正面', '中性', '负面'], ax=axes4[1])
    axes4[1].set_title('最优模型混淆矩阵')

    plt.tight_layout()
    # 保存图片
    save_path4 = os.path.join(IMAGE_SAVE_DIR, '4_model_evaluation.png')
    plt.savefig(save_path4, dpi=300)
    print(f"    已保存: {save_path4}")
    plt.show()


# ================= 主程序 =================

if __name__ == "__main__":
    df_full, X_train, X_test, y_train, y_test = load_data()

    if df_full is not None:
        metrics_summary = {}

        # 方案A
        dict_analyzer = DictionaryAnalyzer()
        dict_analyzer.load_resources()
        y_pred_a = X_test.apply(dict_analyzer.predict)
        metrics_summary['Scheme A'] = {'Accuracy': accuracy_score(y_test, y_pred_a)}
        print(f"方案A 准确率: {metrics_summary['Scheme A']['Accuracy']:.4f}")

        # 方案B
        ml_model, tfidf, y_pred_b = run_ml_models(X_train, y_train, X_test)
        metrics_summary['Scheme B'] = {'Accuracy': accuracy_score(y_test, y_pred_b)}
        print(f"方案B 准确率: {metrics_summary['Scheme B']['Accuracy']:.4f}")

        # 可视化与保存
        best_results = {'y_true': y_test, 'y_pred': y_pred_b}
        visualize_and_save(df_full, best_results, metrics_summary)

        print(f"\n>>> 所有任务完成！图片已保存至目录: ./{IMAGE_SAVE_DIR}/")