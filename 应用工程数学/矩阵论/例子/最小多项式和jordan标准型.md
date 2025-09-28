[[2 Jordan标准型介绍#^b6e7f8]]

矩阵的最小多项式与 Jordan 标准形之间存在紧密的结构关联，核心关系可总结为：

### 1. 特征值与 Jordan 块的 “最高阶” 对应

对于矩阵 A 的某个特征值  $\lambda_0$ ，最小多项式中  $(\lambda - \lambda_0)^k$  的次数 k，**等于  $\lambda_0$  在 Jordan 标准形中对应的 “最大 Jordan 块的阶数”**。

例如：

- 若  $\lambda_0$  对应的 Jordan 块有 1 个 3 阶块、1 个 2 阶块，则最小多项式中  $(\lambda - \lambda_0)$  的次数为 3（最大 Jordan 块的阶数）。
- 若  $\lambda_0$  仅对应 1 阶 Jordan 块（即半单情况），则最小多项式中  $(\lambda - \lambda_0)$  的次数为 1。

### 2. 最小多项式的整体构造

设矩阵 A 有互异特征值  $\lambda_1, \lambda_2, \dots, \lambda_s$ ，且每个  $\lambda_i$  对应的最大 Jordan 块阶数为  $k_i$ ，则 A 的最小多项式为：

 $m_A(\lambda) = (\lambda - \lambda_1)^{k_1} (\lambda - \lambda_2)^{k_2} \cdots (\lambda - \lambda_s)^{k_s}$ 

### 3. 零化性与 Jordan 块的幂次

Jordan 块的结构决定了其 “幂零程度”：对于一个 n 阶 Jordan 块 J（特征值为  $\lambda$ ），有  $(J - \lambda I)^n = 0$ ，但  $(J - \lambda I)^{n-1} \neq 0$ 。因此，单个 Jordan 块的最小多项式是  $(\lambda - \lambda_0)^n$ （n 为块的阶数）。

由于 Jordan 标准形是**分块对角矩阵**，其最小多项式是各 Jordan 块最小多项式的**最小公倍式**。而相似矩阵（A 与它的 Jordan 标准形相似）具有相同的最小多项式，因此 A 的最小多项式由 “各特征值对应的最大 Jordan 块的阶数” 共同决定。

[[2 Jordan标准型介绍#^b6e7f8]]