## 2.1 线性变换的对角矩阵表示

### 1 线性变换的特征值和特征向量

**定义2.1**：设 $T$ 为线性空间 $V$ 上的线性变换，若存在 $\xi \in V$ 和数 $\lambda \in F$, $\xi \neq 0$，使 $T(\xi) = \lambda \xi$，则称数 $\lambda$ 为 $T$ 的特征值，向量 $\xi$ 为线性变换 $T$ 的对应于特征值 $\lambda$ 的特征向量。

**定理 2.1** 设 $V$ 上线性变换 $T$ 在基 $\{\alpha_1, \alpha_2, \cdots, \alpha_n\}$ 下矩阵为 $A$，则 $A$ 的特征值 $\lambda$ 就是变换 $T$ 的特征值；若 $X$ 是 $A$ 的特征向量，则  
$$\xi = (\alpha_1, \alpha_2, \cdots, \alpha_n)X$$  
就是 $T$ 的特征向量。

**注 2.1** 矩阵 $A$ 是和基相关的，若 $T$ 在另一组基下矩阵为 $B$，从 §1.3 知 $B$ 相似于 $A$，即 $B = P^{-1}AP$ 或 $PB = AP$。

$$AX = \lambda X \Longleftrightarrow B(P^{-1}X) = \lambda (P^{-1}X),$$

因此 $B$ 与 $A$ 的特征值是一样的，特征向量不一样。也就是说，**$T$ 的特征值是由 $T$ 决定的，和基的选择无关。**

---

下面从线性空间的角度讨论线性变换T的特征向量的性质。

**定义 2.2** 设 $\lambda$ 为线性变换 $T$ 的特征值，$\xi_1, \xi_2, \cdots, \xi_t$ 是 $T$ 对应于 $\lambda$ 的特征向量的极大线性无关组，则称子空间 $V_\lambda = L\{\xi_1, \xi_2, \cdots, \xi_t\}$ 为 $T$ 关于 $\lambda$ 的特征子空间。
	应该注意到，$V_\lambda - \{0\}$才是线性变换T的全部特征向量；并且$V_\lambda = N(T - \lambda I)$．

如果线性变换T有s个互异的特征值:  $\lambda_1, \lambda_2, \cdots, \lambda_s$ ，它们就会对应$s$个特征子空间: $V_{\lambda_1}, V_{\lambda_2}, \cdots, V_{\lambda_s}$．特征子空间有如下性质:

**定理 2.2** 设 $\lambda_1, \lambda_2, \cdots, \lambda_s$ 是线性空间 $V_n(F)$ 上线性变换$T$的互异的特征值．$V_{\lambda_i}$ 是 $\lambda_i$ 的特征子空间，$i = 1, 2, \cdots, s$，则有
1.  $V_{\lambda_i}$ 是T的不变子空间；
2.  $\lambda_i \neq \lambda_j$ 时， $V_{\lambda_i} \cap V_{\lambda_j} = \{\mathbf{0}\}$ ；
3. 若 $\lambda_i$ 是T的 $k_i$ 重特征值，则 $\dim V_{\lambda_i} \leq k_i$ ．(几何重数<代数重数)

**注 2.2**  
1. 此定理说明
$$ V_{\lambda_1} + V_{\lambda_2} + \cdots + V_{\lambda_s} = V_{\lambda_1} \oplus V_{\lambda_2} \oplus \cdots \oplus V_{\lambda_s} \subseteq V, $$
并且$\sum_{i=1}^s \dim V_{\lambda_i} \leq \dim V = n.$

如果 $\sum_{i=1}^s \dim V_{\lambda_i} = n,$
则有
$$ V_{\lambda_1} \oplus V_{\lambda_2} \oplus \cdots \oplus V_{\lambda_s} = V. $$
2. $\dim V_{\lambda_i}$ 称为特征值 $\lambda_i$ 的几何重数，而 $k_i$ 称为 $\lambda_i$ 的代数重数。
3. 设矩阵 $A \in \mathbb{C}^{m \times n}$ ， $B \in \mathbb{C}^{n \times m}$ ，AB与BA有相同的非零特征值，且代数重数也相同

### 2 线性变换矩阵的对角化

本段主要从变换T和空间分解的角度给出线性变换T有对角矩阵表示的充分必要条件。

**定理 2.3** 线性变换 $T$ 有对角矩阵表示的充分必要条件是 $T$ 有 $n$ 个线性无关的特征向量．


前面已讨论过，特征子空间 $V_{\lambda_i}$ 是 $T$ 的不变子空间，特别地，线性变换 $T$ 在不变子空间 $V_{\lambda_i}$ 上矩阵为对角矩阵 $\lambda_i I_{t_i}$ ．从空间分解的角度，我们有下列定理:

**定理 2.4** 线性变换 $T$ 有对角矩阵表示的充分必要条件是：  
$$V_{\lambda_1} \oplus V_{\lambda_2} \oplus \cdots \oplus V_{\lambda_s} = V. \tag{2 - 2}$$
	**推论 2.1** 若线性变换 $T$ 有 $n$ 个互异的特征值 $\lambda_1, \lambda_2, \cdots, \lambda_n$ ，则必有  $V_{\lambda_1} \oplus V_{\lambda_2} \oplus \cdots \oplus V_{\lambda_n} = V,$  从而 $T$ 有对角线上元素互异的对角矩阵表示．
	**推论 2.2** 线性变换 $T$ 有对角矩阵表示的充分必要条件是 $V$ 可分解成 $T$ 的一维不变子空间的直和．


## 2.2 Jordan矩阵介绍

### 1 Jordan矩阵

**定义 2.3** 形如  
$$J(\lambda) = \begin{pmatrix} \lambda & 1 & & \\ & \lambda & 1 & \\ & & \ddots & 1 \\ & & & \lambda \end{pmatrix} \tag{2-3}$$  
的 $r$ 阶方阵称为一个 $r$ 阶**Jordan块**，其中 $\lambda$ 是复数．由若干个Jordan块 $J_i(\lambda_i)$ 构成的准对角矩阵
 $$J(\lambda) = \begin{pmatrix} J_1(\lambda_1) & & & \\ & J_2(\lambda_2) & & \\ & & \ddots & \\ & & & J_m(\lambda_m) \end{pmatrix}, \tag{2-4}$$
称为**Jordan矩阵**．

Jordan矩阵是准对角形矩阵，呈上三角形，对角线上是它的全部特征值。Jordan块的特点是主对角线上元素相等，紧邻上方元素$a_{ii + 1} = 1$，其余元素为0。Jordan块都是一阶的Jordan矩阵就是对角矩阵。

**定理2.5**：设T是复数域上 $n$ 维线性空间V的一个线性变换，在 $V$ 中必定存在一组基，使 $T$ 在这组基下的矩阵是Jordan矩阵，并且这个Jordan矩阵除去其中Jordan块的排列次序外，是被 $T$ 惟一决定的，它称为 $T$ 的**Jordan标准形**.

定理 2.6 每个 $n$ 阶复数矩阵 $A$ 都与一个Jordan矩阵相似，这个Jordan矩阵除去其中Jordan块的排列次序外，是被矩阵 $A$ 惟一决定的，它称为 $A$ 的Jordan标准形。

### 2 用$\lambda$-矩阵构造Jordan标准形

若一个矩阵的元是 $\lambda$ 的多项式，就称为 $\lambda$-矩阵。

**定义 2.4** 下面的三种变换称为 $\lambda$ -矩阵的**初等变换**: 
1. 矩阵的两行(列)互换位置; 
2. 矩阵的某一行(列)乘以非零的常数$c$; 
3. 矩阵的某一行(列)加另一行(列)的$\varphi(\lambda)$倍，$\varphi(\lambda)$是一个多项式.

定义 2.5 如果λ-矩阵 $A(\lambda)$ 可以经过一系列初等变换化为 $B(\lambda)$ ，则称 $A(\lambda)$ 与 $B(\lambda)$ **等价**。

**定理 2.7** 任意一个非零的 $s \times n$ 的 $\lambda$ -矩阵 $A(\lambda)$ 都惟一地等价于**标准形**: $$ \left( \begin{array}{ccccccc} d_1(\lambda) & & & & & & \\ & d_2(\lambda) & & & & & \\ & & \ddots & & & & \\ & & & d_r(\lambda) & & & \\ & & & & 0 & & \\ & & & & & \ddots & \\ & & & & & & 0 \end{array} \right), $$ 其中 $r \geq 1$ ,  $d_i(\lambda), i = 1, 2, \cdots, r$ 是首项系数为1的多项式,且 $d_i(\lambda)$ 整除 $d_{i+1}(\lambda), i = 1, 2, \cdots, r - 1$ .

**定义 2.6** 标准形的主对角线上非零元 $d_1(\lambda), d_2(\lambda), \cdots, d_r(\lambda)$ 称为 $\lambda$ -矩阵的不变因子。
	矩阵 $A$ 的特征矩阵 $\lambda I - A$ 的不变因子以后就简称为 $A$ 的不变因子.

**定义 2.7** 把矩阵 $A$ (或线性变换 $T$ )的每个次数大于零的不变因子分解成互不相同的一次因式方幂的乘积, 所有这些一次因式方幂(相同的必须按出现的次数计算)称为矩阵 $A$ (或线性变换 $T$ )的**初等因子**.
**定理 2.8** 两个同阶矩阵相似的充分必要条件是它们有相同的初等因子.

---
下面的定理给出了一个求初等因子的方法，它不必事先知道不变因子。

**定理 2.9** 首先用初等变换化特征矩阵 $\lambda I - A$ 为对角形, 然后将主对角线上的元分解成互不相同的一次因式方幂的乘积, 则所有这些一次因式的方幂(相同的按出现的次数计算)就是 $A$ 的全部初等因子.

### 3 用分析确定法构造Jordan标准形

设 $A$ 的特征多项式为 $$|\lambda I - A| = (\lambda - \lambda_1)^{k_1} (\lambda - \lambda_2)^{k_2} \cdots (\lambda - \lambda_s)^{k_s},$$ 其中 $\lambda_i$ 是 $A$ 的 $k_i$ 重特征值.  $\sum_{i=1}^{s} k_i = n$ ,  $\lambda_1, \lambda_2, \cdots, \lambda_s$ 互异, 所以 
$$J(\lambda) = \left( \begin{array}{cccc} J_1(\lambda_1) & & & \\ & J_2(\lambda_2) & & \\ & & \ddots & \\ & & & J_s(\lambda_s) \end{array} \right), \tag{2-9}$$
 $J_i(\lambda_i)$ 是主对角线元素为 $\lambda_i$ 的 $k_i$ 阶Jordan矩阵. 它是同一特征值 $\lambda_i$ 对应的Jordan块放在一起得到的Jordan矩阵.

把相似变换的可逆矩阵 $P$ 依 $(2 - 9)$ 式所示 $J_A$ 的结构, 相应取 $k_1$ 列,  $k_2$ 列,  $\cdots$ ,  $k_s$ 列分块为 $P = (p_1\ p_2\ \cdots\ p_s)$ ,  $p_i \in \mathbb{C}^{n \times k_i}, i = 1, 2, \cdots, s$ ,  $AP = PJ_A$ 可表示为 
$$(Ap_1\ Ap_2\ \cdots\ Ap_s) = (p_1J_1(\lambda_1)\ p_2J_2(\lambda_2)\ \cdots\ p_sJ_s(\lambda_s))$$ 从而有 
$$Ap_i = p_iJ_i(\lambda_i), i = 1, 2, \cdots, s. \tag{2-10}$$

如何求可逆矩阵 $P$ :
1. 求$A$的特征多项式 
$$|\lambda I - A| = (\lambda - \lambda_1)^{k_1} (\lambda - \lambda_2)^{k_2} \cdots (\lambda - \lambda_s)^{k_s},$$
$\lambda_1, \lambda_2, \cdots, \lambda_s$ 互异, 则 $\lambda_i$ 为 $A$ 的 $k_i$ 重特征值, 其代数重数 $k_i$ 决定Jordan矩阵 $J_i(\lambda_i)$ 的阶数为 $k_i$ . 
2. 对于每一个 $\lambda_i$ , 由 $(A - \lambda_i I)X = 0$ , 求出 $A$ 的线性无关的特征向量 $\alpha_1, \alpha_2, \cdots, \alpha_{t_i}$ . 几何重数 $t_i$ 决定 $J_i(\lambda_i)$ 中有 $t_i$ 个Jordan块. 
3. 若 $k_i = t_i$ , 则 $\lambda_i$ 对应的Jordan矩阵为 $k_i$ 阶对角矩阵. 
	若 $t_i < k_i$ , 则选择适当的特征向量 $\alpha_i$ , 由 $(2 - 13)$ 式确定Jordan链的长度 $n_i$ , 从而得到 $J_A$ 的结构.
	上述矩阵等价于如下的方程组 $$ \left\{ \begin{aligned} (A - \lambda_1 I)\alpha_1 &= 0, \\ (A - \lambda_1 I)\beta_2 &= \alpha_1, \\ (A - \lambda_1 I)\beta_3 &= \beta_2, \\ &\vdots \\ (A - \lambda_1 I)\beta_{n_j} &= \beta_{n_j - 1}. \end{aligned} \right. \tag{2-13} $$ 从 $(2 - 13)$ 式求得一组向量 $\{\alpha_1, \beta_2, \cdots, \beta_{n_j}\}$ , 我们称之为Jordan链. 链中的第一个向量 $\alpha_1$ 是特征向量,  $\beta_2, \cdots, \beta_{n_j}$ 称为广义特征向量. 链的长度 $n_j$ 就是 $J_{1j}(\lambda_1)$ 的阶数.  $(2 - 13)$ 式给出一个递归过程, 该过程到线性方程组 $(A - \lambda_1 I)\beta_{n_j + 1} = \beta_{n_j}$ 无解时终止.
4. 所有Jordan链构成的矩阵 $P$ 就是Jordan标准形的变换矩阵.

## 2.3 最小多项式

一个方阵 $A_{n \times n}$ 的Jordan标准形 $J$ 在和矩阵 $A$ 相似的一切矩阵构成的相似类中，是形式最简单的。 $A$ 和 $J$ 又都具有相似不变性，因此，在讨论关于 $A$ 的一些问题时，常利用 $J$ 的简单形式，先就 $J$ 进行讨论，得到 $A$ 的相应结论，这就是Jordan化方法。

这一节用Jordan化方法先讨论 $A$ 的矩阵多项式的计算问题，由此证明Cayley定理，再给出矩阵 $A$ 的最小多项式及其性质。

### 1 矩阵多项式

**定义 2.8** 设 $A \in F^{n \times n}$ ， $a_i \in F$ ， $g(\lambda) = a_m \lambda^m + a_{m-1} \lambda^{m-1} + \cdots + a_1 \lambda + a_0$ 是一个多项式，则称矩阵 $g(A) = a_m A^m + a_{m-1} A^{m-1} + \cdots + a_1 A + a_0 I$ 为 $A$ 的矩阵多项式。

**定理 2.10** 设 $A \in F^{n \times n}$ ， $g(A)$ 是 $A$ 的矩阵多项式，则有如下结果： 
1. 若 $\lambda_0$ 是 $A$ 的特征值，则 $g(\lambda_0)$ 是 $g(A)$ 的特征值。 
2. 如果 $A$ 相似于 $B$ ，即 $P^{-1}AP = B$ ，则 $g(A)$ 相似于 $g(B)$ ，且 $P^{-1}g(A)P = g(B)$ 。
3. 如果 $A$ 为准对角矩阵，则 $g(A)$ 也是准对角矩阵。而且， 若 $A = \begin{pmatrix} A_1 &&& \\ & A_2 && \\ && \ddots & \\ &&& A_k \end{pmatrix}$ ， $A_i$ 为方子块，则 $g(A) = \begin{pmatrix} g(A_1) &&& \\ & g(A_2) && \\ && \ddots & \\ &&& g(A_k) \end{pmatrix}$ 。

---
下面, 我们讨论 $g(A)$ 的计算问题, 即已知方阵 $A$ 和多项式 $g(\lambda)$ , 如何计算方阵 $g(A)$ 的问题. 设 $g(\lambda) = a_m \lambda^m + a_{m-1} \lambda^{m-1} + \cdots + a_1 \lambda + a_0$ .

为计算 $g(A)$ ，用Jordan化方法. 设  
$$A = PJP^{-1}, \quad J = \begin{pmatrix} J_1(\lambda_1) &&& \\ & J_2(\lambda_2) && \\ && \ddots & \\ &&& J_s(\lambda_s) \end{pmatrix},$$
其中 $J_1(\lambda_1), J_2(\lambda_2), \cdots, J_s(\lambda_s)$ 为 $J$ 的全部Jordan块. 由定理2.10知  
$$\begin{align*} g(A) &= Pg(J)P^{-1} \\ &= P \begin{pmatrix} g(J_1(\lambda_1)) &&& \\ & g(J_2(\lambda_2)) && \\ && \ddots & \\ &&& g(J_s(\lambda_s)) \end{pmatrix} P^{-1}, \end{align*}$$ 因此计算 $g(A)$ 的问题转化为对Jordan块 $J_i(\lambda_i)$ 计算 $g(J_i(\lambda_i))$ 的问题.

[[矩阵论例子#g(A)分析]] ^3dfa00

### 2 方阵的化零多项式

对 $n$ 阶方阵 $A$ ，若存在多项式 $g(\lambda)$ ，使矩阵 $g(A) = 0$ ，则称 $g(\lambda)$ 为矩阵 $A$ 的**化零多项式**。 

**定理 2.11 (Cayley-Hamilton)** 设 $A \in F^{n \times n}$ ，则方阵 $A$ 的特征多项式就是 $A$ 的化零多项式。

[[矩阵论例子#定理2.11证明]] ^b68863

### 3 最小多项式

对 $A \in F^{n \times n}$ ，Cayley定理指出它有化零多项式，

事实上， $A$ 有很多化零多项式. 若 $A$ 是 $V$ 上线性变换 $T$ 的矩阵，对 $A$ 的化零多项式 $f(\lambda)$ ，相应地用多项式得到线性变换 $f(T)$ 就是零变换. 所以线性变换 $T$ 的化零多项式也是存在的，而且  
$$\text{线性变换 } f(T) = \mathbf{0} \iff \text{矩阵 } f(A) = 0.$$

---
**定义 2.9** 设 $T$ 是线性空间 $V$ 上的线性变换， $m_T(\lambda)$ 是一个关于文字 $\lambda$ 的多项式，如果 $m_T(\lambda)$ 满足：
1.  $m_T(\lambda)$ 的最高次项系数为1， 
2. $m_T(\lambda)$ 是 $T$ 的一个化零多项式，即 $m_T(T) = \mathbf{0}$ ，
3.  $m_T(\lambda)$ 是 $T$ 的化零多项式中次数最低的多项式， 
则称 $m_T(\lambda)$ 是 $T$ 的**最小多项式**。

**定理 2.12**  $T$ 的特征多项式 $f(\lambda)$ 与最小多项式 $m_T(\lambda)$ 有相同的根(重数不计). 即若 $f(\lambda) = |\lambda I - A| = (\lambda - \lambda_1)^{r_1} (\lambda - \lambda_2)^{r_2} \cdots (\lambda - \lambda_s)^{r_s}$ ，则 $m_T(\lambda) = (\lambda - \lambda_1)^{t_1} (\lambda - \lambda_2)^{t_2} \cdots (\lambda - \lambda_s)^{t_s}$ ， $1 \leq t_i \leq r_i$ ， $i = 1, 2, \cdots, s$ .

**定理 2.13** 设变换 $T$ 的特征多项式为  $$f(\lambda) = (\lambda - \lambda_1)^{r_1} (\lambda - \lambda_2)^{r_2} \cdots (\lambda - \lambda_s)^{r_s},$$  又 $T$ 的Jordan标准形中关于特征值 $\lambda_i$ 的Jordan块的最高阶数为 $\bar{n}_i$ ，则 $T$ 的最小多项式  $$m_T(\lambda) = (\lambda - \lambda_1)^{\bar{n}_1} (\lambda - \lambda_2)^{\bar{n}_2} \cdots (\lambda - \lambda_s)^{\bar{n}_s}.$$ (证明没细看)

**定理 2.14** 线性变换 $T$ 可对角化的充分必要条件是 $T$ 的最小多项式 $m_T(\lambda)$ 是一次因子的乘积.

[[最小多项式和jordan标准型]] ^b6e7f8

