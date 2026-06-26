# 输出位深相关的图像低阶信号统计：信噪比、量化相位与近黑行为

Low-Order Signal Statistics Related to Image Output Bit Depth: Signal-to-Noise Ratio, Quantization Phase, and Near-Black Behavior

作者：姜尧耕（<https://y-g-jiang.github.io>）
日期：2026-06

Update:

https://github.com/y-g-jiang/paper-LOSSRIO

## 摘要

本文研究线性 RAW 输出位深（16/14/12 bit）对信号还原能力的影响，并将其量化为低位深输出与理想线性信号在信噪比、条件均值、误差能量与空间结构上的相似性。近黑 banding 被建模为带高斯读噪的最近邻量化问题，核心对象是输出的条件分布及其各阶条件矩。量化输出的特征函数等于输入特征函数以量化频率为周期的复制求和，因此一切与量化网格相干的统计量，其各阶傅里叶谐波都被输入特征函数在相应格点处的取值加权；对高斯读噪，该权重随读噪与量化步长之比按指数衰减，读噪由此充当一个高效的非减性 dither，将与相位相干的确定性结构指数压低。在此框架下，本文给出一阶条件均值偏差的闭式傅里叶级数并证明其指数衰减，进而推导二阶误差能量的完整傅里叶形式。对大量相机数据进行分析，比较其 14-bit、12-bit 下的信噪比与动态范围差异，并给出相应的低阶统计判据；对富士 GFX100S 的近黑实测读噪（读噪与量化步长之比约 1.18），16→14 bit 的一阶确定性偏置约 10⁻¹² DN16、二阶相位调制约 10⁻¹⁰ DN16²，均远低于常规数值分辨能力，相应的信噪比代价亦可忽略。因此，原始 16→14 bit 量化本身不足以保留可测的低阶量化相位结构，14 bit 与 16 bit 输出对理想信号已充分相似。全部分析限定于线性 RAW 域、最近邻量化与独立同分布高斯读噪模型。

关键词：RAW 位深；16-bit；14-bit；12-bit；PDR；量化噪声；dither；banding；GFX100S；条件矩；傅里叶谐波；特征函数；Sheppard 修正

---

## 1. 引言

编码精度的意义建立在其承载的信号信息上。如果传感器读噪、shot noise、PRNU/DSNU、黑场校正误差和后续链路已经让最低几位变成随机摆动，此时若直接截位到低位，数值上的空格不必然对应稳定的量化相位结构，极端的可以参考类似二值半调的形式。

量化误差是输入的确定性函数，加性独立噪声只在特定条件下对矩成立 [8]。所以本文给出 $s$、读噪 $N$，研究条件分布 $P(Y_b\mid s)$ 和它的各阶条件矩，如此刻画读噪把确定性结构打散到什么程度。

设理想线性 RAW 信号为 $s$，输出位深为 $b$，量化输出为 $Y_b$。我们关心的不是码值本身，而是输出相对理想信号的误差

$$
e_b(s)=Y_b-s .
$$

如果两个位深 $b_1,b_2$ 的误差分布在几个关键统计量上足够接近，也就是

$$
P(Y_{b_1}\mid s)\approx P(Y_{b_2}\mid s),
$$

那么对这个信号来说，$b_1$ 和 $b_2$ 就是相似的，在如下四层：

1. PDR 或目标 SNR 的损失够小；
2. 条件均值偏差 $B_b(s)=\mathbb E[Y_b-s\mid s]$ 够小；
3. 二阶误差能量 $M_{2,b}(s)=\mathbb E[(Y_b-s)^2\mid s]$ 的量化相位调制够小；
4. 在视觉或空间处理的尺度上，残差的相关性、频谱和局部包络没有成稳定结构。

第一层关心损失多少信噪比，第二层关心无限堆栈之后会不会留下均值结构，第三层关心零均值残差的颗粒强弱还会不会随量化相位起伏，第四层关心这些起伏会不会在空间上排成稳定图样。本文把条件均值、条件二阶矩等低阶条件统计量中随 $s\bmod\Delta$ 周期变化的部分，称为低阶量化相位结构；其中一阶结构对应确定性均值偏置，二阶结构对应零均值残差的能量调制。

---

## 2. 经典框架：dither与量化的特征函数

### 2.1 量化输出的特征函数

设最近邻量化器

$$
Q_\Delta(x)=\Delta\operatorname{round}\!\left(\frac{x}{\Delta}\right),
$$

作用在随机变量 $X$ 上。Widrow 的统计量化理论给出量化输出的特征函数等于输入特征函数以 $\Psi=2\pi/\Delta$ 为周期复制求和 [4][5][6]。

$$
\Phi_{Q_\Delta(X)}(u)=\sum_{l=-\infty}^{\infty}\Phi_X(u+l\Psi)\,\operatorname{sinc}\!\Big(\tfrac{\Delta(u+l\Psi)}{2}\Big).
$$

因此，量化误差里凡是随 $s$ 以量化网格为周期、也就是和相位相干的统计量，它的第 $k$ 个傅里叶谐波系数都正比于输入特征函数在频点 $2\pi k/\Delta$ 处的值 $\Phi_X(2\pi k/\Delta)$。当 $X=s+N$、读噪 $N$ 与 $s$ 独立时，

$$
\Phi_X\!\left(\tfrac{2\pi k}{\Delta}\right)=e^{\,i 2\pi k s/\Delta}\,\Phi_N\!\left(\tfrac{2\pi k}{\Delta}\right).
$$

Sripad 和 Snyder 证明，量化误差严格均匀且白的充要条件，正是这些格点值全为零，即 $\Phi_X(2\pi l/\Delta)=0,\ \forall l\neq0$ [7]。一般输入达不到严格为零，但只要这些格点值足够小，量化误差就近似是均匀白噪声 [6][7][9][11]。下文中即为，对高斯读噪，它可化为一个关于读噪与量化步长之比 $r=\sigma/\Delta$ 的容差判据。最坏相位下量化噪声方差对 $\Delta^2/12$ 的主导偏离低于给定容差 $\varepsilon$ 时，所需 $r$ 只随 $\varepsilon$ 对数级增长（§3.1 给出 $\delta(r)=(4r^2+1/\pi^2)e^{-2\pi^2r^2}$ 及其近似阈值，§7 给出相应的二阶相位调制量级）。

### 2.2 高斯读出噪声的dithering行为

对零均值高斯读噪 $N\sim\mathcal N(0,\sigma^2)$，特征函数处处不为零，

$$
\Phi_N\!\left(\tfrac{2\pi k}{\Delta}\right)=\exp\!\left[-\tfrac12\Big(\tfrac{2\pi k}{\Delta}\Big)^2\sigma^2\right]=\exp\!\left[-2\pi^2k^2\Big(\tfrac{\sigma}{\Delta}\Big)^2\right].
$$

这个阻尼因子说明，高斯读噪没法把任何一阶矩调制精确清零（高斯特征函数没有实零点），但能把第 $k$ 个相位相干谐波指数十分有效地降低。记 $r=\sigma/\Delta$，主谐波 $k=1$ 被压低 $e^{-2\pi^2r^2}$。

Bennett 在 1948 年就用 Rice 的特征函数方法，对高斯输入的均匀量化器误差算出了精确谱，其误差相关函数 [2]

$$
G(a)=\frac{k}{2\pi^2}\sum_{n=1}^{\infty}\frac{1}{n^2}\exp\!\Big[-\frac{4\pi^2 n^2(1-a)}{k}\Big],\qquad k=\frac{\Delta^2}{\sigma^2},
$$

每个谐波都带着同型的高斯因子，和本文的阻尼因子同源。一般地，输入分布越平滑、格点特征函数越小，量化误差越接近满足量化定理 [13]；在本文采用的高斯读噪模型中，特征函数本身就是高斯型，因而在格点 $2\pi k/\Delta$ 处给出非常强的指数衰减。所以本文的 $r=\sigma/\Delta$ 指数压低，本质就是高斯特征函数在这些格点处的残值。

---

## 3. 量化步长、均匀量化噪声与 PDR 损失

### 3.1 步长与均匀量化噪声

统一用 16-bit DN（记 DN16）做单位，16-bit 最小步长 $\Delta_{16}=1\,\mathrm{DN16}$，$16\!\to\!14$-bit 的等效步长 $\Delta_{14}=4\,\mathrm{DN16}$。若改用 14-bit ADU 比较 14-bit 与更低位深 $b$，量化步长为 $\Delta_b=2^{14-b}$。当量化误差与输入低协方差、量化相位又被噪声打散时，均匀量化噪声近似给出经典结果：

$$
\sigma_{q,b}^2=\frac{\Delta_b^2}{12},\qquad \sigma_{q,b}=\frac{\Delta_b}{\sqrt{12}} .
$$

这只是信噪比层面的低协方差近似。$\Delta^2/12$ 只在量化相位被噪声充分随机化时才成立。对高斯输入，完整推导见我之前的文章 [27]，这个近似的成立条件可以写成无量纲比 $r=\sigma/\Delta$ 的容差判据：最坏相位下量化噪声方差对 $\Delta^2/12$ 的主导偏离为 $\delta(r)=(4r^2+1/\pi^2)e^{-2\pi^2r^2}$，要它低于容差 $\varepsilon$ 时 $r$ 只需对数级增长；若忽略括号中的常数修正，可得到常用的 Lambert W 近似闭式 $r\gtrsim\sqrt{-W_{-1}(-\pi^2\varepsilon/2)/(2\pi^2)}$。常被引用的 $r\approx0.5$ 并不是一个绝对门槛，而是对应百分之一量级容差的取值：$r=0.5$ 给出主导方差偏离 $\delta\approx7.9\times10^{-3}\,\Delta^2$，要更严的拟合就按上式取更大的 $r$。与本文 §7 二阶矩中相位偏离按 $e^{-2\pi^2k^2(\sigma/\Delta)^2}$ 衰减的结论一致。平滑渐变里的相干结构得回到 §6 至 §8 的条件分布和相位分析。

![图 1. 16-bit 与 14-bit 的 PTC fit DR 对照，用来估计 16→14-bit 对近黑噪声和 PDR 的实际影响。](figures/bitdepth_original_01.png)

因而我们可以快速算出，对 $16\!\to\!14$-bit 的 GFX100S，按近黑目标 SNR 极限估计，最大动态范围差异约
$$
\Delta DR_{\max}=\log_2\frac{4.8766}{4.7378}=0.0417\ \mathrm{stops}.
$$

也就是说只要如此小的曝光补偿，GFX100S 后采样到 14-bit下的结果就能优于 16-bit 的噪声水平，因此可以认为是在 SNR 意义上近乎无差别的。

### 3.2 目标 SNR 与 PDR 损失

设信号电子数 $S$、总读出等效噪声 $N$（电子单位，shot 方差取 $S$），则

$$
\mathrm{SNR}=\frac{S}{\sqrt{S+N^2}} .
$$

给定目标 SNR $T$，令 $S/\sqrt{S+N^2}=T$，平方后得 $S^2-T^2S-T^2N^2=0$，取正根

$$
S_T(N)=\frac{T^2+\sqrt{T^4+4T^2N^2}}{2}.
$$

若 14-bit 与低位深 $b$ 的噪声分别是 $N_{14}(I),N_b(I)$，则 $S_{14}=S_T(N_{14})$、$S_b=S_T(N_b)$，PDR 损失（单位 stop）为

$$
\Delta\mathrm{PDR}_b(I)=\log_2\!\left(\frac{S_b(I)}{S_{14}(I)}\right),
$$

含义就是低位深为了达到同一目标 SNR，需要多给多少曝光。如我在 [27] 中证明的，当读噪达到约半个量化步长（$r=\sigma/\Delta\gtrsim0.5$）时，量化噪声相对理想 $\Delta^2/12$ 的偏离已经很小。在最坏相位下主导方差偏离为 $\delta(0.5)=(1+1/\pi^2)e^{-\pi^2/2}\approx7.9\times10^{-3}\,\Delta^2$，约为理想量化噪声方差 $\Delta^2/12$ 的 9.5%（判据见 §2 与 [27]）。因此可以把低位深量化噪声近似按正交方式并进读噪：

$$
N_b(I)=k_{14}(I)\sqrt{R_{14}(I)^2+\frac{2^{2(14-b)}}{12}},
$$

其中 $R_{14}(I)$ 是 14-bit ADU 单位的读噪，$k_{14}(I)$ 是从 ADU 折算到电子或统一噪声单位的比例因子，也就是转换增益的倒数（见 §4.2）。

![图 2. 14→12-bit 的 ΔPDR 临界 ISO 曲线，多机型在不同 ISO 下的 PDR 损失趋势。](figures/bitdepth_original_02.png)

![图 3. 多机型 14→12-bit PDR 损失曲线，在 0.1 stop 阈值附近的横向分布。](figures/bitdepth_original_03.png)

把图 2、图 3 的多机型趋势落到实拍，可得两条更直接的结论。

其一，在 $12$-bit 模式比 $14$-bit 多得约 $0.2$ 档曝光的条件下，对 ISO 800 以上、除 A7S3 外的一切索尼机型（截至 2025 年；索尼已是当下噪声表现最好的一批，故这大概率覆盖一切相机），$12$-bit 在信噪比上已完全领先 $14$-bit 模式。

其二，更一般地，ISO 400 以上的一切市售相机（截至 2025 年），$14\!\to\!12$-bit 的 PDR 损失都低于 $0.1$ 档。

这里"完全领先"是因为：降位深的 SNR 损失在噪声底最重，越往亮部、信号越强，量化噪声相对信号越小、损失越轻，而 PDR 损失正是这条 SNR–EV 曲线在噪声底处取到的极值。于是只要用这 $0.2$ 档曝光把 $12$-bit 的噪声底补到与 $14$-bit 齐平，就等于补平了损失最重的那一点；其余更亮的 EV 处损失本就更小，同样的曝光增益便会富余，于是 $12$-bit 在全线 SNR–EV 上都占优。

在高模拟增益区，12-bit 额外那点量化噪声通常已经被读噪和 photon noise 淹没，代价就更为细微了。

### 3.3 像素数与位深损失

固定传感器面积，固定单位面积满阱容与单位面积读噪功率，只改变像素数 $P$，即得位深损失随 $P$ 的标度。

取目标 SNR $\to0$ 极限，§3.2 的位深损失化为噪声底之比。

$$
\Delta\mathrm{PDR}=\log_2\frac{\sqrt{R^2+\Delta^2/12}}{R}=\tfrac12\log_2\!\Big(1+\frac{1}{12\,r^2}\Big),\qquad r=\frac{R}{\Delta},
$$

$R$ 为读噪、$\Delta$ 为量化步长（同一单位），故损失仅由 $r=\sigma/\Delta$ 决定。各量随 $P$ 的标度为：

- 单位面积满阱一定 $\Rightarrow$ 每像素满阱 $\mathrm{FW}\propto1/P$；
- 定位深、ADC 满量程对满阱（PLT 模型 [24] 的增益约定 $k=\mathrm{FW}/2^{14}$）$\Rightarrow$ 步长 $\Delta\propto\mathrm{FW}\propto1/P$；
- 单位面积读噪功率一定 $\Rightarrow P\,R^2=\mathrm{const}$，即 $R\propto1/\sqrt P$。

代入得

$$
r=\frac{R}{\Delta}\propto\sqrt P,\qquad \Delta\mathrm{PDR}=\tfrac12\log_2\!\Big(1+\frac{c}{P}\Big),
$$

$c$ 为常数；$\Delta\propto1/P$ 随满阱增大快于 $R\propto1/\sqrt P$，故 $P$ 越大则 $r$ 越大、损失越小。

代入主页中PLT站公式来直观展示（全画幅量级示意参数，单位面积满阱与读噪功率固定，14-bit，基准 ISO）：

| 像素数 (MP) | 满阱 (e⁻) | 读噪 (e⁻) | $r=R/\Delta$ | $\Delta\mathrm{PDR}_{14\to12}$（极限） | $\Delta\mathrm{PDR}_{14\to12}$（$T=20$） |
|---|---|---|---|---|---|
| 100 | 8000 | 1.50 | 0.77 | 0.095 | 0.0011 |
| 50 | 16000 | 2.12 | 0.54 | 0.180 | 0.0041 |
| 25 | 32000 | 3.00 | 0.38 | 0.323 | 0.016 |
| 12.5 | 64000 | 4.24 | 0.27 | 0.546 | 0.057 |

这个标度可依赖下两种约定：定位深、满阱匹配时 $\Delta\propto1/P$，单位面积量化功率 $\propto1/P$；固定增益时（$\Delta$ 不变、位深随满阱增大）单位面积量化功率 $\propto P$。本文与 PLT 取前者。

满阱本身给出与量化无关的另一标度：$\mathrm{FW}\propto1/P$、$R\propto1/\sqrt P$，每像素动态范围 $\mathrm{FW}/R\propto1/\sqrt P$，每减半像素增 $\tfrac12$ stop（表中 $100\to25$ MP，$\mathrm{FW}/R$ 由 $12.4$ 升至 $13.4$ stop）。

---

## 4. banding风险位置

### 4.1 决定比值 $r=\sigma/\Delta$

banding 由量化步长 $\Delta$、读噪标准差 $\sigma$、空间渐变斜率、显示映射和观察尺度一起决定。对本文的信号统计问题，最核心的无量纲比值是 $r=\sigma/\Delta$。$r$ 很小时，输入在量化阈值附近缺少随机抖动，平滑渐变就会和量化网格相干；$r$ 够大时，量化相位被打散，一阶均值结构按 $e^{-2\pi^2r^2}$ 迅速衰减（§2）。

由光子转移特性可知，同一位深下量化步长 $\Delta$ 不随信号改变，而读噪标准差 $\sigma$ 随信号单调上升（§4.2），故 $r=\sigma/\Delta$ 在低码值处最小、向亮部单调增大。这使 $r$ 成为一个可跨亮度比较的保守指标：在亮度、空间频率、观察条件、tone curve 局部导数与 CSF 响应都相同的前提下，若两处的 $r$ 相等，则量化台阶的可见性只随 $\Delta$ 单调增大——步长越大、台阶越粗，越容易越过视觉阈值。而暗部最低码值处 $r$ 最小，且将其提亮到与亮部相同的目标感知范围后 $r$ 仍保持不变（线性提亮对信号与噪声同步缩放，不改变 $\sigma/\Delta$）。因此暗部最低码值处的表现构成整机断层性能最保守的瓶颈：在上述前提近似一致时，暗部最低码值处给出较保守的量化相位风险估计。

![图 4. PTC 曲线示例。随信号升高噪声标准差增大，所以低码值处往往给出最保守、也就是最小的 σ/Δ。](figures/bitdepth_original_04.png)

### 4.2 近黑噪声模型

GFX100S 的近黑噪声标准差以线性信号码值 $m$ 为自变量，由实测 PTC 拟合得到下式：

$$
\sigma_{\mathrm{DN}}(m)=\sqrt{4.729^2+\frac{m}{0.71}+(0.007\,m)^2}\ \ \mathrm{DN16}
$$

$m$ 是线性信号码值。即为传感器噪声的标准三项结构 $\sigma^2=\sigma_{\text{read}}^2+S/g+(\mathrm{PRNU}\cdot S)^2$（Janesick 式 3.16–3.17，EMVA 1288 线性模型）[15][16][17]，三项一一对应：

- 常数项 $4.729^2=22.36$ 是读噪方差（量化项 $\Delta^2/12$ 已正交并进去）；
- 线性项 $m/0.71=1.408\,m$ 是 shot 噪声，系数是转换增益的倒数，$K=0.71\ \mathrm{e^-/DN16}$（等价 $1.408\ \mathrm{DN16/e^-}$）；
- 二次项 $(0.007\,m)^2=4.9\times10^{-5}\,m^2$ 是 $(\mathrm{PRNU}\cdot m)^2$，$\mathrm{PRNU}=0.7\%$，落在 Janesick 给的约 1% 的典型范围里。

因为此结果是我对dpr标版采样得到的（PTP官网测量错误），PRNU可能略高，但PRNU不影响极暗部噪声行为。

由此读噪 $4.729\,\mathrm{DN16}\times0.71\,\mathrm{e^-/DN16}\approx3.36\ \mathrm{e^-}$，近黑 $m=0$ 时

$$
\sigma_{\mathrm{DN}}(0)=4.729\,\mathrm{DN16},\qquad r=\frac{\sigma}{\Delta}=\frac{4.729}{4}\approx1.18 .
$$

这就是 GFX100S 近黑的工作比值，本文结论就建立在这个实测读噪上。图 6 即可见，在相同的码值跨度上，高码值区域的总噪声更大，因此同一量化步长对应的低阶相位结构更容易被随机项淹没。

![图 5. GFX100S 相关的 PTC 拟合图。近黑读噪模型用来给出最保守（最小）的 σ/Δ。](figures/bitdepth_original_05.png)

![图 6. 不同码值范围里的噪声与渐变对照。高码值区自然噪声更强，低码值区才是 banding 的最坏测试点。](figures/bitdepth_original_06.png)

---

## 5. 后期拉伸与非线性变换

设后期操作在局部可微，记为 $g$，对小误差 $e$ 有 Taylor 展开

$$
g(s+e)=g(s)+g'(s)\,e+\tfrac12 g''(s)\,e^2+O(e^3).
$$

一阶项说明，旧读噪和旧量化残差被同一个局部导数 $g'(s)$ 放大，所以比值 $r=\sigma/\Delta$ 在线性缩放、局部白平衡、通道增益、密度线性变换这些操作下一阶不变——决定低阶量化相位结构的不是绝对码值差，而是这个被同步缩放的比值。这和显示侧的结论一致：决定 banding 风险的更多是位的分配方式（编码曲线、gamma、LUT），而不是位数本身 [23]。

二阶项 $\tfrac12 g''(s)\,\mathbb E[e^2\mid s]$ 中，唯有 $\mathbb E[e^2\mid s]$ 随 $s$ 相位相干变化的那部分才可能传播为新的确定性结构；而对原始 $16\!\to\!14$-bit 量化误差，这一相位相干调制已由 §7 算定，近黑 $r\approx1.18$ 处仅约 $9.5\times10^{-11}\,\mathrm{DN16}^2$。对常见的平滑局部变换而言，这一极小二阶调制难以被放大到主导信号统计的程度。

---

## 6. 一阶条件均值偏差：确定性均值结构

### 6.1 模型及期望

设量化前 $X=s+N,\ N\sim\mathcal N(0,\sigma^2)$，量化输出 $Y=Q_\Delta(X)$，误差 $e=Y-s$。固定 $s$ 反复拍、在独立同分布读噪的假设下做无限平均：

$$
\lim_{n\to\infty}\frac1n\sum_{i=1}^n Q_\Delta(s+N_i)=\mathbb E[Q_\Delta(s+N)\mid s].
$$

于是真正在无限堆栈之后还留得下来的，是条件均值偏差

$$
B_\Delta(s)=\mathbb E[Q_\Delta(s+N)\mid s]-s .
$$

如果 $B_\Delta(s)$ 随 $s\bmod\Delta$ 周期变化，就有确定性的均值相位结构；如果它够小，则在条件均值意义上不留下可测的量化相位结构。

### 6.2 闭式傅里叶级数

最近邻量化残差 $q(u)=\operatorname{round}(u)-u$ 是周期函数，标准锯齿展开为

$$
q(u)=\sum_{k=1}^{\infty}\frac{(-1)^k}{\pi k}\sin(2\pi k u).
$$

按 §2 的基本原理对 $X=s+N$ 逐项取期望，第 $k$ 谐波乘上 $\Phi_N(2\pi k/\Delta)=\exp(-2\pi^2k^2r^2)$，得到

$$
B_\Delta(s)=\frac{\Delta}{\pi}\sum_{k=1}^{\infty}\frac{(-1)^k}{k}\exp\!\left[-2\pi^2k^2\Big(\frac{\sigma}{\Delta}\Big)^2\right]\sin\!\Big(2\pi k\frac{s}{\Delta}\Big).
$$

$\sin(2\pi k s/\Delta)$ 是量化相位的周期结构，指数项是高斯读噪对它的抹平。

### 6.3 严格上界及阈值对判据的对数级不敏感

取首谐波得到主导项 $\dfrac{\Delta}{\pi}e^{-2\pi^2r^2}$。用 $1/k\le1$ 和 $k^2\ge1+3(k-1)$ 可以得到严格上界

$$
|B_\Delta(s)|\le\frac{\Delta}{\pi}\sum_{k\ge1}\frac1k e^{-2\pi^2k^2r^2}\le\frac{\Delta}{\pi}\,\frac{e^{-2\pi^2r^2}}{1-e^{-6\pi^2r^2}},
$$

其中修正因子 $1/(1-e^{-6\pi^2r^2})$ 在 $r\gtrsim0.4$ 时与 $1$ 的差小于 $10^{-4}$，故工程上以首谐波作为上界即可。

把 16-bit 写成 $Q_1$、$16\!\to\!14$-bit 写成 $Q_4$（都用 DN16），两者期望差 $B_D(s)=B_4(s)-B_1(s)$，由三角不等式 $|B_D|\le|B_4|+|B_1|$。$\Delta=1$ 的 16-bit 项衰减极快（$r_{16}=\sigma$），主导风险来自 14-bit 的 $\Delta=4$ 项（$r_{14}=\sigma/4$）：

$$
|B_D(s)|\lesssim\frac{4}{\pi}\exp\!\left[-\frac{\pi^2\sigma^2}{8}\right].
$$

设幅度判据是确定性偏置 $<\varepsilon\,\mathrm{DN16}$，解 $\tfrac4\pi e^{-\pi^2\sigma^2/8}=\varepsilon$ 得

$$
\sigma(\varepsilon)=\frac{2\sqrt2}{\pi}\sqrt{\ln\!\Big(\frac{4}{\pi\varepsilon}\Big)} .
$$

 $\sigma(\varepsilon)$ 对判据 $\varepsilon$ 仅有对数级依赖，因而极为稳健：

| 判据 $\varepsilon$ (DN16) | $10^{-2}$ | $10^{-4}$ | $10^{-6}$ | $10^{-12}$ |
|---|---|---|---|---|
| 所需 $\sigma$ (DN16) | $1.98$ | $2.77$ | $3.38$ | $4.75$ |

把判据收紧 $4$ 个数量级（$10^{-2}$ 到 $10^{-6}$），所需读噪也才从 $1.98$ 升到 $3.38\,\mathrm{DN16}$。而 GFX100S 近黑 $\sigma=4.729\,\mathrm{DN16}$，代进去得确定性偏置

$$
|B_D|\approx\frac{4}{\pi}e^{-\pi^2(4.729)^2/8}\approx1.3\times10^{-12}\,\mathrm{DN16},
$$

落在 $\varepsilon=10^{-12}$ 一档——比一个 16-bit 码低约 $12$ 个数量级。因此 $16\!\to\!14$-bit 不留下可测的一阶相位结构这一结论，对阈值的选取几乎不敏感；以阈值的对数稳健性为依据比给一个高有效位的阈值要好。后者（例如 $\sigma=1.982073445$）只是某个近似不等式在 $\varepsilon=0.01$ 处的根，模型与判据都不足以支撑如此多的有效位。

![图 7. 不同 RNDN 下 14-bit 与 16-bit 的 RMS 差异随平均帧数变化。一阶相干偏置被压住之后，堆栈效率主要由随机噪声决定。](figures/bitdepth_original_07.png)

---

## 7. 二阶矩分布

一阶均值偏差为零不代表误差分布没有结构。一个零均值随机变量照样可以有随输入变化的方差、尾部、偏斜、峰度和空间相关性。对低阶量化相位结构来说，第二个最基础的量是二阶误差能量 $M_{2,\Delta}(s)=\mathbb E[(Y-s)^2\mid s]$。

### 7.1 三项分解

归一化 $\theta=s/\Delta,\ r=\sigma/\Delta,\ Z=N/\Delta\sim\mathcal N(0,r^2)$，归一化总误差 $\eta_\theta=(Y-s)/\Delta=\operatorname{round}(\theta+Z)-\theta=Z+q(\theta+Z)$。于是

$$
M_2(\theta)=\mathbb E[\eta_\theta^2]=\underbrace{\mathbb E[Z^2]}_{\text{读噪能量}}+\underbrace{2\,\mathbb E[Z\,q(\theta+Z)]}_{\text{交叉项}}+\underbrace{\mathbb E[q(\theta+Z)^2]}_{\text{量化残差能量}} .
$$

该分解若直接写成读噪方差加 $\Delta^2/12$，便默认交叉项与相位项可以略去；这在 PDR 近似中通常足够，但在二阶相位分析中需要显式保留这些项（见 §7.3）。

### 7.2 傅里叶闭式

量化残差平方的周期展开为

$$
q(u)^2=\frac{1}{12}+\sum_{k=1}^{\infty}\frac{(-1)^k}{\pi^2 k^2}\cos(2\pi k u).
$$

两个高斯期望（第二个可由 Stein 引理 $\mathbb E[Zg(Z)]=r^2\mathbb E[g'(Z)]$ 得到）

$$
\mathbb E[\cos(2\pi k(\theta+Z))]=e^{-2\pi^2k^2r^2}\cos(2\pi k\theta),
$$
$$
\mathbb E[Z\sin(2\pi k(\theta+Z))]=2\pi k r^2 e^{-2\pi^2k^2r^2}\cos(2\pi k\theta).
$$

代回三项分解（$q$ 和 $q^2$ 的傅里叶系数分别是 $O(1/k)$、$O(1/k^2)$，乘上 $|\sin|,|\cos|\le1$ 再经高斯阻尼后级数绝对收敛，由控制收敛定理保证），得到

$$
M_2(\theta)=r^2+\frac{1}{12}+\sum_{k=1}^{\infty}(-1)^k\Big(4r^2+\frac{1}{\pi^2k^2}\Big)e^{-2\pi^2k^2r^2}\cos(2\pi k\theta).
$$

乘回 DN16 单位即 $\mathbb E[(Q_\Delta(s+N)-s)^2\mid s]=\Delta^2 M_2(s/\Delta)$。你可以见到，$r\to0$ 时 $M_2\to q(\theta)^2$（纯量化误差），$r\to\infty$ 时 $M_2\to r^2+\tfrac1{12}$，回到我原文章中的理想状况。

### 7.3 量级及交叉项

二阶相位谐波同样被指数降低。在 GFX100S 近黑 $r\approx1.18$：首谐波阻尼 $e^{-2\pi^2(1.18225)^2}\approx1.04\times10^{-12}$，二阶主相位半幅 $(4r^2+1/\pi^2)e^{-2\pi^2r^2}\approx5.9\times10^{-12}$，乘回 $\Delta^2=16$ 后约 $9.5\times10^{-11}\,\mathrm{DN16}^2$。

交叉项只在小 $r$ 时才显著，而非在本文的近黑工作点。在 $r=0.5$ 处，省掉交叉项会造成约 $0.00719\,\Delta^2$ 的二阶矩误差（首谐波 $4r^2e^{-2\pi^2r^2}$），这时读噪方差加 $\Delta^2/12$ 确实不够。但在 $r\approx1.18$ 处，交叉项和全部谐波都只有 $10^{-11}$ 到 $10^{-12}$ 量级，使得 $\sigma^2+\Delta^2/12$ 精确到约 11 至 12 位。所以保留交叉项，是为了验算数值的需要；而在近黑工作点本身，简化式已经够用。

![图 8. 高斯读噪下的一阶均值偏差和二阶噪声能量调制阈值。](figures/gaussian_noise_moment_modulation_thresholds.png)

![图 9. 实用范围放大：σ/Δ=0.5 已经很小，GFX100S 近黑处更低。](figures/gaussian_practical_moment_modulation_zoom.png)

## 8. 从条件矩到低阶量化相位结构

前两节给出的是条件矩。但 banding 本质上是空间相干现象，因此这里说明逐点条件矩如何约束低阶量化相位结构。

稳定条带至少要求误差里存在随 $s$ 确定性变化、也就是相位相干的分量：确定性依赖 $s$ 的成分才会在相邻像素间排成稳定的等级边界，并在重复实现、多帧平均或空间积分里保留下来；和 $s$ 不相干的部分只是随机噪声，空间上不会攒成条带。由 §2 的基本原理，每一阶条件矩里随 $s$ 周期、相位相干的那部分，第 $k$ 谐波都被乘上 $e^{-2\pi^2k^2r^2}$。在 GFX100S 近黑 $r\approx1.18$，连 $k=1$ 都只有 $e^{-2\pi^2r^2}\approx1.0\times10^{-12}$。也就是说，本文所分析的一阶均值结构与二阶能量调制都被压到极低；在独立读噪假设下，剩余项主要表现为能量接近 $\sigma^2+\Delta^2/12$ 的随机残差。

---

## 9. 关于可视化展示

banding 是视觉和信号统计一起决定的。折线在码值上上下跳，不等于就有 banding；若把它归因于量化相位结构，至少应能指出某个空间位置或灰阶区间出现了稳定的等级边界，而且这条边界在重复实现、平均或空间观察里保持相干。

下方是我对一行进行采样的结果：

![图 10. GFX100S 低码值平滑渐变的单帧与平均对比。测试区取 σ/Δ 最小的近黑端；单帧虽然随机跳动，500 帧平均后回到同一条理想渐变。](figures/bitdepth_original_08.png)

![图 11. 不同平均帧数下的 14-bit 与 16-bit 输出比较。在自然读噪充当 dither 时，低位深未稳定显出相对 16-bit 的阶梯边界。](figures/bitdepth_original_11.png)

![图 12. 单帧 16-bit 与 16→14-bit 的分排对比。](figures/single_frame_rndn_split_16_left_14_right.png)

上图使用单行剖面展示量化后的离散跳变，是为了把逐像素码值跳动直接暴露出来。实际视觉判断发生在二维空间中，那里人视觉处理流程会把一个最小可辨识圆内的多个像素共同积分，这个过程会削弱逐像素量化跳变的显著性；其效果取决于观察尺度、显示采样和局部空间频率。

另外，单行曲线在几个码值之间跳动，只说明数字样本是离散的；它本身并不等价于可见 banding。一个极端的例子比如下面的图像完全由纯白和纯黑组成：



![图 13. 二值半色调示例。局部只有黑白，但空间积分之后形成连续灰度感。](figures/bitdepth_original_09.png)

![图 14. 二值半色调的照片示例。视觉系统会在一定尺度上做空间整合。](figures/bitdepth_original_10.png)

二值半色调很直观地说明局部像素值可以极端离散，只要空间采样和视觉积分之后没有稳定的等级边界，人眼看到的仍然接近连续灰度。这当然不是 RAW 量化模型，但可说明几个码值上下跳并不等于一道可见台阶。

---

## 10. 结论

量化输出的特征函数是输入特征函数以 $2\pi/\Delta$ 为周期的周期化，于是一切与量化网格相干的统计量，都按输入特征函数在格点 $2\pi k/\Delta$ 处的值加权。对高斯读噪，这个权重是 $\exp[-2\pi^2k^2(\sigma/\Delta)^2]$，读噪因此成了一个指数有效的 dither。本文主要使用一阶和二阶两个条件矩：一阶定位确定性均值结构，二阶界定随机残差的能量调制。

在这个框架下，对 $16\!\to\!14$-bit：把确定性均值结构压到多个数量级范围内的低幅度判据以下，所需读噪不超过约 $2$ 到 $3.4\,\mathrm{DN16}$，且门槛对判据只是对数级敏感；GFX100S 近黑读噪 $\sigma=4.729\,\mathrm{DN16}$（$r\approx1.18$）已高于它，一阶确定性偏置约 $1.3\times10^{-12}\,\mathrm{DN16}$，二阶相位主谐波约 $9.5\times10^{-11}\,\mathrm{DN16}^2$，都远低于 1 DN16 量级，也低于本文图示和常规模拟能稳定分辨的量级。PDR 的代价也同样可忽略。$16\!\to\!14$-bit 约 $0.04$ stop，多数机型 $14\!\to\!12$-bit 在高增益区不到 $0.1$ stop，说明动态范围损失主要由噪声底决定，而非由输出位深单独决定。

所以在本文的模型与作用域里，14-bit RAW 与 16-bit RAW 对理想线性信号已经足够相似；16-bit 多出来的两位主要表现为读噪的随机摆动，可以作为浮点处理或平均过程里的数值容器，但不构成单张 RAW 在线性信号统计、PDR 损失或低阶量化相位结构上的实质优势。

本文核心模型是线性 RAW 域、最近邻量化、独立同分布的高斯读噪，以及普通的下采样，这和相机机内切换ad行为有细微偏差，后者后端噪声会略微高些。12-bit 的结论还依赖 ISO、读噪、模拟增益与 ADC 设计，若 $\sigma/\Delta$ 掉进危险区需把数据重新代回本文框架核算。

---

## 参考文献

量化与矩的经典理论

[1] A. Zaman, *The Density of the Fractional Part of a Normal Distribution*, Lahore University of Management Sciences, 2004.

[2] W. R. Bennett, "Spectra of Quantized Signals," *Bell System Technical Journal*, vol. 27, no. 3, pp. 446–472, 1948. doi:10.1002/j.1538-7305.1948.tb01340.x.

[3] W. F. Sheppard, "On the Calculation of the Most Probable Values of Frequency-Constants, for Data Arranged According to Equidistant Divisions of a Scale," *Proc. London Mathematical Society*, vol. 29, pp. 353–380, 1898.

[4] B. Widrow, "A Study of Rough Amplitude Quantization by Means of Nyquist Sampling Theory," *IRE Trans. Circuit Theory*, vol. 3, no. 4, pp. 266–276, 1956.

[5] B. Widrow, I. Kollár, and M.-C. Liu, "Statistical Theory of Quantization," *IEEE Trans. Instrumentation and Measurement*, vol. 45, no. 2, pp. 353–361, 1996. doi:10.1109/19.492748.

[6] B. Widrow and I. Kollár, *Quantization Noise: Roundoff Error in Digital Computation, Signal Processing, Control, and Communications*, Cambridge University Press, 2008.

[7] A. B. Sripad and D. L. Snyder, "A Necessary and Sufficient Condition for Quantization Errors to be Uniform and White," *IEEE Trans. Acoustics, Speech, and Signal Processing*, vol. 25, no. 5, pp. 442–448, 1977.

[8] R. M. Gray and D. L. Neuhoff, "Quantization," *IEEE Trans. Information Theory*, vol. 44, no. 6, pp. 2325–2383, 1998. doi:10.1109/18.720541.

[9] R. M. Gray and T. G. Stockham, "Dithered Quantizers," *IEEE Trans. Information Theory*, vol. 39, no. 3, pp. 805–812, 1993.

[10] A. G. Clavier, P. F. Panter, and D. D. Grieg, "Distortion in a Pulse Count Modulation System," *Trans. AIEE*, vol. 66, pp. 989–1005, 1947.

Dither 理论

[11] L. Schuchman, "Dither Signals and Their Effect on Quantization Noise," *IEEE Trans. Communication Technology*, vol. 12, no. 4, pp. 162–165, 1964. doi:10.1109/TCOM.1964.1088973.

[12] S. P. Lipshitz, R. A. Wannamaker, and J. Vanderkooy, "Quantization and Dither: A Theoretical Survey," *J. Audio Engineering Society*, vol. 40, no. 5, pp. 355–375, 1992.

[13] R. A. Wannamaker, S. P. Lipshitz, J. Vanderkooy, and J. N. Wright, "A Theory of Nonsubtractive Dither," *IEEE Trans. Signal Processing*, vol. 48, no. 2, pp. 499–516, 2000. doi:10.1109/78.823976.（机制细节亦见 R. A. Wannamaker, *The Theory of Dithered Quantization*, Ph.D. thesis, Univ. of Waterloo, 1997.）

[14] J. Vanderkooy and S. P. Lipshitz, "Resolution Below the Least Significant Bit in Digital Systems with Dither," *J. Audio Engineering Society*, vol. 32, no. 3, pp. 106–113, 1984.

传感器噪声与光子传递曲线

[15] J. R. Janesick, *Photon Transfer: DN → λ*, SPIE Press Monograph PM170, 2007. doi:10.1117/3.725073.

[16] EMVA, *Standard 1288: Standard for Characterization of Image Sensors and Cameras, Release 4.0 Linear*, European Machine Vision Association, 2021.

[17] M. Konnik and J. Welsh, "High-level numerical simulations of noise in CCD and CMOS photosensors: review and tutorial," arXiv:1412.4031, 2014.

[18] ISO 15739:2023, *Photography — Electronic still-picture imaging — Noise measurements*（含 Visual Noise / HVS 加权）.

位深感知、动态范围与编码

[19] P. G. J. Barten, *Contrast Sensitivity of the Human Eye and Its Effects on Image Quality*, SPIE Press, 1999.

[20] S. Miller, M. Nezamabadi, and S. Daly, "Perceptual Signal Coding for More Efficient Usage of Bit Codes," *SMPTE Motion Imaging Journal*, vol. 122, no. 4, pp. 52–59, 2013.

[21] DICOM PS3.14, *Grayscale Standard Display Function (GSDF)*, NEMA.

[22] Rec. ITU-R BT.2100 / Report ITU-R BT.2246（HDR 与 contouring 阈值，8/10/12-bit）.

[23] C. Poynton, *Digital Video and HD: Algorithms and Interfaces*, 2nd ed., Morgan Kaufmann, 2012.

本项目材料

[24] 姜尧耕, *14bit RAW Delta PDR Critical ISO Curves*, project page, <https://y-g-jiang.github.io/PLT.html>.

[25] 姜尧耕, *paper_IIARPTC*, GitHub repository, <https://github.com/y-g-jiang/paper_IIARPTC>.

[26] 姜尧耕, *在线 PTC 拟合与相机噪声分析图表材料*, 2026.

[27] Y. Jiang (姜尧耕), *Impact of ISP Affine Re-quantization on Photon Transfer Characteristics and Perceptual Quality*, Zenodo preprint, 2025. doi:10.5281/zenodo.18949317. 代码与更新：<https://github.com/y-g-jiang/paper_IIARPTC>。（给出高斯输入下 $\Delta^2/12$ 量化噪声近似成立的阈值——最坏相位主导偏离 $\delta(r)=(4r^2+1/\pi^2)e^{-2\pi^2r^2}$ 及其 Lambert W 近似阈值——与非整数仿射再量化的相位振荡/相移分析。）
