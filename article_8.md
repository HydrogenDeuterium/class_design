# 八、筛板塔设计计算

## 8.1 精馏段、提馏段设计计算

### 8.1.1 平均压力

取每层塔板的压降为$0.7\mathrm{kPa}$。由2.1知，塔顶表压为 $4\mathrm{kPa}$，因此绝对压力为$101.3+4=105.3\mathrm{kPa}$。

对于加料板，压力为$105.3+0.7\times8=110.9 \mathrm{kPa}$,因此精馏段平均压力为$(105.3+110.9)/2=108.1 \mathrm{kPa}$。

对于塔底，压力为$105.3+0.7\times20=119.3 \mathrm{kPa}$,因此提馏段平均压力为$(110.9+119.3)/2=115.1 \mathrm{kPa}$。

### 8.1.2 平均温度

由7.2，进料板温度为$88.23℃$；由7.3，塔顶温度为$80.29℃$;由7.3，塔底温度为$131.41℃$。

因此，$t_{m,精}=\frac{80.29+88.23}{2}=84.26℃$；$t_{m,提}=\frac{88.23+131.41}{2}=109.82℃$。

### 8.1.3 平均分子量

塔顶$y_1=x_D=0.993$,$x_1=y1*101.3/\mathrm{p}_a^0=0.993*101.3/101.9=0.968$

进料板$x_F=0.728$,$y_F=\mathrm{p}_A^0*x_A/101.3=0.935$

塔底$x_D=0.00288$,$y_D=\mathrm{p}_A^0*x_D/101.3=0.0115$

使用如下的代码计算平均分子量：

```Mathematica
MeanMmass[x_]:=(x*78.11+(1-x)*112.61)
Print@MeanMmass@{(0.968+0.728)/2,(0.993+0.935)/2}
Print@MeanMmass@{(0.728+0.00288)/2,(0.935+0.0115)/2}
```

精馏段液相平均分子量$\mathrm{M_{V,m}}=83.35\mathrm{kg/kmol}$，气相平均分子量$\mathrm{M_{V,m}}=79.35\mathrm{kg/kmol}$

提馏段液相平均分子量$\mathrm{M_{V,m}}=100.00\mathrm{kg/kmol}$，气相平均分子量$\mathrm{M_{V,m}}=96.28\mathrm{kg/kmol}$

### 8.1.4 液相平均密度

苯和氯苯理想混合，可视为理想流体。

塔顶：$t=80.29℃$,$x=0.987$

进料板：$t=88.23℃$，$x=0.728$

塔底：$t=131.41℃$,$x=0.00288$

使用如下的代码进行计算：

```Mathematica
rhoA[t_]:=912.13-1.186t
rhoB[t_]:=1124.4-1.0657t
rhoMix[t_,x_]:=1/(x/rhoA[t]+(1-x)/rhoB[t])
Print[rhoMix[84.26,(0.993+0.728)/2],rhoMix[109.82,(0.993+0.728)/2]]
```

<!--可以得到，精馏段平均密度$\rho_{L,m,精}=837.3\mathrm{kg/m^3}$，提馏段平均密度$\rho_{L,m,提}=807.1\mathrm{kg/m^3}$
这里我密度计算与他有小区别，会导致结果微微偏小。
-->

可以得到，精馏段平均密度$\rho_{L,m,精}=845.95\mathrm{kg/m^3}$，提馏段平均密度$\rho_{L,m,提}=815.6 \mathrm{kg/m^3}$

### 8.1.5 气相平均密度

压力较低，可以使用理想气体方程计算。

- 平均压力：$108.1\mathrm{kPa}$，$115.1\mathrm{kPa}$
- 平均分子量：$79.35\mathrm{kg/kmol}$，$96.28\mathrm{kg/kmol}$
- 平均温度：$84.26℃$，$109.82℃$

$$
\rho_V=\frac{\mathrm{p_mM_{V,m}}}{\mathrm{RT_m}}
   \rho_{V,1}=\frac{108.1\times79.35}{8.314\times(273.15+84.26)}\\
   =2.888\mathrm{kg/m^3}(精馏段)\\
   \rho_{V,2}=\frac{115.1\times96.28}{8.314\times(273.15+109.82)}\\
   =3.48\mathrm{kg/m^3}(提馏段)
$$

### 8.1.6 平均表面张力

使用如下的代码计算：

```Mathematica
{A80,B80}={21.23,25.93}
{A88,B88}={20.27,25.05}
{A131,B131}={15.16,20.34}
sigmaMix[A_,B_,x_]:=A*B/(A*(1-x)+B*x)
(sigmaMix[A80,B80,0.993]+sigmaMix[A88,B88,0.728])/2
(sigmaMix[A88,B88,0.728]+sigmaMix[A131,B131,0.00288])/2
```

可得精馏段平均表面张力$\sigma_{m,1}=21.32 \mathrm{mN/m}$，提馏段平均表面张力为$\sigma_{m},2=20.85\mathrm{mN/m}$。

### 8.1.7 平均粘度

使用如下代码计算

```Mathematica
{A80,B80}={0.307,0.427}
{A88,B88}={0.284,0.400}
{A131,B131}={0.197,0.291}
muMix[A_,B_,x_]:=A*x+B*(1-x)
(muMix[A80,B80,0.997]+muMix[A88,B88,0.728])/2
(muMix[A88,B88,0.728]+muMix[A131,A131,0.00288])/2
```

<!--可得精馏段平均粘度为$0.311 \mathrm{mPa\cdot s}$，提馏段平均粘度为$0.256 \mathrm{mPa\cdot s}$
晕 这小子自己算错了，舍入问题差了0.1-->

可得精馏段平均粘度为$0.312 \mathrm{mPa\cdot s}$，提馏段平均粘度为$0.256 \mathrm{mPa\cdot s}$

## 8.2 气液负荷

- 精馏段气相摩尔流率：$V=(R+1)D=1.492\times45.61=68.05\mathrm{kmol/h}$
- 精馏段气相体积流率：$V_R=\frac{VM_{V,m}}{\rho_{V,m}}=\frac{68.05\times79.35}{2.888}=1870\mathrm{m^3/h}=0.519\mathrm{m^3/s}$
- 精馏段液相摩尔流率：$L=RD=0.492\times45.61=22.44\mathrm{kmol/h}$
- 精馏段液相体积流率：$L_R=\frac{VM_{L,m}}{\rho_{L,m}}=\frac{22.44\times83.35}{845.95}=2.196\mathrm{m^3/h}=0.0006100\mathrm{m^3/s}$
- 冷凝器负荷：$2.0974\times10^6\mathrm{kJ/h}$

- 提馏段气相摩尔流率：$V=68.05+(1-q)F=68.05\mathrm{kmol/h}$
- 提馏段气相体积流率：$V_R=1875\mathrm{m^3/h}=0.52\mathrm{m^3/s}$
- 提馏段液相摩尔流率：$L=22.44+QF=84.71\mathrm{kmol/h}$
- 提馏段液相体积流率：$L_R=\frac{VM_{L,m}}{\rho_{L,m}}=\frac{84.71\times100}{924.0}=9.17\mathrm{m^3/h}=0.00255\mathrm{m^3/s}$
- 再沸器负荷：$2.59\times10^6\mathrm{kJ/h}$

## 8.3 塔径

1. 取塔板间距$H_T=500\mathrm{mm}$，板上液层高度$H_L=47\mathrm{mm}$，则板上空间为$453\mathrm{mm}$。

2. 以Smith法求空塔、泛点气速：

对精馏段：

$$
\left(\frac{\mathrm{L_s}}{\mathrm{V_s}}\right)\left(\frac{\rho_{\mathrm{L}}}{\rho_{\mathrm{V}}}\right)^{0.5}
=\left(\frac{0.00061}{0.519}\right)\left(\frac{845.95}{2.888}\right)^{0.5}=0.0201
$$

查表[^chart]可得，$C_{f20}=0.092$。负荷因子表面张力校正：$C_f=C_{f20}\left(\frac{\sigma}{20}\right)^{0.5}=0.0950$

[^chart]:管国锋等.化工原理（第四版）[M].化学工业出版社，2015:360，图8-25

故泛点气速：$u_{\max}=C\left(\frac{\rho_L-\rho_V}{\rho_V}\right)^{0.5}=0.0258\times17.09=1.623\mathrm{m/2}$

对提馏段：

$$
\left(\frac{\mathrm{L}_{s}}{\mathrm{V}_{\mathrm{s}}}\right)\left(\frac{\rho_{\mathrm{L}}}{\rho_{\mathrm{V}}}\right)^{0.5}
=\left(\frac{0.00255}{0.52}\right)\left(\frac{924.0}{3.48}\right)^{0.5}=0.0797
$$

查表[^chart]可得，$C_{f20}=0.025$。负荷因子表面张力校正：$C_f=C_{f20}\left(\frac{\sigma}{20}\right)^{0.5}=0.0258$

故泛点气速：$u_{\max}=C\left(\frac{\rho_L-\rho_V}{\rho_V}\right)^{0.5}=0.0258\times16.264=0.420\mathrm{m/s}$

3. 综合考虑，泛点气速取为更低的泛点气速$0.420\mathrm{m/s}$的75%，即$0.315\mathrm{m/s}$
4. 精馏段塔径$D=\sqrt{\frac{4V_s}{\pi u}}=1.45\mathrm{m}$

取整为1600mm，此时操作气速为$4\times0.52/(\pi*1.6^2)=0.259\mathrm{m/s}$，塔横截面积为$8.042\mathrm{m^2}$。

## 8.4 塔板工艺结构尺寸的设计与计算

1. 溢流装置

采用单溢流型的平顶弓形溢流堰、弓形降液管、平行受液盘，且不设进口内堰。

   1. 出口堰长$l_w$

   取$l_w=0.7D=1.12\mathrm{m}$，此时,溢流强度$E=L_R/l_w=9.17/1.12=8.19\mathrm{m^2/h}$，其中$L_R$为液相体积流量(提馏段)。

   2. 出口堰高$h_W$

   由$l_w/D=0.7$，$L_R/l_w^{2.5}=2.196/1.12^{2.5}=1.65$，查表[^chart2]可得液流收缩系数$E=1.01$.

   [^chart2]:管国锋等.化工原理（第四版）[M].化学工业出版社，2015.356：图8-19

   于是堰上溢流高度$h_{ow}=8.19\times10^-3\times1.01\left(\frac{L_R}{l_w}\right)^{\frac{2}{3}}=8.27\times 10^{-3}\times8.19^{\frac{2}{3}}=0.0336\mathrm{m}\geq 0.006\mathrm{m}$

   因此出口堰高$h_W=H_L-h_{ow}=0.06-0.0336=0.0264\mathrm{m}=26.4\mathrm{mm}$

   3. 降液管宽度与面积
   由$l_w/D=0.7$，查表[^chart3]可得$W_d/D=0.14$，$A_f/A_T=0.087$，其中$A_T=2.01\mathrm{m^2}$为塔内空间的横截面积.

   [^chart3]:管国锋等.化工原理（第四版）[M].化学工业出版社，2015.354：图8-17

   即：宽度$W_d=0.14D=0.125\mathrm{m}$，面积$A_f=2.01\times0.087=0.175\mathrm{m^2}$

   此时，液体在降液管内停留的时间$\tau=A_fH_T/L_R=0.175\times0.5/0.00255=34.3 \mathrm{s} \geq 5\mathrm{s}$，符合条件。

   4. 降液管底隙高度

   液体通过降液管底隙的流速一般为0.07~0.25 m/s，取液体通过降液管底隙的流速$u_{0,L}=0.08 \mathrm{m\cdot s^{-1}}$，则
   $$
   h_0=\frac{L_{r,s}}{L_Wu_{0,L}}=\frac{0.00255}{1.12\times0.08}=0.02846\mathrm{m}\\
   =28.5\mathrm{mm}\gt 25\mathrm{mm}
   $$

2. 塔板布置

   1. 塔板分布

   由塔径为$1600\mathrm{mm}$，因此塔板宜作 4 块安装。

   2. 边缘区、安定区宽度

   本设计取边缘区宽度$W_C=60\mathrm{mm}$，安定区宽度$W_S=100\mathrm{mm}$。

   3. 开孔区面积

   $$
   r=D/2-W_C=0.8-0.060=0.740\mathrm{m}\\
   x=D/2-W_d-W_S=0.8-0.224-0.100=0.675\mathrm{mm}
   $$

   因此开孔区面积$A_a=2\left[x \sqrt{r^2-x^2}+\frac\pi{180} \arcsin\left(\frac\pi r\right) \right]\\$
   $=2\left[0.476+\frac\pi{180}0.74^2\arcsin\left(\frac{0.476}{} \right) \right]$
   $=1.304\mathrm{m^2}$

   4. 开孔数与开孔率

   取筛孔半径为$d_0=3.5\mathrm{mm}$，正三角形排列；筛板使用碳钢，厚度$\delta=3\mathrm{mm}$，孔的距径比$t/d_0=4$，即$t=3.5*4=14\mathrm{mm}$。

   每层塔板的开孔数$n=\left(\frac{1158\times10^3}{t^2}\right)A_a=\left(\frac{1158\times10^3}{14^2}\right)1.304=7704$(个)

   塔板开孔率$\varphi=\frac{0.907}{(t/d_0)^2}=\frac{0.907}{4^2}=0.0567\in$[5%,15%]，满足条件。

   每层塔板的开孔面积$A_0=\varphi A_s=0.0567\times 1.304=0.0739\mathrm{m^2}$

   气体通过筛孔的气速$u_{0,V}=V_{R,s}/A_0=0.52/0.0567=9.17\mathrm{m\cdot s^{-1}}$

3. 塔高
    - 精馏段：$Z_1=(N_{p1}-1)H_T=7\times0.5=3.5\mathrm{m}$
    - 提馏段：$Z_2=(N_{p2}-1)H_T=11\times0.5=5.5\mathrm{m}$
    - 总高：$Z=Z_1+Z_2=9\mathrm{m}$

## 8.5 塔板上的流体力学验算

### 8.5.1 气体通过筛板压降的验算

1. 气体通过干板的压降

由$\delta/d_0=0.86$，开孔率$\varphi=0.0567$，查表[^chart851]可知，孔流系数$C_0=0.78$。

[^chart851]:管国锋等.化工原理（第四版）[M].化学工业出版社，2015.357：图8-21

$$
h_d=0.051\left(\frac{u_{0,V}}{C_0}\right)^2\frac{\rho_v}{\rho_L}\\
h_{d,1}=0.051\left(\frac{9.17}{0.78}\right)^2\frac{3.01}{843.9}
       =0.0251\mathrm{m}(精馏段)\\
h_{d,2}=0.051\left(\frac{9.17}{0.78}\right)^2\frac{3.48}{943.0}
       =0.0260\mathrm{m}(提馏段)
$$

2. 气体通过板上液层的压降

$h_L=\beta(h_w+h_{ow})=\beta H_L$，其中$H_L$为板上的液层高度。

对于单流型塔，有效截面的空气速$u_a=\frac{V_{R,s}}{A_T-2A_f}=\frac{0.52}{2.01-2\times0.175}=0.313\mathrm{m\cdot s^{-1}}$

因此，动能因子$
F_a=u_a\sqrt{\rho_V}$

$$
F_{a,1}=0.313\times \sqrt{3.01}=0.543\\
F_{1,2}=0.313\times \sqrt{3.48}=0.584
$$

查表[^chart8512]求液层充气系数$\beta$，如下：$\beta_1=0.68,\beta_2=0.66$。

[^chart8512]:管国锋等.化工原理（第四版）[M].化学工业出版社，2015.358：图8-23

因此，通过板上液层的压降
$$
h_f=\beta H_L\\
h_{f,1}=0.0326\mathrm{m}\\
h_{f,2}=0.0351\mathrm{m}
$$

3. 单板压降

$$
h_f=h_d+h_L\\
h_{f,1}=0.0251+0.0326=0.0577\mathrm{m}\\
h_{f,2}=0.0260+0.0351=0.0611\mathrm{m}
$$

因此压降$\Delta P_f=\rho_L g h_f$，也即：

$$
\Delta P_{f,1}=843.9\cdot g \cdot 0.0420=477.7\mathrm{Pa}\\
\Delta P_{f,2}=924.0\cdot g \cdot 0.0448=553.8\mathrm{Pa}
$$

$\Delta P_f<0.7\mathrm{kPa}$，符合设计条件。

### 8.5.2 雾沫夹带量的验算

塔内气速$u_n=\frac{V_{R,s}}{A_T-A_f}=\frac{0.52}{2.01-0.175}=0.28\mathrm{m\cdot s^{-1}}$

取板上泡沫层厚度为液层的2.5倍，即$H_f=2.5H_L=60\times2.5=150\mathrm{mm}$

则雾沫夹带量$e_V=\frac{5.7\times10^{-6}}{\sigma_{m}}\cdot \left(
   \frac{u_n}{H_T-H_f}
\right)$，故有：

$$
e_{v,1}=\frac{5.7\times10^{-6}}{21.26}\cdot
\left(
   \frac{0.28}{0.5-0.15}
\right)^{3.2}\\
=1.28\times10^{-7}\mathrm{kg(Liquid)\cdot kg^{-1}(Gas)}\\

e_{v,2}=\frac{5.7\times10^{-6}}{22.85}\cdot
\left(
   \frac{0.28}{0.5-0.15}
\right)^{3.2}\\
=1.22\times10^{-7}{kg(Liquid)\cdot kg^{-1}(Gas)}
$$

计算结果不大于$0.1\mathrm{kg(Liquid)\cdot kg^{-1}(Gas)}$，满足要求，不发生过量液沬夹带。

### 8.5.3 漏液限计算

按如下经验公式计算漏液点气速：
$$
u_{om}=4.4C_0\sqrt{(0.0056+0.13H_L-h_\sigma)\rho_L/\rho_V}\\
$$

其中，$h_\sigma$为克服筛孔处界面张力产生的压降（以清液柱高度计算），计算公式为：$h_\sigma=\frac{4\times10^{-3}\sigma_m}{\rho_L d_0}$：

$$
h_{\sigma,1}=\frac{4\times10^{-3}\times21.26}{843.9\times9.81\times 0.004}=0.00293\mathrm{m}\\

h_{\sigma,2}=\frac{4\times10^{-3}\times22.85}{924.0\times9.81\times 0.004}=0.00288\mathrm{m}\\
$$

带入可得：

$$
u_{om,1}=4.4\times0.78\sqrt{(0.0056+0.13\times0.06-0.00293)843.9/3.01}\\=5.88\mathrm{m\cdot s^{-1}}\\

u_{om,1}=4.4\times0.78\sqrt{(0.0056+0.13\times0.06-0.00288)924.0/3.48}\\=5.74\mathrm{m\cdot s^{-1}}\\
$$

此时$u_0=9.17\mathrm{m\cdot s^{-1}}=1.5\times6.11\mathrm{m\cdot s^{-1}}$，符合条件。

### 8.5.4 液泛的计算

为防止液泛，应该使得降液管清液高度$H_d\le(H_T+h_W)$，而$H_d=h_f+H_L+\Sigma h_f$。

降液管阻力：$\Sigma H_f=0.153\left(\frac{L_{R,s}}{L_Wh_0}\right)^2$，故：

$$
\Sigma H_{f,1}=0.153\left(\frac{0.00061}{1.12\times0.02846}\right)
=3.62\times10^{-11}\mathrm{m}\\

\Sigma H_{f,1}=0.153\left(\frac{0.00255}{1.12\times0.02846}\right)
=6.33\times10^{-10}\mathrm{m}
$$

因此，管内清液层高度$H_d=h_f+H_L+\Sigma h_f$，有：
$$
H_{d,1}=0.0577+0.06+0.00293=0.121\mathrm{m}\\
H_{d,2}=0.0611+0.06+0.0120=0.133\mathrm{m}
$$

取泡沫相对密度$\Phi=0.5$，故允许的最大清液层高度$\Phi(H_T+h_W)=0.2632\mathrm{m}$

因此$H_d\le \Phi(H_T+h_W)$成立，故，不会发生液泛。

经过以上的流体力学验算，可以认为，精馏塔塔径与塔板工艺尺寸合适。

## 8.6 精馏、提馏段塔板负荷性能图

负荷性能图的绘制过程如下：

### 8.6.1 雾沫夹带线

液沬夹带量按照如下的经验公式计算：

$$
e_V=\frac{5.7\times10^{-6}}{\sigma_{m}}\cdot \left(
   \frac{u_n}{H_T-H_f}
\right)^{3.2}
$$

其中，$u_n=\frac{V_{R,s}}{A_T-A_F}=0.545V_{R,s}$。

取塔上泡沫层厚度$H_f$为板上液层厚度的2.5倍，则：

<!--下面算出来的东西只和L_W有关，我是这么理解的-->

$$
\mathrm{H}_{\mathrm{f}} =2.5\mathrm{H}_{\mathrm{L}}=2.5\left(
   \mathrm{h}_{\mathrm{w}}+\mathrm{h}_{\mathrm{ow}}\right)
   =2.5 \times\left[0.0264+0.00284\left(\frac{3600 \mathrm{L}_{\mathrm{R,s}}}{1.12}\right)^{\frac23}\right]
   \\
   =0.121+1.546L^{\frac23}\mathrm{m}
$$

取$e=0.1$，化简有:

$$
V_1=8.299-33.877 L_1^{\frac23}\\
\\\\
V_2=8.114-33.112 L_2^{\frac23}
$$

### 8.6.2 液泛线（气相负荷上限线）

液泛条件：$\Phi(H_T+h_w)=h_f+h_w+h_{ow}+\Sigma H_f$，其中泡沫相对密度$\Phi=0.5$

板压降$h_f=h_d+h_L$，其中$h_d=0.051\left(\frac{u_{0,V}}{C_0}\right)^2\frac{\rho_v}{\rho_L}$：
$$
h_{d1}=0.051\left(\frac{V_1}{0.78\times0.0739}\right)^2
   \frac{3.01}{843.9}=0.0547V_1^2\\
h_{d2}=0.051\left(\frac{V_2}{0.78\times0.0739}\right)^2
   \frac{3.48}{943.0}=0.0566V_2^2
$$

液泛条件可改写如下：$$\Phi(H_T+h_w)=h_d+(\beta+1)(h_w+h_ow)+\Sigma H_f$$

而$h_L=\beta(h_w+h_{ow})=\beta H_L$,有:

$h_w=0.0264\mathrm{m}$，$\beta可取0.58$。

$h_{ow}=0.00284 E \left(\frac{3600L}{L_W}\right)^{\frac{2}{3}}$，其中液流收缩系数$E$不妨取$1$，有：$h_{ow}=0.618L^{\frac{2}{3}}=0.6185L^{\frac{2}{3}}$

$\Sigma H_f=0.153\left(\frac{L}{L_Wh_0}\right)^2=0.153\left(\frac{L}{1.12\times 0.02846}\right)^2=150.6L^2$

综上，有：
$$
0.5\times(0.5+0.0264)=0.0547V_1^2+1.58(0.0264+0.6185L_1^\frac{2}{3})+150.6L_1^2\\
0.5\times(0.5+0.0264)=0.0566V_2^2+1.58(0.0264+0.6185L_2^\frac{2}{3})+150.6L_2^2
$$

化简，得：
$$
0.041712+0.97723 L_1^\frac{2}{3}+150.6L_1^2+0.0547V_1=0.2632\\
0.041712+0.97723 L_2^\frac{2}{3}+150.6L_2^2+0.0566V_2=0.2632\\
V_1 = 4.049 - 17.87 L_1^\frac{2}{3} - 2753 L_1^2\\
V_2 = 3.913 - 17.27 L_2^\frac{2}{3} - 2660 L_2^2\\
$$

### 8.6.3 漏液线

漏液点气速为：$=4.4C_0\sqrt{(0.0056+0.13h_L)\frac{\rho_L}{\rho_V}}$(筛孔界面张力忽略不计)

其中，$h_L=h_W+h_{ow}=0.0264+0.6185L^{\frac{2}{3}}$

化简可得：
$$
V_1^2=0.102+ 0.0145L_1^{\frac{2}{3}}\\
V_2^2=0.096+ 0.0137L_2^{\frac{2}{3}}
$$

### 8.6.4 液相负荷下限线

取液流收缩系数$E=1.0$，堰上液流高度$h_{0w}=0.006\mathrm{m}$，即有：

$0.61855L^{\frac{2}{3}}=0.06$，解之可得：$L=9.55\times10^{-4}\mathrm{m^{3}\cdot s^{-1}}$

### 8.6.5 液相负荷上限线

取$\tau_{\min}=5\mathrm{s}$，则
$$
L_{s,\max}=\frac{H_TAf}{\tau}=\frac{0.5\times0.175}{5}=0.0175\mathrm{m^{3}\cdot s^{-1}}
$$

### 8.6.6 操作线与操作弹性

由8.2可知，精馏段操作气液比$V_R/L_R=85.2$，提馏段$V_R/L_R=204$。

绘图代码如下：

```WL
wm := 8.299 - 33.877 L^(2/3)
yf := 4.049 - 17.87 L^(2/3) - 2753 L^2
ly := 0.102 + 0.0145 L^(2/3)
op := 85.2 L
yx := Line[{{9.55*10^(-4), 0}, {9.55*10^(-4), 8}}]
ys := Line[{{0.0175, 0}, {0.0175, 8}}]
Legends =Placed[{"雾沫夹带线", "液泛线", "漏液线", "操作线", "液相下限线", "液相上限线"}, Above];
styleColor = {Blue, Red, Green, Purple, White, White};

Plot[{wm, yf, ly, op, 0, 8.5}, {L, 0, 0.02},
 PlotLabels -> Legends,
 PlotStyle -> styleColor,
 PlotRange -> Full,
 Filling -> {3 -> {2}},
 AspectRatio -> 1,
 Epilog -> {{{Cyan, yx}, {Orange, ys}}}]
```

# 参考文献
