# Boltzmann entropy code

## 1. Calculating Q

$$
Q_k = \sum_l B(k,l) f_l f_{k-l}
$$



1d example

二维数组计算效率更高，或者，按i,j遍历A(i,j)会更快

`int`也比`u_int`快

``k->l f[k-l]`代表逆序遍历`f`, `f[]`代表顺序遍历

`vector<>, a[]` `vector<vector<double> >, a[][]` 

```cpp
for k
  for l
    Q[k] += B[k*N+l]*f[l]*f[k-l]
    
for l  0:N
  for k 0:N
    Q[k] += B[l*N+k]*f[l]*f[k-l+n+N]
    
 k-l < 0 +N, >= N -N
```

```sh
g++ -o main -O3 
```

```
k = -n,...,0,...,n, 0,...n,...2n, 0,...,6n+2. N = 2n+1, 3N,   ki = k+n
l = -n,...,0,...,n, li = l+n, 
m = k-l= (ki-n)-(li-n),    m1 =m+n =  ki-li+n
```



| n=20001/CPUTime(s),O3(O0)                | `k->l（f[k-l]）` | `l->k（f[k-l]）` |
| ---------------------------------------- | ---------------- | ---------------- |
| 二维数组存储                             | 0.5 (4.46)       | 9.72(30.42)      |
| 一维数组按列存储                         | 0.50(3.89)       | 5.26(11.42)      |
| 一维数组按对应行或列存储                 | 0.50(3.89)       | 0.39(3.89)       |
| 一维数组按列存储+删去`k-l`的判断         | 0.49(3.58)       | 5.18(11.20)      |
| 一维数组按对应行或列存储+删去`k-l`的判断 | 0.49(3.58)       | 0.30(3.59)       |
| 一维数组按对应行或列存储+3倍f            |                  | 0.30(3.53)       |

结论：的确，遍历`k,l`的顺序，严重决定运行的时间，差别在优化后甚至有数十倍

2d example

````
B[i1,i2,j1,j2] B[i,j]
i1 = i/N; i2 = i %N; i2=i - i1*N;
````



| 1d数组+行或列存储 | `i->j（f[n-j]）` | `j->i（f[i]）`+4层loop | `j->i（f[i]）`+2层 | $3^2$倍f |
| ----------------- | ---------------- | ---------------------- | ------------------ | -------- |
| n=121*121         | 0.29             | 0.26                   | 1                  | 0.18     |
| n=151*151         | 0.69(6.12)       | 0.63(6.14)             | 2.38(7.4)          | 0.41     |
| n=161*161         | 0.9              | 0.81                   | 3.06               | 0.54     |
| n=181*181         | 1.43             | 1.28                   | 4.91               | 0.85     |
|                   | 100%             | 90%                    | 340%               | 60%      |



3d example

| 1d数组+行或列存储 | `i->j（f[n-j]）` | `j->i（f[i]）`+2层 | $3^3$倍f   |
| ----------------- | ---------------- | ------------------ | ---------- |
| n=$25^3$          | 0.33             | 0.32               | 0.24       |
| n=$27^3$          | 0.53(4.82)       | 0.51(4.81)         | 0.38(4.56) |
| n=$29^3$          | 0.81             | 0.77               | 0.58       |
| n=$31^3$          | 1.29             | 1.15               | 0.86       |
| n=$33^3$          | 1.77             | 1.75               | 1.24       |
| n=$35^3$          | 2.50             | 2.43               | 1.72       |
|                   | 100%             | >96%               | 70%        |



### 1.1 配合fftw的顺序

$$
Q_k = \sum_l B(k,l) f_l f_{k-l}
$$

1d数组+2倍长度

```sh
指数i 0,1,...,n,  -n,...,-1,   0,   1,...,   n,   -n,...,-1
指标k 0,1,...,n, n+1,...,2n,2n+1,2n+2,...,3n+1, 3n+2,...,4n+1
关系 k = i + N; N =2n+1
指数m = （i - i'） % N = （k - k'）%N
指标m = m + N =  k - k' + N
```

```cpp
for l n+1:3n+1
  for k n+1:3n+1
    Q[k] += B[(l-n-1)*N+(k-n-1)] * f[l] * f[k-l+N]
```
最后需要补充完整Q

## 2. Evaluation of B

[**G. Dimarco and L. Pareschi 2014**] Section 5.1
$$
\hat{Q}_{k}=\sum_{l, m=-N \atop l+m=k}^{N} \hat{f}_{l} \hat{f}_{m} \hat{\beta}(l, m), \quad k=-N, \ldots, N
$$
where the Boltzmann kernel modes $\hat{\beta}(l, m)=\hat{B}(l, m)-\hat{B}(m, m)$ are now given by 
$$
\hat{B}(l, m)=\int_{\mathcal{B}_{0}(2 \lambda \pi)} \int_{\mathbb{S}^{2}}|q| \sigma(|q|, \cos \theta) e^{-i\left(l \cdot q^{+}+m \cdot q^{-}\right)} d \omega d q
$$
In the VHS case, $|q| \sigma(|q|, \cos \theta)=C_{\alpha}|q|^{\alpha}$ where we can chose $C_{\alpha}=\left((4 \pi)^{2}(2 \lambda \pi)^{3+\alpha}\right)^{-1}$. Furthermore, the integral can be reduces to a one-dimensional integral (Pareschi and Russo 2000b)
$$
\hat{B}(l, m)=\int_{0}^{1} r^{2+\alpha} \operatorname{Sinc}(\xi r) \operatorname{Sinc}(\eta r) d r=F_{\alpha}(\xi, \eta)
$$
where $\xi=|l+m| \lambda \pi, \eta=|l-m| \lambda \pi$, $\lambda = 2/(3+\sqrt{2})$. We report the expressions for the case of Maxwell molecules $\alpha = 0$ and hard spheres $\alpha = 1$
$$
\begin{aligned} F_{0}(\xi, \eta) &=\frac{p \sin (q)-q \sin (p)}{2 \xi \eta p q} \\ F_{1}(\xi, \eta) &=\frac{q \sin (q)+\cos (q)}{2 \xi \eta q^{2}}-\frac{p \sin (p)+\cos (p)}{2 \xi \eta p^{2}}-\frac{2}{p^{2} q^{2}} \end{aligned}
$$
where $p=(\xi+\eta), q=(\xi-\eta)$.





