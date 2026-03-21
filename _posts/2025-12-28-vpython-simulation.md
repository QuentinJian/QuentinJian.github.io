---
layout: post
title:  "Rk4-In-Vpython"
date:   2025-12-28
categories: numerical-method python
---

在高中物理課時我曾學過基本的 vpython 做物理模擬，而之後便再也沒碰過。一直到普通物理課時老師在講解作業提到了用 runge kutta4 integrator 畫 phase diagram可以加分，因此寫了這個 vpython 程式。

## 4th Order Runge Kutta

Runge Kutta Method 是一個解非線性常微分方程近似的一系列方法，這些方法中以 4th order runge kutta (RK4)最常被使用，因為利用RK4 我們可以用電腦算出許多難以計算的微分方程的近似解，也因此可以被用來做物理模擬。解微分方程就可以模擬物理是因為描述物理行為的數學式通常以微分式表示，像是 ` F=ma ` 可以寫成 `F = m(dv/dt) => dv = dt * F/m 而 v = dx/dt` 只要知道當下的受力與速度，就可以知道下一個瞬間的速度，知道當下的速度，就可以知道下一個瞬間的位置
利用這個想法，我們可以寫出一個十分直觀的程式:

```python
x = 0
t = 0
dt = 0.001
v = 1 #(m/s)
F = 10 #(N)
m = 1 #(kg)
while (t < 10):
  x = x + v * dt
  v = v + (F/m) * dt
```

來簡單估算 10 秒後物體的位置、速度等參數，而上面的這個簡單方法稱為歐拉法(Euler method)，當我們給的dt值越小，計算的結果也就越精確。

### 為何歐拉法不夠好

利用歐拉法，大部分基本的物理模擬都能十分貼近現實的表現，過程產生的誤差相對十分渺小，但當我們遇到要模擬複雜一點的系統時就會出問題，例如模擬考慮摩擦力的彈簧振盪時會發現彈簧的能量竟然不減反增。

```python
while (t &lt; 100):
  rate(1000)
  square.a=(-K*(spring.length-spring.L0)*spring.axis.norm() - (0.1)*square.v) / m
  square.pos += square.v*dt
  square.v += square.a*dt
  
  
  spring.axis= square.pos-spring.pos
  
  plot_graph(square)
  update_arrow(square, 5)
  t += dt
```

![test result]({{ site.baseurl }}/res/image.png)

### Runge Kutta Methods

Runge Kutta 法除了這篇文章主要介紹的 RK4 外，其實Runge Kutta Methods 指的是一系列求微分方程數值解的方法，而他基本上也是前面 Euler Method 的一種延伸，相較於前面只針對 `t` 與 `t + dt` 便求得下一瞬間的 `x` 值，Runge Kutta 法是利用 `t` 與 `t + dt` 之間的不同時間計算出結果後進行加權得到的，因此直觀上精準度也隨著提高，而背後大概的原理來自於利用已知的數值對函式進行模擬泰勒展開後取其中前幾項進行逼近，基本的Euler Method 屬於取第一項，有點像只看起點的斜率，而RK4則是進行四次取樣因此可以更接近的重現原本函數的曲線。

前面的程式利用RK4重寫後便成為這樣

```python
def evaluate(initial: State, t: float, dt: float, d: Derivative) -> Derivative:
    state = State()
    state.x = initial.x + d.dx * dt
    state.v = initial.v + d.dv * dt
    
    
    output = Derivative()
    output.dx = state.v
    
    F = -K*(state.x-spring.L0)*spring.axis.norm()
    
    output.dv = F/m
    
    return output

def RK4(state: State, t: float, dt: float) -> State:
    k1 = evaluate(state, t, 0, Derivative())
    k2 = evaluate(state, t + dt/2, dt/2, k1)
    k3 = evaluate(state, t + dt/2, dt/2, k2)
    k4 = evaluate(state, t + dt, dt, k3)
    
    state.x += (dt/6.0) * (k1.dx + 2*k2.dx + 2*k3.dx + k4.dx)
    state.v += (dt/6.0) * (k1.dv + 2*k2.dv + 2*k3.dv + k4.dv)
    
    return state

while t &lt; t_max:
    rate(1000)
    
    state = RK4(state, t, dt)

    square.pos.x = state.x
    
    spring.axis= square.pos-spring.pos

    t += dt
```

換成這個寫法後，模擬出來有比前面的結果好了一些。

![alt text]({{ site.baseurl }}/res/RK4Damper.png)
