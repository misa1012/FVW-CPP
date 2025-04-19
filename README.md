# FVW-CPP

## 代码实现说明
代码的主要实现在main.cpp

### 第一部分、设置参数
代码的主要实现，首先是对参数进行设置。
参数分为
1. 模拟相关的参数simParams：包括dt和totalTime
2. Turbine相关的参数turbineParams


### 第二部分、初始化
初始化包括
1. 对叶片几何进行初始化：
包括了对shed节点和trailing节点的初始化：r, chord和twist
以及计算叶片前缘，lifting line和后缘的shed和trailing的坐标 （这里的坐标是叶片坐标系下的坐标）

2. 读取翼型dataset的数据

3. 对position进行初始化
计算前缘，lifting line和后缘的shed和trailing，hub和platform的全局坐标
该部分主要实现的是DCMRot算法

（本来有删除了的原4.） 对velocity进行初始化 -- 产生了一个疑问，这个真的需要吗？-- 感觉好像不需要啊。删掉了。

4. 初始化performance
初始化aoa, cl, cd的空的列表

-- 以上已进行验证 --
验证结论是：interpolate得到的数值和python中的不太一样，其它验证是一致的

### 第三部分、BEM - 用于找到初始的aoa和cl、cd
该文件中有两个函数：计算aoa和计算BEM
BEM的原理在于它要同时满足动量理论和叶素理论：
动量理论：诱导因子决定了速度场
叶素理论：cl和cd取决于aoa

所以计算步骤为，需要迭代计算：
1. 猜测一个induction factor诱导因子
2. 根据该值计算phi（入流角度）和angle of attack
3. 根据aoa的值，插值airfoil data找到cl, cd （叶素理论）
4. 优化：tip loss correction 
5. 根据cl, cd计算出推力系数，进行计算出induction factor（动量理论）
6. 根据该新计算的induction factor计算之前的（需要用weight factor更新）

BEM的作用：通过BEM计算出的CL，计算出初始环量
所以BEM更新的变量是perf


### 第四部分、更新wake （还在开发中）

## 变量解释