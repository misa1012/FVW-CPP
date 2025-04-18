# FVW-CPP

## 代码实现说明
代码的主要实现在main.cpp

第一部分、设置参数
代码的主要实现，首先是对参数进行设置。
参数分为
1. 模拟相关的参数simParams：包括dt和totalTime
2. Turbine相关的参数turbineParams


第二部分、初始化
初始化包括
1. 对叶片几何进行初始化：
包括了对shed节点和trailing节点的初始化：r, chord和twist
以及计算叶片前缘，lifting line和后缘的shed和trailing的坐标 （这里的坐标是叶片坐标系下的坐标）

2. 读取翼型dataset的数据

3. 对position进行初始化
计算前缘，lifting line和后缘的shed和trailing，hub和platform的全局坐标
该部分主要实现的是DCMRot算法

-- 以上已进行验证 --
验证结论是：interpolate得到的数值和python中的不太一样，其它验证是一致的

（本来有删除了）. 对velocity进行初始化 -- 产生了一个疑问，这个真的需要吗？-- 感觉好像不需要啊。删掉了。






## 变量解释