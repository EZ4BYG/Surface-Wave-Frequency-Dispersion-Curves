# 2D Surface Wave Modeling

**Contents**: 2D Forward Modeling of Rayleigh and Love Wave MATLAB program (Basic program). Surface waves play an extremely important role in shallow surface geophysical exploration. Surface waves account for more than 80% of the total energy of all elastic waves on the shallow surface (50m above)! It's really important to study them in depth.

**Tips**:
- The absorption boundaries in the program are PML and free surfaces
- Many parameters in the program can be shown in diagrams, just need a simple modification. Plus, the velocity parameter V is written into the Excel in the snapshot result of 130ms wave field
- Snapshots of the wave field at different times are also in the list

# Forward Simulation of Surface Wave Dispersion Curve

**Contents**: Forward modeling of Love wave and Rayleigh wave dispersion curves based on Haskell and Knopoff algorithms. 

**Intro**: Surface wave is a kind of seismic wave propagating along the free interface. It has the characteristics of low speed, low frequency, high amplitude and high frequency dispersion. Dispersion means that surface waves of different frequencies travel at different speeds! Moreover, the calculation of dispersion curve of surface wave in multilayered media is very important for surface wave inversion. At present, surface wave exploration mainly USES the dispersion characteristics of surface wave.

**Tips**: 
- The programs are in the **Frequency_Dispersion** folder, detailed information can be found in its readme.md file
- The Matlab version completes the basic functions, and the Python version are aim to parallel computing!
- Reference: *High-Frequency Surface-Wave Methods*, Jianghai Xia, China University of Geosciences Press

----

# 二维面波正演模拟

内容：本文中提供了Rayleigh和Love波的二维正演程序。面波在浅地表地球物理勘探中发挥极其重要的作用！浅地表(50m以上)所有弹性波中面波占总能量的80%以上！

**注意**：
- 程序中的吸收边界是PML和自由表面；
- 程序中很多参数分量都可以做图展示，只需做简单修改即可；另外速度参数V在130ms的波场快照结果写入excel表格中。
- 不同时间的波场快照图也在列表中

# 面波频散正演模拟

内容：基于Haskell和Knopoff算法的Love波和Rayleigh波频散曲线正演模拟 

介绍：面波是一种沿自由界面传播的地震波，具有低速、低频、高振幅、高频散的特点。频散是指不同频率的面波的传播速度不同！并且多层层状介质中面波频散曲线的计算对于面波的反演至关重要。目前面波勘探主要利用的是面波频散特性。

**注意**: 
- 所有相关程序在**Frequency_Dispersion**文件夹中，详细信息查看其中readme.md文件
- Matlab版本的程序完成所有基本功能，Python版本程序完成并行计算！
- 参考书：《高频面波方法》，夏江海，中国地质大学出版社
