# polar-coded-SCMA
***注：本程序只供学习交流使用，请勿用于商业目的。***

## 程序说明

1. 程序实现了polar码译码和SCMA多用户检测的联合检测译码问题，其中polar码译码采用SCAN软译码算法，SCMA采用MPA算法
2. 主函数包括四个子函数
   - initPC：用于polar码编译码的初始化工作   
   - pencode：polar码编码
   - scmaenc：SCMA编码映射
   - JIDD：联合检测译码主函数
     - SCMA检测--->polar码译码(SCAN算法)--->SCMA检测--->polar码译码(SCAN算法)-->......
   - 参数: 
     - polar_N: polar码码长
     - polar_K: polar码信息位长度
     - alpha：论文当中的$\alpha$ 因子
     - isInterleaver: 是否增加交织器
3. constructedCode文件夹下为构造的polar码码字
4. 仿真结果：

![image](https://github.com/ZhipengPan/polar-coded-SCMA/blob/master/result.bmp)

## 参考文献

[1] Z. Pan, E. Li, L. Wen, J. Lei, and C. Tang, “Joint iterative detection and decoding receiver for polar coded SCMA system,” in 2018 IEEE International Conference on Communications Workshops (ICC Workshops), May 2018, pp. 1–6.

[2] Z. Pan, E. Li, L. Zhang, J. Lei, and C. J. Tang, “Design and optimization of joint iterative detection and decoding receiver for uplink polar coded scma system,” IEEE Access, pp. 1–1, 2018.    



## 联系作者

潘志鹏

湖南，长沙

邮件：panzhipeng10@nudt.edu.cn

