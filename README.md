# EIS_Toolbox: A Physically-Informed Automated Robust Framework

EIS_Toolbox is a high-throughput automated analysis framework designed for Electrochemical Impedance Spectroscopy (EIS) data. It provides a robust solution for the batch analysis of impedance spectra.

## Applicability
Strictly applicable to polymer solid electrolyte blocking systems and similar ionic conductor systems.

## Algorithmic Robustness
This framework is engineered to screen and isolate non-ideal polarization scenarios encountered in experimental environments:
- Strong Random Noise: Identified and suppressed via effective smoothing operators.
- High-Frequency Response Truncation: Managed by topology recognition mechanisms.
- Low-Frequency Non-Stationary Drift: Isolated via segmented governance strategies and asymmetric geometric truncation to ensure analysis on valid data segments.

## Key Features
- Robust Data Screening: Multi-level adaptive preprocessing based on Kramers-Kronig (K-K) consistency verification to isolate non-intrinsic signals.
- Smart Topology Recognition: Hierarchical recognition integrating geometric priors and statistical decision-making.
- High-Throughput Inversion: Automated thermodynamic and kinetic analysis extending to temperature-dependent sequences.
- User-Friendly: Includes both a Python CLI and a GUI application (EIS_Toolbox.exe).

## Prerequisites
This project requires Python 3.8+ and the following dependencies:
pip install -r requirements.txt

## Quick Start
Running the Application

**GUI Version :** No installation required. Click the links below to download directly, and double-click to run:
- [Download English Version (EIS_Toolbox_v1.5_EN.exe)](https://github.com/ChemCode-Chang/EIS-Toolbox/releases/download/v1.5/EIS_Toolbox_v1.5_EN.exe)
- [Download Chinese Version (EIS_Toolbox_v1.5_CN.exe)](https://github.com/ChemCode-Chang/EIS-Toolbox/releases/download/v1.5/EIS_Toolbox_v1.5_CN.exe)

- Script Version: Run the main program via terminal: `python gui_app.py`

## Citation
(Citation details will be added upon formal paper publication.)

## License
This project is licensed under the MIT License.

## Contact
For questions or feedback, please contact: 2014488020@qq.com

----------------------------------------------------------------------------------------------------------------
# EIS_Toolbox: 高通量、物理依从的自动化 EIS 分析框架

EIS_Toolbox 是一个为电化学阻抗谱 (EIS) 数据设计的高通量自动化分析框架，旨在为阻抗谱数据提供鲁棒的批处理分析方案。

## 适用范围
仅适用于聚合物固态电解质阻塞体系及其类似离子导体体系。

## 算法鲁棒性
本框架针对实验中常见的非理想极化场景进行了深度优化，具备极强的鲁棒性，能够客观地识别并处理干扰信号：
- 强随机噪声： 通过有效的滤波算子进行筛选与抑制。
- 高频响应截断：通过拓扑识别机制进行特征识别。
- 低频非稳态漂移： 采用分段治理策略与非对称几何截断，对漂移数据进行隔离，确保仅在有效数据区间内进行分析。

## 主要功能
- 鲁棒的数据筛选： 基于 Kramers-Kronig (K-K) 一致性校验的多级自适应预处理，用于甄别异常信号。
- 智能拓扑识别：融合几何先验与统计判决的分层识别机制，消除了人工拟合偏差。
- 高通量反演： 支持从单点分析扩展至全温度序列的热力学与动力学自动化计算。
- 用户友好： 提供 Python 脚本接口以及配套的 GUI 桌面程序 (EIS_Toolbox.exe)。

## 安装说明
本项目需要 Python 3.8+ 环境，请运行以下命令安装依赖：
pip install -r requirements.txt

## 快速开始
如何运行程序
**GUI 软件版:** 无需配置 Python 环境，直接点击下方链接下载，双击即可运行：
- [点击下载 中文版 (EIS_Toolbox_v1.5_CN.exe)](https://github.com/ChemCode-Chang/EIS-Toolbox/releases/download/v1.5/EIS_Toolbox_v1.5_CN.exe)
- [点击下载 英文版 (EIS_Toolbox_v1.5_EN.exe)](https://github.com/ChemCode-Chang/EIS-Toolbox/releases/download/v1.5/EIS_Toolbox_v1.5_EN.exe)

- 源码脚本版: 在终端中运行主程序: `python gui_app.py`

## 引用说明
(本软件相关论文正在审稿中，引用信息将在论文正式发表后更新。)

## 开源协议
本项目采用 MIT 开源协议。

## 联系方式
如有问题或建议，请联系：2014488020@qq.com
