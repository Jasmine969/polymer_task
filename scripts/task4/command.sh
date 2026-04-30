#!/usr/bin/env bash

# 任务四：纯水体系 production 轨迹分析命令。

# cd systems/water/ANALYSIS
cd systems/water_ethanol/ANALYSIS

# TPR = prod.tpr
TPR=prod_15ns.tpr

# 为聚合物单链建立索引组。后续 gmx gyrate 会使用这个索引组，
# 从而只计算聚合物本身的回转半径，而不把溶剂包含进去。
gmx_mpi select -s ../PROD/$TPR \
  -select '"CLS_SINGLE" atomnr 1 to 212' \
  -on cls.ndx

# 在几何分析前先处理周期性边界条件，使跨盒子的分子尽量保持完整。
# 这样可以减少聚合物被周期性边界切开而导致的分析误差。
printf "0\n" | gmx_mpi trjconv \
  -s ../PROD/$TPR \
  -f ../PROD/prod.xtc \
  -o prod_whole.xtc \
  -pbc whole

# 计算聚合物单链的回转半径 Rg(t)。
# 输出文件：systems/water/PROD/rg.xvg
printf "CLS_SINGLE\n" | gmx_mpi gyrate \
  -s ../PROD/$TPR \
  -f prod_whole.xtc \
  -n cls.ndx \
  -o rg.xvg

# 计算两个链端原子之间的末端距 Ree(t)。
# 这里主链一端的环碳（atom 2）到主链另一端的对应环碳（atom 195）作为链段
# 输出文件：systems/water/PROD/ree.xvg
gmx_mpi distance \
  -s ../PROD/$TPR \
  -f prod_whole.xtc \
  -select '"Ree" atomnr 2 plus atomnr 195' \
  -oall ree.xvg
