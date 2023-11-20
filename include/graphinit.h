#pragma once

#include "_readdata.h"

constexpr int V = 22;

constexpr int E = 1000;
// extern int h[V], e[E], ne[E], w[E];
// extern int cnt1;
// extern int vis[V];
// extern int dp[V][1ll << 21];
// 重新初始化Net数据
void Init_Graph();
// 进行建图
void Add_Edge(int u, int v, int w);
// 求解最小斯坦纳树
void min_stntree();
// 得到每一条边的实时权值
int Edge_Val(int curNum);
// 得到fpga之间的权值
int Wire_Val();
// dijkstra进行dp转移
void dijkstra(int s);
// 处理过程
void solve();
// 输出一组net的边
void printedge(int a, int b);
// 输出路径
void printnet_path(int netid);
// 计算延时
int calcdelay(int netid);
// 调整路径
void adjustpath();
// 初始化生成树
void Init_Tree(int netid);
// Die的数量
int Set_DieNum();
//分配wire信息（***方向）
void Assign_wire_info();
//输出布线结果
void Print_Layout_Res();
//输出TDM分配结果
void Print_Tdm_Res();
//并查集find函数
int find(int x);
int Find(int x);
//建立DFS树
void Build_Dfs_Tree(int netid, int u,int fa);
//计算每一个负载的delay
void calc_load_delay();




