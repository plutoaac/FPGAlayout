#include "graphinit.h"

#include <array>
#include <cstring>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>
#include <algorithm>

#include "_readdata.h"

// Prework
//----------------------------------------------------------------------

#define RG ReadDataSource ::getInstance()

// 全局最大延时
int Glo_Max_Delay;

// 每一组的最大延时

// Die的数量
int DieNum;

// 连边信息
int h[V], e[E], ne[E], w[E], fro[E];
// 向前星数组
int cnt1;
// 最短路的访问数组
int st[V];
// DP数组
int dp[V][1 << 21];
// wire之间的信号数量
int wiresig[V][V];

// 记录每一次斯坦纳树使用过的边
std::vector<int> line;

// Die之间的连线使用情况
int use[V][V];

//每一条路径的时延
std::unordered_map<int,std::unordered_map<int,double>>Net_Die_Path_delay;

// 每一个net的结点的路径
std::unordered_map<int, std::unordered_map<int, std::vector<int> > >Net_Die_Path;

//Wire的最大的netid 以及最大权值
std::unordered_map<std::pair<int,int>, std::unordered_map<int,double> >Net_on_Wire_Max_Val;

//每一个Wire上通过的net
std::unordered_map<std::pair<int, int>, std::vector<std::vector<int> >> Net_on_Wire;



// netid到斯坦纳树边的集合的映射
std::unordered_map<int, std::vector<int> > linemp;

// 输出路径的辅助数组
std::array<int, 3> pre[V][1 << 21];

// 将生成树进行重新构图
std::vector<int> G[V];

// BFS生成树使用的辅助数组
int vis2[V];

// dij转移使用的堆
std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >,
                    std::greater<std::pair<int, int> > >
    q;

// 每一组net的Die上的所包含的结点信息
std::unordered_set<int> Net_DieToNode[V];

//每一条跨fpga 所经过的netid的 max delay val的权值
struct node{
  int id;
  double  Max_Val{0.0};
};
std::vector<node>Net_on_Wire_SS;


//临时Wire的集合 （考虑方向） BFS树不可能有重边
std::vector<std::pair<int,int> >Tmp_Net_on_Wire_S;

//全局的Wire集合
std::unordered_set<std::pair<int,int> >Glo_Wire_Set;

//---------------------------------------------------------------------------

// 每一次处理Die时需要初始化
void Init_Graph() {
  cnt1 = 0;
  // memset(h, -1, sizeof(h));
  for (int i = 0; i < ReadDataSource::getInstance().DieNum; i++) {
    h[i] = -1;
  }
}

// 进行图的连接
void Add_Edge(int u, int v, int c) {
  fro[cnt1] = u;

  e[cnt1] = v;

  ne[cnt1] = h[u];

  w[cnt1] = c;

  h[u] = cnt1++;
}

// 如果当前可以使用的sll数量越大权值越小
//  权值越小斯坦纳树选到该边的概率越小
int Edge_Val(int CurNum) {
  return -1 * CurNum + ReadDataSource::getInstance().Max_SLL + 1;
}

// 记住改写
int Wire_Val() { return 1; }

int Set_DieNum() { DieNum = RG.DieNum; }

// 最短路进行dp转移
void dijkstra(int s) {
  for (int i = 0; i < ReadDataSource::getInstance().DieNum; i++) {
    st[i] = 0;
  }
  while (!q.empty()) {
    int u = q.top().second;
    q.pop();

    if (st[u]) continue;
    st[u] = 1;
    for (int i = h[u]; ~i; i = ne[i]) {
      int v = e[i];

      if (dp[v][s] > dp[u][s] + w[i]) {
        dp[v][s] = dp[u][s] + w[i];

        pre[v][s] = {1, u, i};

        q.push({dp[v][s], v});
      }
    }
  }
}
void printedge(int u, int A) {
  // DebuG
  //  std::cout << "op u  A" << pre[u][A][0] << u << " " << A << std::endl;

  if (pre[u][A][0] == 1) {
    line.push_back((pre[u][A][2]) >> 1);

    // 进行重新构图
    G[fro[pre[u][A][2]]].push_back(e[pre[u][A][2]]);

    G[e[pre[u][A][2]]].push_back(fro[pre[u][A][2]]);
    // Debug
    // std::cout << fro[pre[u][A][2]] << " " << e[pre[u][A][2]] << std::endl;

    printedge(pre[u][A][1], A);
  }

  if (pre[u][A][0] == 2) {
    printedge(u, pre[u][A][1]);

    printedge(u, pre[u][A][2]);
  }
}
void solve() {
  Set_DieNum();

  for (auto net : RG.Task) {
    Init_Graph();

    int netid = net.first;

    std::vector<int> key;

    int vis[20] = {0};

    for (auto node : net.second) {
      vis[RG.NodeToDie_Map[node]] = 1;
    }

    for (int i = 0; i < DieNum; ++i) {
      if (vis[i]) key.push_back(i);
    }

    // 对所有点进行连边

    for (int i = 0; i < DieNum; ++i) {
      for (int j = i + 1; j < DieNum; ++j) {
        if (RG.datatmp[i][j] == 0) continue;  // 不连通则跳出
        // 在同一个fpga
        if (RG.DieToFpga_Map[i] == RG.DieToFpga_Map[j]) {
          Add_Edge(i, j, Edge_Val(RG.datatmp[i][j]));

          Add_Edge(j, i, Edge_Val(RG.datatmp[i][j]));
          //  std::cout<<cnt1-1<<" "<<i<<" "<<j<<std::endl;
        }
        // 跨fpga
        else {
          Add_Edge(i, j, Wire_Val());

          Add_Edge(j, i, Wire_Val());
          //  std::cout << cnt1-1 << " " << i << " " << j << std::endl;
        }
      }
    }
    // 状态数量
    int mx = (1 << (key.size()));

    for (int i = 0; i < DieNum; i++) {
      for (int j = 0; j < mx; j++) {
        pre[i][j] = {0};
      }
    }

    for (int i = 0; i < DieNum; ++i) {
      for (int j = 0; j < mx; j++) {
        dp[i][j] = 1 << 30;
      }
    }
    for (int i = 0; i < key.size(); i++) {
      dp[key[i]][1 << i] = 0;
    }

    for (int sta = 1; sta < mx; sta++) {
      for (int i = 0; i < DieNum; i++) {
        for (int j = sta & (sta - 1); j; j = sta & (j - 1))

          if (dp[i][j] + dp[i][j ^ sta] < dp[i][sta]) {
            dp[i][sta] = dp[i][j] + dp[i][j ^ sta];

            pre[i][sta] = {2, j, j ^ sta};
          }

        if (dp[i][sta] != (1 << 30)) q.push({dp[i][sta], i});
      }

      dijkstra(sta);
    }
    // 打印该轮求解斯坦纳树使用的边
    printedge(key[0], mx - 1);

    // DeBug
    std::cout << "------------------" << std::endl;
    //  std::cout << dp[key[0]][mx - 1] << std::endl;
    //  std::cout << "-----------" << std::endl;

    for (auto &&it : line) {
      int a = fro[it];

      int b = e[it];

      if (RG.DieToFpga_Map[a] == RG.DieToFpga_Map[b]) {
        RG.datatmp[a][b]--, RG.datatmp[b][a]--;
        use[a][b]++;
        use[b][a]++;
      }
    }

    linemp[netid] = std::move(line);

    line.clear();

    printnet_path(netid);

    // 清空生成树连边信息
    for (int i = 0; i < DieNum; i++) {
      G[i].clear();
    }
    // Debug
    // return;
  }
}

void Init_Tree(int netid) {
  // Debug
  // std::cout << "begin" << std::endl;
  for (int i = 0; i < DieNum; i++) {
    Net_DieToNode[i].clear();
  }

  for (auto &&it : RG.Task[netid]) {
    // Debug
    // std::cout << RG.NodeToDie_Map[it] << " " << it << std::endl;

    Net_DieToNode[RG.NodeToDie_Map[it]].insert(it);
  }
  // Debug
  // std::cout << "end" << std::endl;

  for (int i = 0; i < DieNum; i++) {
    vis2[i] = 0;
  }
}

void printnet_path(int netid) {
  // Debug
  //  std::cout<<"netid"<<" "<<netid <<std::endl;
  //  std::cout<<"-------------------"<<std::endl;

  // BFS使用的队列
  std::queue<std::pair<int, std::vector<int> > > q;

  // 初始化BFStree
  Init_Tree(netid);

  int tmp = RG.NodeToDie_Map[RG.Task[netid][0]];

  // 处理属于同一个Die的情况
  for (auto &it : Net_DieToNode[tmp]) {
    if (it != tmp) {
    }
  }

  std::vector<int> v;
  v.push_back(tmp);

  // 将驱动节点所属于的Die扔进队列进行BFS
  q.push({tmp, v});

  vis2[tmp] = 1;
  while (!q.empty()) {
    int u = q.front().first;

    std::vector<int> ve = q.front().second;

    q.pop();

    for (auto v : G[u]) {
      if (!vis2[v]) {
        vis2[v] = 1; 

        ve.push_back(v);

        // 不在同一个FPGA
        if (RG.DieToFpga_Map[u] != RG.DieToFpga_Map[v]) {
        
           
          wiresig[u][v]++;
        }

        for (auto &it : Net_DieToNode[v]) {
          Net_Die_Path[netid][it] = ve;
        }
        q.push({v, ve});
      }
    }
  }
  // Debug
  std::unordered_map<int, std::vector<int> > tmp2;

  tmp2 = Net_Die_Path[netid];

  std::cout << "print node ----  path" << std::endl;
  for (auto &&it : tmp2) {
    std::cout << it.first << "----";

    for (auto &&v : it.second) {
      std::cout << v << " ";
    }

    std::cout << std::endl;
  }
  std::cout << "-----------------------------------------" << std::endl;
}

// 后续需要改
void Assign_wire_info() {
  // for (auto &&it : RG.Not_Fpga_Die) {
  //   int a = it.first;
  //   int b = it.second;
  //   // double tt= wiresig[a][b]*1.0/wiresig[b][a];
  //   int ass_0 = RG.data[a][b] *
  //               (1.0 * (wiresig[a][b]) / (wiresig[a][b] + wiresig[b][a]));
  //   int ass_1 = RG.data[a][b] - ass_0;

  //   // int p=
  // }

   
  //将一条wire上的所经过的net的 max vir val 排序    Glo

  for(auto &&it:Glo_Wire_Set){
      
      //wire的两端  Die的编号
      int a= it.first;
      int b= it.second;
   
     for(auto &&it2:  Net_on_Wire_Max_Val [{a,b}]){
          
          //netid和其对应的最大的权值 Glo
          int Netid=it2.first;
          double val=it2.second; 

          Net_on_Wire_SS.push_back(node{Netid,val});
     }
  }

  //将全局的每一条Wire的最大权值进行排序
  sort(Net_on_Wire_SS.begin(),Net_on_Wire_SS.end(),
   [](node a,node b)->bool{
        return a.Max_Val>b.Max_Val;
   } 
  );



  //1.优先队列合并
  //2.线段树如何合并

  for(auto &&it:Glo_Wire_Set){
      
      //wire的两端  Die的编号
      int a= it.first;

      int b= it.second;
    
      std::cout<<"---"<<std::endl;
  }  
  

}

int calcdelay(int netid) {
  //下标从1开始  0为驱动结点
  for (int i = 1; i < RG.Task[netid].size(); i++) {
    //拿到负载节点
    int id = RG.Task[netid][i];
    // int fi = *Net_Die_Path[netid][id].begin();
    //初始化每一组nie的结点的延时为0
    double Delay_ans=0;
    for (int j = 1; j < Net_Die_Path[netid][id].size(); j++) {
      
      int next2 = Net_Die_Path[netid][id][j];
      
      int next1 = Net_Die_Path[netid][id][j - 1];
      
      if (RG.DieToFpga_Map[next1] != RG.DieToFpga_Map[next2]) {
        // 本次处理net所经过的跨越FPGA的Die 的集合
        Tmp_Net_on_Wire_S.push_back({next1, next2});

        //将该条Wire扔进全局的Wire集合
        Glo_Wire_Set.insert(std::make_pair(next1,next2));  

        Net_on_Wire[{next1, next2}][netid].push_back(id);

        // （2*x+1）/2   val% 4.5==0   首先按照只消耗了最低的delay考虑
        // 后续继续动态调整
        Delay_ans += 4.5;  
      
      }else{

        Delay_ans++;
      
      }

    }

    Net_Die_Path_delay[netid][id]=Delay_ans;
    
    for(auto &&it:Tmp_Net_on_Wire_S){
      
      int a=it.first;

      int b=it.second;
      
      Net_on_Wire_Max_Val[{a,b}][netid]=std::max(Delay_ans, Net_on_Wire_Max_Val[{a,b}][netid]);

   
    }

    //每次处理一条驱动节点到负载节点 清空临时的Wire集合
    Tmp_Net_on_Wire_S.clear();
  }

}

void adjustpath() {







}
