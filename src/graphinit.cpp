#include "graphinit.h"

#include <assert.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>

#include "_readdata.h"

// Prework
//----------------------------------------------------------------------

#define RG ReadDataSource ::getInstance()

#define __TIMIT 3

using ll = long long;

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

// 全局最大延时
double Glo_Max_Delay;

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
//  SLL到达瓶颈标志
int flag{0};
int flag1{0};

// 每一个net的结点的路径
std::unordered_map<int, std::unordered_map<int, std::vector<int>>> Net_Die_Path;

// 输出路径的辅助数组
std::array<int, 3> pre[V][1 << 21];

// 将生成树进行重新构图
std::vector<int> G[V];

// BFS生成树使用的辅助数组
int vis2[V];

// dij转移使用的堆
std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                    std::greater<std::pair<int, int>>>
    q;

// 每一组net的Die上的所包含的结点信息
std::unordered_map<int, std::unordered_map<int, std::vector<int>>>
    Net_DieToNode;

// 临时Wire的集合 （考虑方向） BFS树不可能有重边
std::vector<std::pair<int, int>> Tmp_Net_on_Wire_S;

// 全局的Wire集合
std::vector<std::pair<int, int>> Glo_Wire_Set;

// 每一组Net的最大时延
std::vector<std::pair<int, double>> Net_Max_Delay;

// Wire的使用结果
std::unordered_map<std::pair<int, int>, std::map<int, std::vector<int>>,
                   pair_hash>
    Wire_Res;

// 存储临时的超限的信号
std::vector<std::pair<int, double>> Tmp_S;

// 并查集的fa结点
std::unordered_map<int, int> fa;

// 最小生成树存边
std::vector<std::array<int, 3>> ee;

// kru的并查集
std::unordered_map<int, int> F;

// DFS序 使用的数据
std::unordered_map<int, std::map<int, int>> l;
std::unordered_map<int, std::map<int, int>> r;
std::unordered_map<int, std::map<int, int>> id;

int idx{0};

// 每一组net的Die的权值
std::unordered_map<int, std::unordered_map<int, double>> Net_Die_Val;

// 每一条wire 所经过的net id 集合
std::unordered_map<std::pair<int, int>, std::vector<int>, pair_hash> Wire_Net_S;

std::vector<int> G1[22];

// 存储优先队列里的结点
struct pri_node {
  double ma{0};
  int zuhao;
  int sz;
};

struct CompareDelay {
  bool operator()(const pri_node &node1, const pri_node &node2) {
    // 按照 delay 从小到大排序
    return node1.ma > node2.ma;
  }
};

//---------------------------------------------------------------------------

// 每一个net需要维护一个线段树
struct SGT {
  // 线段树的node结点
  struct node {
    int l, r;
    double sum{0};
    double maxv{0};
#define ls i * 2
#define rs i * 2 + 1
  } tr[22 * 4];
  // 需要改
  std::map<int, double> a;
  // push操作维护
  void pushup(int i) {
    tr[i].sum = tr[ls].sum + tr[rs].sum;
    tr[i].maxv = std::max(tr[ls].maxv, tr[rs].maxv);
  }

  void build(int i, int l, int r) {
    tr[i] = {l, r, 0.0, 0.0};
    if (l == r) {
      tr[i].sum = a[l];
      return;
    }
    int mid = l + r >> 1;
    build(ls, l, mid);
    build(rs, mid + 1, r);
    pushup(i);
  }
  // 区间修改
  void change(int i, int l, int r, double val) {
    if (tr[i].l == tr[i].r) {
      tr[i].sum += val;
      tr[i].maxv += val;
      return;
    }
    int mid = (tr[i].l + tr[i].r) >> 1;
    if (l <= mid) change(ls, l, r, val);
    if (r > mid) change(rs, l, r, val);
    pushup(i);
  }
  // 查询区间和
  double query_sum(int i, int x) {
    if (tr[i].l == tr[i].r) return tr[i].sum;
    int mid = (tr[i].l + tr[i].r) >> 1;
    if (x <= mid)
      return query_sum(ls, x);
    else
      return query_sum(rs, x);
  }

  // 查询区间max
  double query_max(int i, int l, int r) {
    if (l <= tr[i].l && tr[i].r <= r) return tr[i].sum;
    double mx = 0.0;
    int mid = (tr[i].l + tr[i].r) >> 1;
    if (l <= mid) mx = std::max(mx, query_max(ls, l, r));
    if (r > mid) mx = std::max(mx, query_max(rs, l, r));
    return mx;
  }
};

//-------------------------------------------------------------------------

// 存储每一个线段树
std::map<int, SGT> Sgt_Trees;

// 合并TDM结束后的每一组    pair存储他是正向的还是反向的
std::unordered_map<int, std::vector<std::pair<int, int>>> Merge_S;

// 每个wire

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
  // Debug

  if (flag <= __TIMIT || !flag1) return 8;

  return 1;
  // return -1 * CurNum + ReadDataSource::getInstance().Max_SLL + 2;
}

// 记住改写
int Wire_Val(int CurNum) {
  if (flag <= __TIMIT || !flag1) return 1;
  // if (flag <= __TIMIT || !flag1) return 20340 - static_cast<int>(CurNum) *
  //  0.2; 0.6;
  //   0.6;

  return 4;
}

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
int Find(int x) { return x == F[x] ? x : F[x] = Find(F[x]); }
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
    // std::cout << "-------------------" << std::endl;

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

    // std::cout<<netid<<std::endl;

    std::vector<int> key;

    int vis[22] = {0};

    for (auto node : net.second) {
      vis[RG.NodeToDie_Map[node.first]] = 1;
    }

    for (int i = 0; i < DieNum; ++i) {
      if (vis[i]) key.push_back(i);
    }
    if (key.size() == DieNum) {
      /////
      std::cout << "kru" << std::endl;
      ee.clear();
      for (int i = 0; i < DieNum; i++) F[i] = i;
      for (int i = 0; i < DieNum; i++) {
        for (int j = i + 1; j < DieNum; j++) {
          if (RG.data[i][j] == 0) continue;
          if (use[i][j] == RG.data[i][j]) continue;
          if (RG.DieToFpga_Map[i] == RG.DieToFpga_Map[j])
            ee.push_back({Edge_Val(use[i][j]), i, j});
          else {
            ee.push_back({Wire_Val(use[i][j]), i, j});
          }
        }
      }
      sort(ee.begin(), ee.end());
      for (int i = 0; i < ee.size(); i++) {
        int w = ee[i][0];
        int u = ee[i][1];
        int v = ee[i][2];
        int a = Find(u);
        int b = Find(v);

        if (a == b)
          continue;
        else {
          G[u].push_back(v);
          G[v].push_back(u);
          F[a] = b;
          // use[u][v]++;
          // use[v][u]++;
        }
      }
    }
    // else if(key.size()>=RG.DieNum/2){
    //     //斯坦纳森林
    //     //先拿一半的点去做斯坦纳
    //     // 不一定最优
    // }

    else {
      // 对所有点进行连边

      for (int i = 0; i < DieNum; ++i) {
        for (int j = i + 1; j < DieNum; ++j) {
          if (RG.data[i][j] == 0) continue;  // 不连通则跳出
          // 在同一个fpga
          assert(use[i][j] == use[j][i]);
          if (use[i][j] == RG.data[i][j]) continue;
          if (RG.DieToFpga_Map[i] == RG.DieToFpga_Map[j]) {
            Add_Edge(i, j, Edge_Val(use[i][j]));
            Add_Edge(j, i, Edge_Val(use[j][i]));
            //  std::cout<<cnt1-1<<" "<<i<<" "<<j<<std::endl;
          }
          // 跨fpga
          else {
            Add_Edge(i, j, Wire_Val(use[i][j]));

            Add_Edge(j, i, Wire_Val(use[j][i]));
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
      // std::cout << "------------------" << std::endl;
      //  std::cout << dp[key[0]][mx - 1] << std::endl;
      //  std::cout << "-----------" << std::endl;
    }

    printnet_path(netid);
    calcdelay(netid);
    Net_Die_Val[netid][RG.NodeToDie_Map[net.second[0].first]] = 0;

    SGT tmp;
    Sgt_Trees[netid] = std::move(tmp);

    idx = 0;
    Build_Dfs_Tree(netid, RG.NodeToDie_Map[net.second[0].first], -1);

    Sgt_Trees[netid].build(1, 1, idx);

    // 清空生成树连边信息
    for (int i = 0; i < DieNum; i++) {
      G[i].clear();
    }
    cnt1 = 0;

    // Debug
    // return;
  }

  Assign_wire_info();
  calc_load_delay();
  adjustpath();
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

    // 第一个结点是驱动结点 我们这里只存储 每一组net的Die上所包含的负载结点
    if (it.first == RG.Task[netid][0].first) continue;

    Net_DieToNode[netid][RG.NodeToDie_Map[it.first]].push_back(it.first);
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
  std::queue<std::pair<int, std::vector<int>>> q;

  // 初始化BFStree
  Init_Tree(netid);

  // 驱动结点所属于的Die
  int tmp = RG.NodeToDie_Map[RG.Task[netid][0].first];

  // 处理属于同一个Die的情况
  for (auto &it : Net_DieToNode[netid][tmp]) {
    // 如果在同一个Die里面 那么他的路径只有驱动结点所在的Die   Delay为0

    Net_Die_Path[netid][it].push_back(tmp);
  }

  std::vector<int> v;
  // Debug  输出驱动节点的DIE
  // std::cout << tmp << std::endl;
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
        std::vector<int> ve1 = ve;
        ve1.push_back(v);

        // 不在同一个FPGA
        if (RG.DieToFpga_Map[u] != RG.DieToFpga_Map[v]) {
          Wire_Net_S[{u, v}].push_back(netid);
          use[u][v]++;
          use[v][u]++;
          if (use[u][v] >= RG.data[u][v] * 0.5) {
            flag++;
          }
          if (use[u][v] >= RG.data[u][v] * 0.85) {
            flag1++;
          }
        } else {
          use[v][u]++;
          use[u][v]++;
        }

        for (auto &it : Net_DieToNode[netid][v]) {
          Net_Die_Path[netid][it] = ve1;
        }
        q.push({v, ve1});
      }
    }
  }

  //  Debug
  // 输出 负载节点-- 路径 所经过的结点
  std::unordered_map<int, std::vector<int>> tmp2;

  std::cout << netid << std::endl;
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

void Build_Dfs_Tree(int netid, int u, int fa) {
  id[netid][u] = ++idx;
  l[netid][u] = idx;
  double tmp = Net_Die_Val[netid][u];
  Sgt_Trees[netid].a[idx] = tmp;
  for (auto v : G[u]) {
    if (v == fa) continue;
    Build_Dfs_Tree(netid, v, u);
  }
  r[netid][u] = idx;
}

// 并查集find函数
int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }

void Assign_wire_info() {
  // RG.Not_Fpga_Die.clear();
  // RG.Not_Fpga_Die.push_back({3,7});
  // RG.Not_Fpga_Die.push_back({0,4});

  for (auto &it : RG.Not_Fpga_Die) {
    int a = it.first;
    int b = it.second;

    // a,b 表示 wire的两端  Die的编号
    if (Wire_Net_S[{a, b}].size() + Wire_Net_S[{b, a}].size() == 0) continue;

    int all = (Wire_Net_S[{a, b}].size() + 3) / 4 +
              (Wire_Net_S[{b, a}].size() + 3) / 4;

    if (all <= RG.data[a][b]) {
      int h = 0;
      for (const int netid : Wire_Net_S[{a, b}]) {
        Wire_Res[{a, b}][(h++) / 4 + 1].push_back(netid);
      }

      // 不能紧接着 因为 不同方向不能再同一条Wire
      if (h % 4) {
        h = ((h - 1) / 4 + 1) * 4;
      }
      for (const int netid : Wire_Net_S[{b, a}]) {
        Wire_Res[{a, b}][(h++) / 4 + 1].push_back(netid);
      }

      std::ofstream outputFile("./design.tdm.out", std::ios::app);
      if (outputFile.is_open()) {
        outputFile << "[Die" << a << ",Die" << b << ']' << std::endl;
        for (int i = 1; i <= all; i++) {
          outputFile << '[';
          for (auto x : Wire_Res[{a, b}][i]) {
            if (x != Wire_Res[{a, b}][i].back())
              outputFile << x << ',';
            else
              outputFile << x << ']';
          }
          outputFile << " 4" << std::endl;
        }

        outputFile.close();
      } else {
        // 错误处理
        std::cerr << "无法打开文件" << std::endl;
      }

      continue;
    }
    std::cout << "mergebegin" << std::endl;

    std::vector<std::pair<int, double>> x1;
    for (const int &netid : Wire_Net_S[{a, b}]) {
      x1.push_back(
          {netid, Sgt_Trees[netid].query_max(1, l[netid][b], r[netid][b])});
    }
    sort(x1.begin(), x1.end(),
         [](std::pair<int, double> a, std::pair<int, double> b) -> bool {
           return a.second > b.second;
         });

    std::vector<std::pair<int, double>> x2;
    for (const int &netid : Wire_Net_S[{b, a}]) {
      x2.push_back(
          {netid, Sgt_Trees[netid].query_max(1, l[netid][a], r[netid][a])});
    }
    sort(x2.begin(), x2.end(),
         [](std::pair<int, double> a, std::pair<int, double> b) -> bool {
           return a.second > b.second;
         });
    // return ;
    //  时延 组号 大小
    std::priority_queue<pri_node, std::vector<pri_node>, CompareDelay> Merge_q1;

    std::priority_queue<pri_node, std::vector<pri_node>, CompareDelay> Merge_q2;

    // 初始化并查集结点的fa
    for (int i = 1; i <= all; i++) {
      fa[i] = i;
    }

    int h{0};

    int idd{0};

    for (const std::pair<int, double> &itt : x1) {
      if (h % 4 == 0) {
        Merge_q1.push({itt.second, ++idd, 1});
      }
      Wire_Res[{a, b}][(h++) / 4 + 1].push_back(itt.first);
    }

    int zheng = idd;

    h = idd * 4;

    for (const std::pair<int, double> &itt : x2) {
      if (h % 4 == 0) {
        Merge_q2.push({itt.second, ++idd, 1});
      }
      Wire_Res[{b, a}][(h++) / 4 + 1].push_back(itt.first);
    }
    // DEbug
    // std::cout << Merge_q1.size() << std::endl;
    // std::cout << Merge_q2.size() << std::endl;

    int Merge_cnt = (Wire_Net_S[{a, b}].size() + 3) / 4 +
                    (Wire_Net_S[{b, a}].size() + 3) / 4 - RG.data[a][b];
    // DEbug
    // std::cout << Merge_cnt << std::endl;

    int cur_cnt{0};

    bool flag1 = 0, flag2 = 0;
    pri_node t1, t2;
    pri_node t1_, t2_;
    double val1, val1_, val2, val2_;
    int id1, id1_, id2, id2_;
    int sz1, sz1_, sz2, sz2_;
    cur_cnt = 0;

    // std::cout << "size" << Merge_q1.size() << " " << Merge_q2.size()
    //           << std::endl;
    while ((Merge_q1.size() >= 2 || flag1 == 1) &&
           (Merge_q2.size() >= 2 || flag2 == 1)) {
      if (flag1 == 0) {
        t1 = Merge_q1.top();
        Merge_q1.pop();

        val1 = t1.ma;
        id1 = t1.zuhao;
        sz1 = t1.sz;

        t1_ = Merge_q1.top();
        Merge_q1.pop();

        val1_ = t1_.ma;
        id1_ = t1_.zuhao;
        sz1_ = t1_.sz;
        flag1 = 1;
      }
      //--------------------------------------
      if (flag2 == 0) {
        t2 = Merge_q2.top();
        Merge_q2.pop();

        val2 = t2.ma;
        id2 = t2.zuhao;
        sz2 = t2.sz;

        t2_ = Merge_q2.top();
        Merge_q2.pop();

        val2_ = t2_.ma;
        id2_ = t2_.zuhao;
        sz2_ = t2_.sz;
        flag2 = 1;
      }

      if (val1_ + sz1 * 4.0 < val2_ + sz2 * 4.0) {
        // 左合并
        int a = find(id1);
        int b = find(id1_);
        fa[a] = b;
        double tmp_val = val1_ + sz1 * 4.0;

        Merge_q1.push({tmp_val, b, sz1 + sz1_});
        flag1 = 0;
      } else {
        // 右合并
        int a = find(id2);
        int b = find(id2_);
        fa[a] = b;
        double tmp_val = val2_ + sz2 * 4.0;

        Merge_q2.push({tmp_val, b, sz2 + sz2_});

        flag2 = 0;
      }
      //++
      cur_cnt++;
      if (cur_cnt == Merge_cnt) break;
    }

    if (cur_cnt < Merge_cnt) {
      while (Merge_q1.size() >= 2 || flag1 == 1) {
        if (flag1 == 0) {
          t1 = Merge_q1.top();
          Merge_q1.pop();

          val1 = t1.ma;
          id1 = t1.zuhao;
          sz1 = t1.sz;

          t1_ = Merge_q1.top();
          Merge_q1.pop();

          val1_ = t1_.ma;
          id1_ = t1_.zuhao;
          sz1_ = t1_.sz;
          flag1 = 1;
        }

        // 左合并
        int a = find(id1);
        int b = find(id1_);
        fa[a] = b;
        double tmp_val = val1_ + sz1 * 4.0;

        Merge_q1.push({tmp_val, b, sz1 + sz1_});
        flag1 = 0;
        cur_cnt++;
        if (cur_cnt == Merge_cnt) break;
      }

      while (Merge_q2.size() >= 2 || flag2 == 1) {
        if (flag2 == 0) {
          t2 = Merge_q2.top();
          Merge_q2.pop();

          val2 = t2.ma;
          id2 = t2.zuhao;
          sz2 = t2.sz;

          t2_ = Merge_q2.top();
          Merge_q2.pop();

          val2_ = t2_.ma;
          id2_ = t2_.zuhao;
          sz2_ = t2_.sz;
          flag2 = 1;
        }
        int a = find(id2);
        int b = find(id2_);
        fa[a] = b;
        double tmp_val = val2_ + sz2 * 4.0;

        Merge_q2.push({tmp_val, b, sz2 + sz2_});
        // flag2 = 0;
        // if (cur_cnt == Merge_cnt) break;
        // 右合并
        flag2 = 0;
        cur_cnt++;
        if (cur_cnt == Merge_cnt) break;
      }
    } else {
      if (flag1 != 0) {
        // 左扔回去
        Merge_q1.push(t1);
        Merge_q1.push(t1_);
      } else {
        // 右扔回去
        Merge_q2.push(t2);
        Merge_q2.push(t2_);
      }
    }
    // std::cout << "kaishi" << std::endl;
    // std::cout << Merge_q1.size() + Merge_q2.size() << std::endl;
    while (!Merge_q1.empty()) {
      auto t = Merge_q1.top();
      //  std::cout << t.sz << std::endl;
      Merge_q1.pop();
    }

    while (!Merge_q2.empty()) {
      auto t = Merge_q2.top();
      // std::cout << t.sz << std::endl;
      Merge_q2.pop();
    }
    for (int i = 1; i <= idd; i++) {
      Merge_S[i].clear();
    }
    for (int i = 1; i <= idd; i++) {
      int a = find(i);
      if (i <= zheng)
        Merge_S[a].push_back({i, 0});
      else
        Merge_S[a].push_back({i, 1});
    }

    for (int i = 1; i <= idd; i++) {
      int t = find(i);
      if (i == t) {
        std::vector<int> anspri;

        int ff = 0;
        for (int i = 0; i < Merge_S[t].size(); i++) {
          auto x = Merge_S[t][i];
          int zuhao = x.first;

          if (x.second == 0) {
            ff = 1;
            for (int j = 0; j < Wire_Res[{a, b}][zuhao].size(); j++) {
              // std::cout << Wire_Res[{a, b}][zuhao][j] << std::endl;
              anspri.push_back(Wire_Res[{a, b}][zuhao][j]);
            }
          } else {
            for (int j = 0; j < Wire_Res[{b, a}][zuhao].size(); j++) {
              anspri.push_back(Wire_Res[{b, a}][zuhao][j]);
            }
          }
        }
        int extra_val = (anspri.size() + 3) / 4 - 1;
        if (!ff) {
          for (auto x : anspri) {
            Sgt_Trees[x].change(1, l[x][b], r[x][b], extra_val * 4);
          }
        } else {
          for (auto x : anspri)
            Sgt_Trees[x].change(1, l[x][a], r[x][a], extra_val * 4);
        }
      }
    }

    //--------------------
    // 输出TDM Ratio

    /*std::cout << idd << std::endl;
    for (int i = 1; i <= idd; i++) {
      if (i == find(i)) std::cout << "i" <<find(i) << std::endl;
    }*/
    std::ofstream outputFile("./design.tdm.out", std::ios::app);

    if (outputFile.is_open()) {
      outputFile << "[Die" << a << ",Die" << b << ']' << std::endl;

      for (int i = 1; i <= idd; i++) {
        int t = find(i);
        if (i == t) {
          std::vector<int> anspri;
          std::vector<int> anspri1;
          for (int i = 0; i < Merge_S[t].size(); i++) {
            auto x = Merge_S[t][i];
            int zuhao = x.first;

            if (x.second == 0) {
              for (int j = 0; j < Wire_Res[{a, b}][zuhao].size(); j++) {
                // std::cout << Wire_Res[{a, b}][zuhao][j] << std::endl;
                anspri.push_back(Wire_Res[{a, b}][zuhao][j]);
              }
            } else {
              for (int j = 0; j < Wire_Res[{b, a}][zuhao].size(); j++) {
                anspri1.push_back(Wire_Res[{b, a}][zuhao][j]);
              }
            }
          }
          outputFile << '[';

          for (int j = 0; j < anspri1.size(); j++) {
            if (j == anspri1.size() - 1) {
              outputFile << anspri1[j];
            } else {
              outputFile << anspri1[j] << ",";
            }
          }

          for (int j = 0; j < anspri.size(); j++) {
            if (j == anspri.size() - 1) {
              outputFile << anspri[j];
            } else {
              outputFile << anspri[j] << ",";
            }
          }
          if (4 * ((anspri.size() + 3) / 4))
            outputFile << "]zheng" << 4 * ((anspri.size() + 3) / 4)
                       << std::endl;
          else {
            outputFile << "]fan" << 4 * ((anspri1.size() + 3) / 4) << std::endl;
          }
        }
      }
      outputFile.close();
    } else {
      // 错误处理
      std::cerr << "无法打开文件" << std::endl;
    }
  }
}

int calcdelay(int netid) {
  // 下标从1开始  0为驱动结点
  for (int i = 1; i < RG.Task[netid].size(); i++) {
    // 负载节点
    int id = RG.Task[netid][i].first;
    // int fi = *Net_Die_Path[netid][id].begin();
    // 初始化每一组nie的结点的延时为0
    double Delay_ans = 0;

    // 遍历驱动节点到负载结点的路径  存的是结点信息 所经过的结点
    for (int j = 1; j < Net_Die_Path[netid][id].size(); j++) {
      // 当前结点
      int next2 = Net_Die_Path[netid][id][j];
      // 前一个结点
      int next1 = Net_Die_Path[netid][id][j - 1];

      // 如果不在同一个FPGA
      if (RG.DieToFpga_Map[next1] != RG.DieToFpga_Map[next2]) {
        // 将该条Wire扔进全局的Wire集合
        Glo_Wire_Set.push_back({next1, next2});

        // （2*x+1）/2   val% 4.5==0   首先按照只消耗了最低的delay考虑
        // 后续继续动态调整
        Delay_ans += 4.5;
      } else {
        Delay_ans += 1.0;
      }
    }
    //  std::cout << Delay_ans << std::endl;
    RG.Task[netid][i].second = Delay_ans;
    Net_Die_Val[netid][RG.NodeToDie_Map[id]] = Delay_ans;
  }
}

void adjustpath() {
  // Debug
  for (int i = 0; i < RG.DieNum; i++) {
    for (int j = 0; j < RG.DieNum; j++) {
      std::cout << use[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "*******************************" << std::endl;
  for (int i = 0; i < RG.DieNum; i++) {
    for (int j = 0; j < RG.DieNum; j++) {
      std::cout << RG.data[i][j] << " ";
    }
    std::cout << std::endl;
  }
  for (int i = 0; i < RG.DieNum; i++) {
    for (int j = i + 1; j < RG.DieNum; j++) {
      if (RG.DieToFpga_Map[i] != RG.DieToFpga_Map[j]) continue;
      if (use[i][j] > RG.data[i][j]) {
        std::cout << "#" << i << " " << j << " " << use[i][j] - RG.data[i][j]
                  << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

void calc_load_delay() {
  for (auto &&it : RG.Task) {
    for (int i = 1; i < it.second.size(); i++) {
      int loader = it.second[i].first;
      int tmp = id[it.first][RG.NodeToDie_Map[loader]];
      double curdelay = Sgt_Trees[it.first].query_sum(1, tmp);
      it.second[i].second = curdelay;
    }
  }
}

void Print_Layout_Res() {
  std::ofstream outputFile("./design.route.out");

  if (outputFile.is_open()) {
    std::unordered_map<int, std::vector<std::pair<int, double>>> Task1 =
        RG.Task;
    for (auto &&it : RG.Task) {
      // 输出netid
      // outputFile << "[" << it.first << "]" << std::endl;

      std::vector<std::pair<int, double>> outve = it.second;

      outve.erase(outve.begin());

      sort(outve.begin(), outve.end(),
           // 时延排序的的lambda函数 从大到小
           [](std::pair<int, double> a, std::pair<int, double> b) -> bool {
             return a.second > b.second;
           });

      // 保存一下使用过的netid
      Net_Max_Delay.push_back({it.first, outve[0].second});

      // 重新移动到Task  已经排好序 并且已经去点了驱动节点
      it.second = std::move(outve);
    }
    // 按照所有net的最大值排序
    sort(Net_Max_Delay.begin(), Net_Max_Delay.end(),

         // 每组net的时延排序的函数 从大到小
         [](std::pair<int, double> a, std::pair<int, double> b) -> bool {
           return a.second > b.second;
         });

    Glo_Max_Delay = std::max(Glo_Max_Delay, Net_Max_Delay[0].second);

    //------------------------------------------------
    // it.first为netid
    for (const std::pair<int, double> &it : Net_Max_Delay) {
      // 输出netid
      int netid = it.first;
      outputFile << "[" << netid << "]" << std::endl;
      for (int i = 1; i < Task1[netid].size(); i++) {
        const std::pair<int, int> &ii = Task1[netid][i];

        outputFile << '[';
        // Net_Die_Path[it][id]
        int nodetmp;
        int loader = ii.first;
          
        if (Net_Die_Path[netid][loader].back() !=
            *Net_Die_Path[netid][loader].begin()) {
          for (const auto node : Net_Die_Path[netid][loader]) {
            if (node != Net_Die_Path[it.first][ii.first].back())
              outputFile << node << ",";
            else {
              nodetmp = node;
              outputFile << node << ']';
            }
          }
        } else {
          nodetmp = Net_Die_Path[netid][ii.first].back();
          outputFile << Net_Die_Path[netid][ii.first].back() << ']';
        }

        //-----------
        outputFile << '[';
        // std::cout << nodetmp << std::endl;
        double tmp = Sgt_Trees[netid].query_sum(1, id[netid][nodetmp]);

        if (tmp != (int)tmp)
          outputFile << std::fixed << std::setprecision(1) << tmp;
        else {
          outputFile << (int)tmp;
        }
        outputFile << ']' << std::endl;
        //-----
      }
    }
    std::ofstream outputFile1("./design.MAX.out");
    outputFile1 << std::fixed << std::setprecision(1) << Glo_Max_Delay;
    outputFile1.close();

    outputFile.close();
  } else {
    // 错误处理
    std::cerr << "无法打开文件" << std::endl;
  }
}
