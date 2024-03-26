#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
class ReadDataSource {
 public:
  static ReadDataSource &getInstance() { return readDataSource; }
  ReadDataSource(int x) {
    Max_SLL = x;
    _read_die_network();
    _read_die_position();
    _readnet();
    _read_fpga();
    Max_SLL = get_MaxSLL();
    soo();
  }
  // Wire
  std::vector<std::pair<int, int>> Not_Fpga_Die;
  // 最大的Die间SLL
  int Max_SLL;
  // Die间的连接信息
  std::vector<std::vector<int>> data;
  int DieNum{0};
  // Die 到fpga的映射 每一个Die属于哪个FPGA
  std::unordered_map<int, int> DieToFpga_Map;
  // 结点到Die的映射 每一个结点属于哪个Die
  std::unordered_map<int, int> NodeToDie_Map;
  // 驱动节点的数量
  int DriverNum = 0;
  // 求解的每一组net 驱动到负载结点
  std::map<int, std::vector<std::pair<int, double>>> Task;

  std::vector<std::pair<int, std::vector<std::pair<int, double>>>> Tmp;

  // 读入design.net文件
  bool _readnet();
  // 读入design.die.position文件
  bool _read_die_position();
  // 读入design.dir.network 文件
  bool _read_die_network();
  // 读入design.fpga.die文件
  bool _read_fpga();
  // 求最大的DIE间SLL
  int get_MaxSLL();

  void soo();
  ~ReadDataSource() = default;
  static ReadDataSource readDataSource;
  ReadDataSource(const ReadDataSource &) = delete;
  ReadDataSource &operator=(const ReadDataSource &) = delete;
};
