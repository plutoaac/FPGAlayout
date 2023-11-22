#include <bits/stdc++.h>

using namespace std;

#define pii pair<int, int>
int n, m, k;
int wire1;
int wire2;
int wire3;
int wire4;
int a[5];
vector<string> ss;
std::map<pii, pair<int, int>> mpp1;
int ma = 0;
std::map<int, std::map<pii, int>> mpp;
int x, y;
std::map<int, int> mp;
set<pair<int, int>> mp1[10000010];
map<int, vector<pair<int, int>>> vv;
set<int> v1;
// Die 到fpga的映射 每一个Die属于哪个FPGA
std::map<int, int> DieToFpga_Map;
int main() {
  ifstream input1("./design.tdm.out");
  if (input1.is_open()) {
    string line;
    while (getline(input1, line)) {
      if (line[1] == 'D') {
        int i = 0;
        string tmpp;
        for (i = 4; i < line.size(); i++) {
          if (line[i] == ',') {
            break;
          }

          tmpp += line[i];
        }
        x = atoi(tmpp.c_str());
        i += 4;
        tmpp.clear();
        for (; i < line.size(); i++) {
          if (line[i] == ']') {
            break;
          }

          tmpp += line[i];
        }
        y = atoi(tmpp.c_str());
        // cout<<"xy"<<x<<" "<<y<<endl;
        continue;
      }
      string ss;
      string t;
      int i;
      int fff = 0;
      for (i = 1; i < line.size(); i++) {
        if (line[i] == 'g' || line[i] == 'f') {
          if (line[i] == 'f') fff = 1;
          break;
        }
      }
      i++;
      for (; i < line.size(); i++) {
        t += line[i];
      }

      int val = atoi(t.c_str());
      ma = max(ma, val);
      // cout<<val<<endl;
      t.clear();

      int netid;
      for (int i = 1; i < line.size(); i++) {
        if (line[i] == ',' || line[i] == ']') {
          netid = atoi(ss.c_str());
          // cout<<netid <<x<<" "<<y<<" "<<val<<endl;
          if (!fff) {
            mpp[netid][{x, y}] = val;
            mp[netid]++;
            vv[netid].push_back({x, y});
          } else {
            mpp[netid][{y, x}] = val;
            mp[netid]++;
            vv[netid].push_back({y, x});
          }
          ss.clear();
          continue;
        }
        ss += line[i];
      }
      ss.clear();
    }
    input1.close();
  } else {
    cout << "failed" << endl;
  }
  ifstream inputFile("design.route.out");
  if (inputFile.is_open()) {
    string line;
    int netid;
    while (std::getline(inputFile, line)) {
      int aa = 0;
      int bb = 0;
      for (auto &ch : line) {
        if (ch == '[')
          aa++;
        else if (ch == ']')
          bb++;
      }
      if (aa == 1 && bb == 1) {
        string sss;
        for (int i = 1; i < line.size() - 1; i++) {
          sss += line[i];
        }
        netid = atoi(sss.c_str());

        sss.clear();

        a[1] += wire1;
        a[2] += wire2;
        a[3] += wire3;
        a[4] += wire4;
        //   cout << wire1 << " " << wire2 << " " << wire3 << " " << wire4 <<
        //   endl;
        string tt;
        if (wire1 != 0 || wire2 != 0 || wire3 != 0 || wire4 != 0) {
          for (int i = 1; i < line.size() - 1; i++) {
            tt += line[i];
          }
          ss.push_back(tt);
        }

        wire1 = 0;
        wire2 = 0;
        wire3 = 0;
        wire4 = 0;
        // cout<<tt<<endl;
        continue;
      }
      std::vector<int> ans;
      string tmp;
      int tmpidx;
      for (int i = 1; i < line.size(); i++) {
        if (line[i] == ']') {
          ans.push_back(atoi(tmp.c_str()));
          tmp.clear();
          tmpidx = i;
          break;
        }
        if (line[i] == ',') {
          // cout<<atoi(tmp.c_str())<<endl;
          ans.push_back(atoi(tmp.c_str()));
          tmp.clear();
          continue;
        }
        if ('0' <= line[i] && line[i] <= '9') {
          tmp += line[i];
        }
      }
      double Delay = 0.0;

      int flag = 0;
      // cout << ans.size() << endl;
      double ma = 0;
      if (ans.size() >= 1) {
        for (int i = 1; i < ans.size(); i++) {
          if (mpp[netid][{ans[i - 1], ans[i]}] == 0)
            Delay += 1.0;
          else {
            flag++, mp1[netid].insert({ans[i - 1], ans[i]}),
                Delay += 1.0 * (mpp[netid][{ans[i - 1], ans[i]}]);
          }
        }
      }
      for (int i = tmpidx + 2; i < line.size() - 1; i++) {
        // cout << line[i] << " ";
        tmp += line[i];
      }
      // cout << endl;
      //  cout << tmp << endl;
      double DD = atof(tmp.c_str());
      tmp.clear();
      if (flag) Delay += flag * 0.5;
      // cout << Delay << " " << DD << endl;
      // cout << mp[netid] << " " << mp1[netid].size() << endl;
      v1.insert(netid);
      if (Delay != DD) {
        cout << Delay << " " << DD << endl;
      }
      assert(Delay == DD);
    }
    cout << "okk" << endl;
  } else {
    cout << "failed" << endl;
  }
  for (auto netid : v1) {
    if (mp[netid] != mp1[netid].size()) {
      cout << netid << endl;
      cout << mp[netid] << " " << mp1[netid].size() << endl;
      for (auto it : vv[netid]) {
        cout << it.first << " " << it.second << endl;
      }
    }
    assert(mp[netid] == mp1[netid].size());
  }

  std::cout << ma << std::endl;
  cout << endl;
}
