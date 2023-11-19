#include <bits/stdc++.h>
using namespace std;

int n, m, k;
int wire1;
int wire2;
int wire3;
int wire4;
int a[5];
vector<string> ss;
int main() {
  ifstream inputFile("design.route.out");
  if (inputFile.is_open()) {
    string line;
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
        a[1] += wire1;
        a[2] += wire2;
        a[3] += wire3;
        a[4] += wire4;

        //   cout << wire1 << " " << wire2 << " " << wire3 << " " << wire4 <<
        //   endl;

        string tt;
        if(wire1!=0||wire2!=0||wire3!=0||wire4!=0){
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
          tmpidx = i;
          break;
        }
        if ('0' <= line[i] && line[i] <= '9') {
          tmp += line[i];
          ans.push_back(atoi(tmp.c_str()));
          tmp.clear();
        }
      }
      double Delay = 0.0;
      for (int i = 1; i < ans.size(); i++) {
        if (ans[i] == 0 && ans[i - 1] == 4)
          Delay += 4.5, wire1 = 1;
        else if (ans[i] == 4 && ans[i - 1] == 0)
          Delay += 4.5, wire2 = 1;
        else if (ans[i] == 7 && ans[i - 1] == 3)
          Delay += 4.5, wire3 = 1;
        else if (ans[i] == 3 && ans[i - 1] == 7)
          Delay += 4.5, wire4 = 1;
        else
          Delay += 1.0;
      }

      for (int i = tmpidx + 2; i < line.size() - 1; i++) {
        // cout << line[i] << " ";
        tmp += line[i];
      }
      // cout << endl;
      //  cout << tmp << endl;
      double DD = atof(tmp.c_str());
      tmp.clear();
      //  cout << Delay << " " << DD << endl;
      //  assert(Delay == DD);
    }

  } else {
    cout << "failed" << endl;
  }
  for (int i = 1; i <= 4; i++) {
    cout << a[i] << endl;
  }
  for (auto it : ss) {
    cout << it << " ";
  }
  


  cout<<endl;
}

/*
5
1.2 2.3  1.3 4.5 2.3
4.2 5.3  0.3


*/