#ifndef STREAMMEASUREMENTSYSTEM_EMFSD_H
#define STREAMMEASUREMENTSYSTEM_EMFSD_H

#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using std::vector;

class EMFSD {
private:
  vector<vector<int>> *bts_result;
  vector<double> *p1;
  double *ps;
  bool has_buf;

public:
  uint32_t w;
  vector<int> counter_dist;
  vector<double> dist_old, dist_new;

  vector<double> ns;
  double n_sum;
  double card_init;
  bool inited = false;

  double n_old, n_new;

  struct BetaGenerator {
    int sum;
    int now_flow_num;
    int flow_num_limit;
    vector<int> now_result;

    explicit BetaGenerator(int _sum) : sum(_sum) {
      now_flow_num = 0;
      if (sum > 2000)
        flow_num_limit = 3;
      else if (sum > 300)
        flow_num_limit = 4;
      else if (sum > 50)
        flow_num_limit = 5;
      else
        flow_num_limit = 7;
    }

    bool get_new_comb() {
      for (int j = now_flow_num - 2; j >= 0; --j) {
        int t = ++now_result[j];
        for (int k = j + 1; k < now_flow_num - 1; ++k) {
          now_result[k] = t;
        }
        int partial_sum = 0;
        for (int k = 0; k < now_flow_num - 1; ++k) {
          partial_sum += now_result[k];
        }
        int remain = sum - partial_sum;
        if (remain >= now_result[now_flow_num - 2]) {
          now_result[now_flow_num - 1] = remain;
          return true;
        }
      }

      return false;
    }

    bool get_next() {
      while (now_flow_num < flow_num_limit) {
        switch (now_flow_num) {
        case 0:
          now_flow_num = 1;
          now_result.resize(1);
          now_result[0] = sum;
          return true;
        case 1:
          now_flow_num = 2;
          now_result[0] = 0;
          // fallthrough
        default:
          now_result.resize(now_flow_num);
          if (get_new_comb()) {
            //                        for (int t = 0; t < now_flow_num; ++t) {
            //                            cout << now_result[t] << " ";
            //                        }
            //                        cout << endl;
            return true;
          } else {
            now_flow_num++;
            for (int i = 0; i < now_flow_num - 2; ++i) {
              now_result[i] = 1;
            }
            now_result[now_flow_num - 2] = 0;
          }
        }
      }

      return false;
    }
  };

  static int factorial(int n) {
    if (n == 0 || n == 1)
      return 1;
    return factorial(n - 1) * n;
  }

  double get_p_from_beta(BetaGenerator &bt, double lambda,
                         vector<double> &now_dist, double now_n) {
    std::map<uint32_t, uint32_t> mp;
    for (int i = 0; i < bt.now_flow_num; ++i) {
      mp[bt.now_result[i]]++;
    }

    // double ret = std::exp(-lambda);
    double ret = 1.0;
    double div_n_w = now_n / w;
    for (auto &kv : mp) {
      uint32_t fi = kv.second;
      uint32_t si = kv.first;
      double lambda_i = now_dist[si] * div_n_w;
      for (uint32_t k = 1; k <= fi; ++k)
        ret *= lambda_i / k;
    }

    return ret;
  }

  int collect_counters(uint32_t *counters) {
    // collect counter information as the dist init
    uint32_t max_counter_val = 0;
    for (uint32_t i = 0; i < w; ++i) {
      max_counter_val = std::max(max_counter_val, counters[i]);
    }
    counter_dist.resize(max_counter_val + 1);
    std::fill(counter_dist.begin(), counter_dist.end(), 0);
    for (uint32_t i = 0; i < w; ++i) {
      counter_dist[counters[i]]++;
    }
    /*
            std::cout << w << endl;
            std::cout << counter_dist[0] << endl;
            std::cout << max_counter_val << endl;
    */
    return max_counter_val;
  }

  int collect_counters(uint16_t *counters) {
    // collect counter information as the dist init
    uint16_t max_counter_val = 0;
    for (uint32_t i = 0; i < w; ++i) {
      max_counter_val = std::max(max_counter_val, counters[i]);
    }
    counter_dist.resize(max_counter_val + 1);
    std::fill(counter_dist.begin(), counter_dist.end(), 0);
    for (uint32_t i = 0; i < w; ++i) {
      counter_dist[counters[i]]++;
    }
    std::cout << w << endl;
    std::cout << counter_dist[0] << endl;
    std::cout << max_counter_val << endl;

    return max_counter_val;
  }

public:
  EMFSD() {
    //        BetaGenerator a(15);
    //        while (a.get_next()) {
    //            for (int j = 0; j < a.now_flow_num; ++j) {
    //                cout << a.now_result[j] << " ";
    //            }
    //            cout << endl;
    //        }
    has_buf = false;
  }

  void set_counters(uint32_t _w, uint32_t *counters) {
    inited = true;
    w = _w;
    int max_cnt = collect_counters(counters);
    n_new = w - counter_dist[0];
    //        std::cout << "w: " << w << std::endl;
    //        std::cout << "dist0: " << counter_dist[0] << std::endl;
    dist_new.resize(counter_dist.size());
    ns.resize(counter_dist.size());
    for (uint32_t i = 1; i < counter_dist.size(); ++i) {
      dist_new[i] = counter_dist[i] / double(w - counter_dist[0]);
      ns[i] = counter_dist[i];
    }
  }

  void set_counters(uint32_t _w, uint16_t *counters) {
    inited = true;
    w = _w;
    int max_cnt = collect_counters(counters);
    n_new = w - counter_dist[0];
    dist_new.resize(counter_dist.size());
    ns.resize(counter_dist.size());
    for (uint32_t i = 1; i < counter_dist.size(); ++i) {
      dist_new[i] = counter_dist[i] / double(w - counter_dist[0]);
      ns[i] = counter_dist[i];
    }
    card_init = (w * std::log(w / double(counter_dist[0])));
  }

  void next_epoch() {
    dist_old = dist_new;
    n_old = n_new;

    double lambda = n_old / double(w);

    std::fill(ns.begin(), ns.end(), 0);

    for (uint32_t i = 1; i < counter_dist.size(); ++i) {
      // enum how to form val:i
      if (counter_dist[i] == 0) {
        continue;
      }
      BetaGenerator bts1(i), bts2(i);
      double sum_p = 0;
      while (bts1.get_next()) {
        double p = get_p_from_beta(bts1, lambda, dist_old, n_old);
        sum_p += p;
      }
      while (bts2.get_next()) {
        double p = get_p_from_beta(bts2, lambda, dist_old, n_old);
        p = counter_dist[i] * p / sum_p;
        for (int j = 0; j < bts2.now_flow_num; ++j) {
          ns[bts2.now_result[j]] += p;
        }
      }
    }

    n_new = 0;
    for (uint32_t i = 1; i < counter_dist.size(); i++) {
      //            ns[i] = int(ns[i]); // TODO
      n_new += ns[i];
    }
    for (uint32_t i = 1; i < counter_dist.size(); i++) {
      //            ns[i] = int(ns[i]);
      dist_new[i] = ns[i] / n_new;
    }

    n_sum = n_new;
    //        std::cout << ns[1] << std::endl;
  }

  static double get_p_from_beta_1(const vector<int> &bt) {
    std::map<uint32_t, uint32_t> mp;
    for (auto i : bt) {
      mp[i]++;
    }

    double ret = 1.0;
    for (auto &kv : mp) {
      uint32_t fi = kv.second;
      uint32_t si = kv.first;
      ret *= EMFSD::factorial(fi);
    }

    return 1.0 / ret;
  }

  static double get_p_from_beta_2(const vector<int> &bt, const vector<double> &now_dist, double now_n, int w) {
    double ret = 1.0;
    double r = now_n / w;
    for (auto i : bt) {
      ret *= now_dist[i];
      ret *= r;
    }
    return ret;
  }

  void alloc_buf() {
    size_t max_size = counter_dist.size();
    bts_result = new vector<vector<int>>[max_size + 1];
    p1 = new vector<double>[max_size + 1];
    size_t max_bts_size = 0, total_bts_size = 0;
    for (size_t i = 1; i <= counter_dist.size(); i++) {
      if (counter_dist[i] == 0) {
        continue;
      }
      EMFSD::BetaGenerator bts(i);
      while (bts.get_next()) {
        vector<int> tmp(bts.now_result.begin(), bts.now_result.begin() + bts.now_flow_num);
        p1[i].push_back(get_p_from_beta_1(tmp));
        bts_result[i].push_back(move(tmp));
      }
      max_bts_size = max(max_bts_size, bts_result[i].size());
    }
    ps = new double[max_bts_size];
    has_buf = true;
  }

  void next_epoch_save() {
    dist_old = dist_new;
    n_old = n_new;

    double lambda = n_old / double(w);

    std::fill(ns.begin(), ns.end(), 0);

    size_t max_size = counter_dist.size();
    if (!has_buf) {
      std::cerr << "no buffer!!" << std::endl;
      alloc_buf();
    }
    for (uint32_t i = 1; i < counter_dist.size(); ++i) {
      // enum how to form val:i
      if (counter_dist[i] == 0) {
        continue;
      }
      double sum_p = 0;
      for (int j = 0; j < bts_result[i].size(); j++) {
        double p = get_p_from_beta_2(bts_result[i][j], dist_old, n_old, w);
        p *= p1[i][j];
        sum_p += p;
        ps[j] = p;
      }
      for (int j = 0; j < bts_result[i].size(); j++) {
        double p = counter_dist[i] * ps[j] / sum_p;
        for (auto k : bts_result[i][j]) {
          ns[k] += p;
        }
      }
    }

    n_new = 0;
    for (uint32_t i = 1; i < counter_dist.size(); i++) {
      //            ns[i] = int(ns[i]); // TODO
      n_new += ns[i];
    }
    for (uint32_t i = 1; i < counter_dist.size(); i++) {
      //            ns[i] = int(ns[i]);
      dist_new[i] = ns[i] / n_new;
    }

    n_sum = n_new;
    //        std::cout << ns[1] << std::endl;
  }
};

#endif // STREAMMEASUREMENTSYSTEM_EMFSD_H
