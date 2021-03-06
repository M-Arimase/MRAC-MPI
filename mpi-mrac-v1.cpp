#include <chrono>
#include <map>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "MRAC.h"
using namespace std;
using namespace chrono;

#define START_FILE_NO 1
#define END_FILE_NO 1

struct FIVE_TUPLE {
  char key[13];
};
typedef vector<FIVE_TUPLE> TRACE;
TRACE traces[END_FILE_NO - START_FILE_NO + 1];

void ReadInTraces(const char *trace_prefix) {
  for (int datafileCnt = START_FILE_NO; datafileCnt <= END_FILE_NO;
       ++datafileCnt) {
    char datafileName[100];
    sprintf(datafileName, "%s%d.dat", trace_prefix, datafileCnt - 1);
    FILE *fin = fopen(datafileName, "rb");

    FIVE_TUPLE tmp_five_tuple;
    traces[datafileCnt - 1].clear();
    while (fread(&tmp_five_tuple, 1, 13, fin) == 13) {
      traces[datafileCnt - 1].push_back(tmp_five_tuple);
    }
    fclose(fin);

    printf("Successfully read in %s, %ld packets\n", datafileName,
           traces[datafileCnt - 1].size());
  }
  printf("\n");
}

void split(vector<int> &nb_split, int max_size) {
  nb_split.resize(max_size + 1);
  vector<vector<int>> dp;

  dp.resize(max_size + 1);
  for (int i = 0; i <= max_size; i++) {
    dp[i].resize(7);
    for (int j = 0; j < 7; j++) {
      dp[i][j] = 0;
    }
  }

  for (int i = 0; i <= max_size; i++) {
    dp[i][1] = 1;
    nb_split[i] = 1;

    int limit = 0;
    if (i > 2000) {
      limit = 3;
    } else if (i > 300) {
      limit = 4;
    } else if (i > 50) {
      limit = 5;
    } else {
      limit = 7;
    }

    for (int j = 2; j < limit; j++) {
      if (i - 1 >= j - 1) {
        dp[i][j] += dp[i - 1][j - 1];
      }
      if (i - j >= j) {
        dp[i][j] += dp[i - j][j];
      }
      nb_split[i] += dp[i][j] * j;
    }
  }
}

double get_p_from_beta(EMFSD::BetaGenerator &bt, double *now_dist, double now_n, int w) {
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

void mrac_worker(int world_rank, int world_size) {
  int *split = new int[world_size];

  MPI_Bcast(split, world_size, MPI_INT, 0, MPI_COMM_WORLD);

  int max_size = split[world_size - 1];
  int w;
  int lend = split[world_rank - 1] + 1;
  int rend = split[world_rank];

  // cout << "processor " << world_rank << " [" << lend << " " << rend << "]"
  //      << endl;

  int *counter_dist = new int[max_size + 1];
  double *dist_new = new double[max_size + 1];
  double *dist_old = new double[max_size + 1];
  double n_new, n_old, n_ns;
  double *ns = new double[max_size + 1];

  MPI_Bcast(counter_dist, max_size + 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dist_new, max_size + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_new, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&w, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int epoch = 1; epoch <= 20; epoch += 1) {
    /*
    cout << "processor " << world_rank << " epoch " << epoch << " begin"
         << endl;
    auto t1 = system_clock::now();
    */

    swap(dist_old, dist_new);
    n_old = n_new;

    double lambda = n_old / double(w);

    memset(ns, 0, (max_size + 1) * sizeof(double));

    for (uint32_t i = lend; i <= rend; ++i) {
      if (counter_dist[i] == 0) {
        continue;
      }

      EMFSD::BetaGenerator bts1(i), bts2(i);
      double sum_p = 0;
      while (bts1.get_next()) {
        double p = get_p_from_beta(bts1, dist_old, n_old, w);
        sum_p += p;
      }
      while (bts2.get_next()) {
        double p = get_p_from_beta(bts2, dist_old, n_old, w);
        p = counter_dist[i] * p / sum_p;
        for (int j = 0; j < bts2.now_flow_num; ++j) {
          ns[bts2.now_result[j]] += p;
        }
      }
    }

    n_ns = 0;
    for (uint32_t i = 1; i <= max_size; i++) {
      n_ns += ns[i];
    }

    /*
    auto t2 = system_clock::now();
    auto duration = duration_cast<microseconds>(t2 - t1).count();
    cout << "processor " << world_rank << " epoch " << epoch << " finish"
         << " duration " << double(duration) / 1000000 << " " << sum << endl;
    */

    MPI_Allreduce(&n_ns, &n_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // cout << "processor " << world_rank << " " << n_ns << " " << n_new <<
    // endl;
    MPI_Allreduce(ns, dist_new, max_size + 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    for (uint32_t i = 1; i <= max_size; i++) {
      dist_new[i] /= n_new;
      // cout << "processor " << world_rank << " " << i << " " << ns[i] << " "
      //      << dist_new[i] << endl;
    }
  }

  delete[] split;
  delete[] counter_dist;
  delete[] dist_new;
  delete[] dist_old;
  delete[] ns;
}

void mrac_controller(int world_size) {
  ReadInTraces("data/");

#define TOT_MEM_IN_BYTES (600 * 1024)
  MRAC<4, TOT_MEM_IN_BYTES> *mrac = NULL;

  for (int datafileCnt = START_FILE_NO; datafileCnt <= END_FILE_NO;
       ++datafileCnt) {
    map<string, int> Real_Freq;
    map<int, int> Real_Dist;
    mrac = new MRAC<4, TOT_MEM_IN_BYTES>();

    int packet_cnt = (int)traces[datafileCnt - 1].size();
    for (int i = 0; i < packet_cnt; ++i) {
      mrac->insert((uint8_t *)(traces[datafileCnt - 1][i].key));

      string str((const char *)(traces[datafileCnt - 1][i].key), 4);
      Real_Freq[str]++;
    }

    for (auto &it : Real_Freq) {
      auto size = it.second;
      Real_Dist[size] += 1;
    }

    mrac->collect_fsd();
    int max_size = mrac->em_fsd_algo->counter_dist.size() - 1;

    vector<int> nb_split;
    split(nb_split, max_size);

    int64_t sum_split = 0;

    for (int i = 0; i <= max_size; i++) {
      if (mrac->em_fsd_algo->counter_dist[i] > 0) {
        sum_split += nb_split[i];
      }
    }

    int64_t split_per_rank = sum_split / (world_size - 1);
    int *split = new int[world_size];

    split[0] = 0;
    for (int i = 1; i < world_size; i++) {
      split[i] = max_size;
    }
    for (int64_t i = 0, j = 0, sum = 0; i <= max_size; i++) {
      if (mrac->em_fsd_algo->counter_dist[i] > 0) {
        sum += nb_split[i];
      }
      if (sum >= split_per_rank) {
        split[++j] = i;
        // cout << "assign " << j << " " << sum << endl;
        sum = 0;
      }
    }

    MPI_Bcast(split, world_size, MPI_INT, 0, MPI_COMM_WORLD);

    int *counter_dist = new int[max_size + 1];
    for (int i = 0; i <= max_size; i++) {
      counter_dist[i] = mrac->em_fsd_algo->counter_dist[i];
    }
    MPI_Bcast(counter_dist, max_size + 1, MPI_INT, 0, MPI_COMM_WORLD);

    double *dist_new = new double[max_size + 1];
    double *dist_old = new double[max_size + 1];
    for (int i = 0; i <= max_size; i++) {
      dist_new[i] = mrac->em_fsd_algo->dist_new[i];
    }
    MPI_Bcast(dist_new, max_size + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    memset(dist_new, 0, (max_size + 1) * sizeof(double));
    memset(dist_old, 0, (max_size + 1) * sizeof(double));

    double n_new = mrac->em_fsd_algo->n_new;
    double n_old = 0;
    MPI_Bcast(&n_new, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int w = MRAC<4, TOT_MEM_IN_BYTES>::w;
    MPI_Bcast(&w, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int epoch = 1; epoch <= 20; epoch += 1) {
      auto t1 = system_clock::now();
      MPI_Allreduce(&n_old, &n_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(dist_old, dist_new, max_size + 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      auto t2 = system_clock::now();

      auto duration = duration_cast<microseconds>(t2 - t1).count();
      double WMRE_1 = 0;
      double WMRE_2 = 0;
      for (int i = 0; i <= max_size; ++i) {
        WMRE_1 += fabs(dist_new[i] - Real_Dist[i]);
        WMRE_2 += (dist_new[i] + Real_Dist[i]) / 2;
      }
      cout << epoch << " " << WMRE_1 / WMRE_2 << " "
           << double(duration) / 1000000 << endl;
    }

    delete mrac;
    delete[] split;
    delete[] counter_dist;
    delete[] dist_new;

    Real_Freq.clear();
    Real_Dist.clear();
    nb_split.clear();
  }
}

int main() {
  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    mrac_controller(world_size);
  } else {
    mrac_worker(world_rank, world_size);
  }

  MPI_Finalize();
}
