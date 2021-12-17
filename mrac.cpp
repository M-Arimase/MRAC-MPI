#include <chrono>
#include <map>
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

int main() {
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

    int max_size = 0;
    for (auto [flow, size] : Real_Freq) {
      Real_Dist[size] += 1;
      max_size = max(max_size, size);
    }

    mrac->collect_fsd();

    printf("%d.dat: ", datafileCnt - 1);
    printf("flow size distribution: <flow size, count>\n");

    for (int epoch = 1; epoch <= 20; epoch += 1) {
      auto t1 = system_clock::now();
      mrac->next_epoch();
      auto t2 = system_clock::now();
      auto duration = duration_cast<microseconds>(t2 - t1).count();

      vector<double> dist;
      mrac->get_distribution(dist);

      double WMRE_1 = 0;
      double WMRE_2 = 0;
      for (int i = 0, j = 0; i <= max_size; ++i) {
        WMRE_1 += fabs(dist[i] - Real_Dist[i]);
        WMRE_2 += (dist[i] + Real_Dist[i]) / 2;
      }
      cout << epoch << " " << WMRE_1 / WMRE_2 << " "
           << double(duration) / 1000000 << endl;
    }

    delete mrac;
    Real_Freq.clear();
    Real_Dist.clear();
  }
}
