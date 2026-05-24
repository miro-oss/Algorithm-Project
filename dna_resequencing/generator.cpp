#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <array>
#include <tuple>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <climits>

using namespace std;
/*
 *  시뮬레이션 시나리오
 *  ──────────────────
 *  1) Reference genome R (길이 N) 랜덤 생성
 *  2) R에 NUM_SNPS개 SNP를 도입 → Sample genome S (= 정답)
 *  3) S의 랜덤 위치에서 L짜리 read를 M개 추출
 *     → 각 염기에 SEQ_ERROR_RATE 확률로 sequencing error 삽입
 *
 *  read를 Reference에 정렬하면 관찰되는 mismatch =
 *     (read가 걸치는 SNP 수) + (sequencing error 수)
 *
 *  N=10000, 80 SNPs → SNP 밀도 0.8 %
 *  L=50 → read 당 평균 SNP ≈ 0.4개
 *  SEQ_ERROR_RATE 1 % → read 당 평균 error ≈ 0.5개
 *  ⇒ 평균 total mismatch ≈ 0.9 → D=2이면 약 92 %의 read가 매핑 가능
 *  Coverage = 5000 × 50 / 10000 = 25×
 */

constexpr int    N              = 10000;
constexpr int    M              = 5000;
constexpr int    L              = 50;
constexpr int    NUM_SNPS       = 80;
constexpr double SEQ_ERROR_RATE = 0.01;

int main() {
    mt19937 rng(42);
    const char B[] = "ACGT";

    /* ── 1. Reference genome ── */
    string ref(N, 'A');
    for (int i = 0; i < N; i++)
        ref[i] = B[rng() % 4];

    /* ── 2. Sample genome (SNP 도입) ── */
    string sample = ref;
    set<int> used;
    uniform_int_distribution<int> pdist(0, N - 1);

    while ((int)used.size() < NUM_SNPS) {
        int p = pdist(rng);
        if (used.count(p)) continue;
        used.insert(p);
        char orig = sample[p];
        char mut;
        do { mut = B[rng() % 4]; } while (mut == orig);
        sample[p] = mut;
    }

    /* ── 3. Reads 생성 ── */
    uniform_int_distribution<int> rdist(0, N - L);
    uniform_real_distribution<double> prob(0.0, 1.0);

    vector<string> reads(M);
    for (int i = 0; i < M; i++) {
        int st = rdist(rng);
        reads[i] = sample.substr(st, L);

        for (int j = 0; j < L; j++) {
            if (prob(rng) < SEQ_ERROR_RATE) {
                char orig = reads[i][j];
                char err;
                do { err = B[rng() % 4]; } while (err == orig);
                reads[i][j] = err;
            }
        }
    }

    /* ── 4. 파일 출력 ── */
    system("mkdir -p data");

    ofstream(  "data/reference.txt") << ref    << "\n";
    ofstream f("data/reads.txt");
    f << M << "\n";
    for (auto& r : reads) f << r << "\n";
    ofstream(  "data/answer.txt")    << sample << "\n";

    /* ── 5. 통계 출력 ── */
    int snps = 0;
    for (int i = 0; i < N; i++)
        if (ref[i] != sample[i]) snps++;

    cout << "=== Data Generation Complete ===\n"
         << "Reference:  " << N << " bp\n"
         << "SNPs:       " << snps << "\n"
         << "Reads:      " << M << " × " << L << " bp\n"
         << "Coverage:   " << (double)M * L / N << "×\n"
         << "Error rate: " << SEQ_ERROR_RATE * 100 << " %\n"
         << "Output:     data/{reference,reads,answer}.txt\n";
}
