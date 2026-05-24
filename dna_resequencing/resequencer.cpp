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
using u64 = uint64_t;

/* ═══════════════════════════════════════════
 *  튜닝 파라미터
 *  ─────────────
 *  K         : k-mer 길이 (seed 정밀도 ↔ 민감도 trade-off)
 *  W         : minimizer window (작을수록 seed 밀도 ↑, index 크기 ↑)
 *  D         : 최대 허용 edit distance
 *  TOP_K     : chain 후 검증할 상위 후보 수
 *  MIN_CHAIN : chain에 필요한 최소 seed 수
 *
 *  Pigeonhole: D=2 → read를 3등분(≥16bp)하면 최소 1구간 error-free
 *  → 16-K+1 ≥ W 일 때 seed 보장. K=15,W=5 → 2 < 5 이므로
 *    이론적 보장은 없지만, 25× coverage에서는 충분히 동작함.
 *    민감도를 높이려면 K=11,W=5 또는 MIN_CHAIN=1로 조정.
 * ═══════════════════════════════════════════ */

constexpr int K         = 11;
constexpr int W         = 5;
constexpr int D         = 2;
constexpr int TOP_K     = 50;
constexpr int MIN_CHAIN = 1;
constexpr int SEARCH_RADIUS = D;

/* ──── 인코딩 / 디코딩 ──── */
int ENC[256];

void init_encoding() {
    ENC['A'] = 0;  ENC['C'] = 1;
    ENC['G'] = 2;  ENC['T'] = 3;
}

inline int  enc(char c) { return ENC[(unsigned char)c]; }
inline char dec(int  v) { return "ACGT"[v]; }


/* ═══════════════════════════════════════════
 *  1) Minimizer Index
 * ═══════════════════════════════════════════ */

// DNA k-mer → 2-bit 정수 해싱 (k ≤ 32이면 충돌 없는 perfect hash)
inline u64 hash_kmer(const string& s, int pos, int k) {
    u64 h = 0;
    for (int i = 0; i < k; i++)
        h = (h << 2) | enc(s[pos + i]);
    return h;
}

using Mini = pair<u64, int>;   // (hash, position)

// 서열에서 minimizer 목록 추출
// k-mer hash를 미리 계산한 뒤 sliding window로 최솟값 선택
vector<Mini> extract_minimizers(const string& seq, int k, int w) {
    int n = (int)seq.size();
    if (n < k + w - 1) return {};

    int nk = n - k + 1;
    vector<u64> H(nk);
    for (int i = 0; i < nk; i++)
        H[i] = hash_kmer(seq, i, k);

    vector<Mini> res;
    int last = -1;

    for (int i = 0; i <= n - k - w + 1; i++) {
        u64 mh = UINT64_MAX;
        int  mp = -1;
        for (int j = i; j < i + w; j++) {
            if (H[j] < mh) { mh = H[j]; mp = j; }
        }
        if (mp != last) {
            res.emplace_back(mh, mp);
            last = mp;
        }
    }
    return res;
}

// Reference minimizer index 구축
// key = minimizer hash, value = reference 내 등장 위치 리스트
unordered_map<u64, vector<int>> build_index(const string& ref) {
    unordered_map<u64, vector<int>> idx;
    for (auto& [h, p] : extract_minimizers(ref, K, W))
        idx[h].push_back(p);
    return idx;
}


/* ═══════════════════════════════════════════
 *  2) Seed 수집
 * ═══════════════════════════════════════════ */

using Hit = pair<int, int>;   // (read_pos, ref_pos)

vector<Hit> collect_seeds(const string& read,
                          const unordered_map<u64, vector<int>>& idx) {
    vector<Hit> hits;
    for (auto& [h, rp] : extract_minimizers(read, K, W)) {
        auto it = idx.find(h);
        if (it != idx.end())
            for (int rfp : it->second)
                hits.emplace_back(rp, rfp);
    }
    return hits;
}


/* ═══════════════════════════════════════════
 *  3) Seed Chaining
 *  ─────────────────
 *  read가 ref position p에 매핑 → seed (r, p+r)
 *  ⇒ diagonal = ref_pos - read_pos = p (일정)
 *  같은 diagonal에 seed가 MIN_CHAIN개 이상 모이면 진짜 후보
 * ═══════════════════════════════════════════ */

vector<int> chain_seeds(const vector<Hit>& hits) {
    if (hits.empty()) return {};

    unordered_map<int, int> dcnt;              // diagonal → seed 개수
    for (auto& [rp, rfp] : hits)
        dcnt[rfp - rp]++;

    vector<pair<int,int>> cands;               // (seed_count, ref_start)
    for (auto& [d, c] : dcnt)
        if (c >= MIN_CHAIN)
            cands.emplace_back(c, d);

    sort(cands.begin(), cands.end(), greater<>());

    vector<int> res;
    int lim = min((int)cands.size(), TOP_K);
    for (int i = 0; i < lim; i++)
        res.push_back(cands[i].second);
    return res;
}


/* ═══════════════════════════════════════════
 *  4) Myers Bit-Parallel 검증
 *  ──────────────────────────
 *  L=50 ≤ 64이므로 uint64_t 1워드로 처리 (word-level parallelism)
 *  
 *  Pv, Mv : positive / negative delta 비트벡터
 *  Peq[c] : read에서 문자 c가 등장하는 위치 bit mask
 *  score  : 현재까지의 edit distance
 *
 *  각 chain 후보 rs 주변에서 가능한 시작점 st를 검사한다.
 *  매 시작점마다 Myers 상태를 새로 초기화하여 read와
 *  reference[st : st + L]의 global edit distance를 계산한다.
 * ═══════════════════════════════════════════ */

struct Match { int pos, score; };

int myers_distance_at(const string& ref, int st, int L,
                      const array<u64, 4>& Peq, u64 top, u64 mask) {
    u64 Pv = mask, Mv = 0;
    int sc = L;

    for (int j = 0; j < L; j++) {
        u64 Eq = Peq[enc(ref[st + j])];

        u64 Xv = Eq | Mv;
        u64 Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        u64 Ph = Mv | ~(Xh | Pv);
        u64 Mh = Pv & Xh;

        Ph &= mask;
        Mh &= mask;

        if      (Ph & top) sc++;
        else if (Mh & top) sc--;

        Ph = ((Ph << 1) | 1ULL) & mask;
        Mh = (Mh << 1)          & mask;
        Pv = (Mh | ~(Xv | Ph))  & mask;
        Mv = (Ph & Xv)          & mask;
    }

    return sc;
}

vector<Match> myers_verify(const string& ref, const string& read,
                           const vector<int>& cands, int d) {
    int L    = (int)read.size();
    int rlen = (int)ref.size();
    if (L == 0 || L > 64 || rlen < L) return {};

    u64 top  = 1ULL << (L - 1);
    u64 mask = (L == 64) ? ~u64{0} : (1ULL << L) - 1;

    // Peq 전처리 — read의 각 문자 위치를 bit mask로 기록
    array<u64, 4> Peq{};
    for (int i = 0; i < L; i++)
        Peq[enc(read[i])] |= 1ULL << i;

    vector<Match> res;

    for (int rs : cands) {
        int lo = max(0, rs - SEARCH_RADIUS);
        int hi = min(rlen - L, rs + SEARCH_RADIUS);
        for (int st = lo; st <= hi; st++) {
            int sc = myers_distance_at(ref, st, L, Peq, top, mask);
            if (sc <= d) {
                res.push_back({st, sc});
            }
        }
    }
    return res;
}


/* ═══════════════════════════════════════════
 *  5) Consensus 복원
 *  ─────────────────
 *  각 position에 대해 매핑된 read들이 투표
 *  가중치 = 1/(edit_distance + 1)  →  정확한 매핑일수록 높은 영향력
 * ═══════════════════════════════════════════ */

string consensus(const string& ref,
                 const vector<tuple<int,int,int>>& maps,
                 const vector<string>& reads) {
    int n = (int)ref.size();
    vector<array<double,4>> votes(n, {0.0, 0.0, 0.0, 0.0});

    for (auto& [ri, rs, sc] : maps) {
        double w = 1.0 / (sc + 1);
        const string& rd = reads[ri];
        for (int i = 0; i < (int)rd.size(); i++) {
            int p = rs + i;
            if (p >= 0 && p < n)
                votes[p][enc(rd[i])] += w;
        }
    }

    string res(n, 'N');
    for (int i = 0; i < n; i++) {
        int best = 0;
        for (int b = 1; b < 4; b++)
            if (votes[i][b] > votes[i][best]) best = b;
        res[i] = (votes[i][best] > 0.0) ? dec(best) : ref[i];
    }
    return res;
}


/* ═══════════════════════════════════════════
 *  6) Main Pipeline
 * ═══════════════════════════════════════════ */

string resequence(const string& ref, const vector<string>& reads) {
    auto idx = build_index(ref);
    vector<tuple<int,int,int>> all_maps;     // (read_idx, ref_start, score)
    int mapped = 0;

    for (int i = 0; i < (int)reads.size(); i++) {
        auto hits  = collect_seeds(reads[i], idx);
        auto cands = chain_seeds(hits);
        if (cands.empty()) continue;

        auto matches = myers_verify(ref, reads[i], cands, D);
        if (!matches.empty()) {
            auto best = *min_element(matches.begin(), matches.end(),
                [](auto& a, auto& b) { return a.score < b.score; });
            all_maps.emplace_back(i, best.pos, best.score);
            mapped++;
        }
    }

    cerr << "Mapped: " << mapped << " / " << reads.size() << " reads\n";
    return consensus(ref, all_maps, reads);
}

int main() {
    init_encoding();
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    /* ── 입력 ── */
    ifstream fref("data/reference.txt");
    string ref;
    getline(fref, ref);
    fref.close();

    ifstream fr("data/reads.txt");
    int M;
    fr >> M;
    vector<string> reads(M);
    for (int i = 0; i < M; i++) fr >> reads[i];
    fr.close();

    cerr << "Reference: " << ref.size() << " bp, Reads: " << M << "\n";

    /* ── 실행 ── */
    auto t0 = chrono::steady_clock::now();
    string result = resequence(ref, reads);
    auto t1 = chrono::steady_clock::now();

    double sec = chrono::duration<double>(t1 - t0).count();
    cerr << "Elapsed: " << fixed << setprecision(3) << sec << " s\n";

    /* ── 결과 저장 ── */
    ofstream("data/result.txt") << result << "\n";

    /* ── 정답 비교 (answer.txt 존재 시) ── */
    ifstream fa("data/answer.txt");
    if (fa.is_open()) {
        string ans;
        getline(fa, ans);
        fa.close();

        int diff_ref = 0, diff_res = 0;
        for (int i = 0; i < (int)ref.size(); i++) {
            if (ref[i]    != ans[i]) diff_ref++;
            if (result[i] != ans[i]) diff_res++;
        }

        cerr << "\n=== Evaluation ===\n"
             << "Ref    ↔ Answer : " << diff_ref << " diffs (original SNPs)\n"
             << "Result ↔ Answer : " << diff_res << " diffs (reconstruction error)\n"
             << "Accuracy        : " << fixed << setprecision(2)
             << 100.0 * ((int)ref.size() - diff_res) / (int)ref.size() << " %\n"
             << "SNP recovery    : " << fixed << setprecision(2)
             << 100.0 * (diff_ref - diff_res) / max(diff_ref, 1) << " %\n";
    }
}
