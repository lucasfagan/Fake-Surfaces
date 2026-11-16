/* fastscan.c -- C kernel for the sigma-assignment sweep in
 * fakesurfaces/enumerate.py.
 *
 * Reads a skeleton spec (precomputed in Python by fakesurfaces/enumerate.py),
 * sweeps all 6^E sheet-gluing assignments with an odometer, walks the
 * corner-chaining orbits, and prints every assignment that yields exactly
 * t+1 orbits (one line per hit: disks as signed edge words, ';' between
 * disks).  Acyclicity and isomorphism dedupe stay in Python.
 *
 * usage: fastscan spec.txt [--fix d0 d1 ...] [--count-only]
 *   --fix pins the digits of the LAST edges (E-1, E-2, ...) for chunking.
 *
 * Spec format (whitespace-separated ints):
 *   t E
 *   selfmem: E lines of 6 ints   msheet slot of side_corners[m][s][k]
 *   oedge:   E lines of 6 ints   (m, arr, cpos) -> other edge
 *   oside:   E lines of 6 ints   (m, arr, cpos) -> other side
 *   omem:    E lines of 6 ints   (m, arr, cpos) -> msheet slot at other
 *   everts:  E lines of 2 ints   skeleton vertex of the x and y ends
 * (each line ordered s=0 k=0..2 then s=1 k=0..2 / arr=0 cpos=0..2 then
 *  arr=1 cpos=0..2)
 *
 * With --acyclic, hits are additionally filtered by the paper's
 * determinant test (collapse a maximal forest, |det| of the non-tree
 * boundary matrix == 1), exactly mirroring FakeSurface.is_acyclic for
 * surfaces living on this connected skeleton.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* must match itertools.permutations(range(3)) order */
static const int P3[6][3] = {
    {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
};

#define MAXE 32

static int E, t, target;
static int sm[MAXE][2][3];      /* selfmem */
static int oe[MAXE][2][3];      /* other edge */
static int os_[MAXE][2][3];     /* other side */
static int om[MAXE][2][3];      /* other msheet slot */
static int ev[MAXE][2];         /* skeleton vertex of each edge end */
static int msheet[6 * MAXE];
static int assign[MAXE];

/* non-tree edges of a maximal forest (fixed per skeleton) */
static int nontree[MAXE], n_nontree, col_of_edge[MAXE];

static void build_forest(void)
{
    int parent[MAXE];
    for (int v = 0; v < t; v++) parent[v] = v;
    n_nontree = 0;
    for (int m = 0; m < E; m++) {
        int a = ev[m][0], b = ev[m][1];
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        while (parent[b] != b) b = parent[b] = parent[parent[b]];
        col_of_edge[m] = -1;
        if (a == b) {
            col_of_edge[m] = n_nontree;
            nontree[n_nontree++] = m;
        } else {
            parent[a] = b;
        }
    }
}

/* |det| == 1 of the (t+1) x (t+1) boundary matrix over the non-tree
 * edges; Bareiss fraction-free elimination, entries stay tiny */
static int det_is_unit(const int *words, const int *dstart,
                       const int *dend, int ndisks)
{
    long long mat[MAXE][MAXE];
    int n = n_nontree;
    if (ndisks != n) return 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) mat[i][j] = 0;
    for (int i = 0; i < ndisks; i++)
        for (int j = dstart[i]; j < dend[i]; j++) {
            int x = words[j];
            int m = (x > 0 ? x : -x) - 1;
            int c = col_of_edge[m];
            if (c >= 0) mat[i][c] += x > 0 ? 1 : -1;
        }
    long long prev = 1;
    for (int col = 0; col < n; col++) {
        int piv = -1;
        for (int r = col; r < n; r++)
            if (mat[r][col]) { piv = r; break; }
        if (piv < 0) return 0;
        if (piv != col) {
            for (int j = 0; j < n; j++) {
                long long tmp = mat[col][j];
                mat[col][j] = mat[piv][j];
                mat[piv][j] = -tmp;
            }
        }
        for (int r = col + 1; r < n; r++)
            for (int j = col + 1; j < n; j++)
                mat[r][j] = (mat[r][j] * mat[col][col] -
                             mat[r][col] * mat[col][j]) / prev;
        prev = mat[col][col];
    }
    long long d = mat[n - 1][n - 1];
    return d == 1 || d == -1;
}

static void set_edge(int m, int pi)
{
    const int *p = P3[pi];
    for (int k = 0; k < 3; k++) {
        msheet[sm[m][0][k]] = k;
        msheet[sm[m][1][p[k]]] = k;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s spec.txt [--fix d ...] [--count-only]\n",
                argv[0]);
        return 2;
    }
    int nfix = 0, fixes[MAXE], count_only = 0, acyclic = 0;
    for (int a = 2; a < argc; a++) {
        if (!strcmp(argv[a], "--count-only")) {
            count_only = 1;
        } else if (!strcmp(argv[a], "--acyclic")) {
            acyclic = 1;
        } else if (!strcmp(argv[a], "--fix")) {
            for (a++; a < argc && argv[a][0] != '-'; a++)
                fixes[nfix++] = atoi(argv[a]);
            a--;
        }
    }

    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
    if (fscanf(f, "%d %d", &t, &E) != 2 || E > MAXE) {
        fprintf(stderr, "bad spec header\n");
        return 2;
    }
    int *blocks[4] = {&sm[0][0][0], &oe[0][0][0], &os_[0][0][0],
                      &om[0][0][0]};
    for (int b = 0; b < 4; b++)
        for (int i = 0; i < E * 6; i++)
            if (fscanf(f, "%d", blocks[b] + i) != 1) {
                fprintf(stderr, "bad spec body\n");
                return 2;
            }
    for (int i = 0; i < E * 2; i++)
        if (fscanf(f, "%d", &ev[0][0] + i) != 1) {
            fprintf(stderr, "bad spec everts\n");
            return 2;
        }
    fclose(f);
    build_forest();

    target = t + 1;
    for (int m = 0; m < E; m++) { assign[m] = 0; set_edge(m, 0); }
    for (int i = 0; i < nfix; i++) {
        int m = E - 1 - i;
        assign[m] = fixes[i];
        set_edge(m, fixes[i]);
    }
    int sweepE = E - nfix;

    char used[3 * MAXE];
    int words[3 * MAXE];
    int dstart[MAXE + 2], dend[MAXE + 2];
    unsigned long long processed = 0, hits = 0;
    char obuf[1 << 16];
    setvbuf(stdout, obuf, _IOFBF, sizeof obuf);

    for (;;) {
        processed++;
        memset(used, 0, 3 * E);
        int ok = 1, ndisks = 0, wlen = 0;
        for (int seed = 0; seed < 3 * E; seed++) {
            if (used[seed]) continue;
            int m0 = seed / 3, k0 = seed % 3;
            int m = m0, k = k0, d = 0;
            int start = wlen;
            for (;;) {
                int pid = 3 * m + k;
                if (used[pid]) { ok = 0; break; }
                used[pid] = 1;
                words[wlen++] = d ? -(m + 1) : (m + 1);
                /* arrival side, corner position there, then hop through
                 * the corner to the other membership's edge/side/sheet */
                int arr = d ? 0 : 1;
                int cpos = arr ? P3[assign[m]][k] : k;
                int m2 = oe[m][arr][cpos];
                int s2 = os_[m][arr][cpos];
                int k2 = msheet[om[m][arr][cpos]];
                m = m2; k = k2; d = s2;
                if (m == m0 && k == k0 && d == 0) break;
            }
            if (!ok) break;
            dstart[ndisks] = start;
            dend[ndisks] = wlen;
            ndisks++;
            if (ndisks > target) { ok = 0; break; }
        }
        if (ok && ndisks == target &&
            (!acyclic || det_is_unit(words, dstart, dend, ndisks))) {
            hits++;
            if (!count_only) {
                for (int i = 0; i < ndisks; i++) {
                    if (i) fputs(";", stdout);
                    for (int j = dstart[i]; j < dend[i]; j++)
                        printf(j > dstart[i] ? " %d" : "%d", words[j]);
                }
                fputs("\n", stdout);
            }
        }

        int pos = 0;
        while (pos < sweepE) {
            assign[pos]++;
            if (assign[pos] < 6) { set_edge(pos, assign[pos]); break; }
            assign[pos] = 0;
            set_edge(pos, 0);
            pos++;
        }
        if (pos == sweepE) break;
    }
    fflush(stdout);
    fprintf(stderr, "processed %llu hits %llu\n", processed, hits);
    return 0;
}
