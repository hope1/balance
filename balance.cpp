#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cctype>
using namespace std;

namespace Cirno {
    int sgn(int x) {
        return x == 0 ? 0 : (x < 0 ? -1 : 1);
    }

    int abs(int x) {
        return x < 0 ? -x : x;
    }

    int gcd(int a, int b) {
        return b == 0 ? a : gcd(b, a % b);
    }

    struct mat {
        int row, col, acol;
        vector<vector<int> > a;
    };

    void print_mat(const mat& A) {
        for(int i = 0; i < A.row; ++i)
            for(int j = 0; j < A.acol; ++j)
                printf("%d%s", A.a[i][j], j == A.col ? "\n" : " ");
    }

    void swap_row(mat& A, int i, int j) {
        if(i == j)
            return;
        for(int k = 0; k < A.acol; ++k)
            swap(A.a[i][k], A.a[j][k]);
    }

    // row[k] = zk * row[k] + zi * row[i]
    void add_row(mat& A, int k, int i, int zk, int zi) {
        for(int j = 0; j < A.acol; ++j)
            A.a[k][j] = zk * A.a[k][j] + zi * A.a[i][j];
    }

    void gauss_jordan(mat& A) {
        vector<int> first;
        first.resize(A.row, A.col);
        // row echelon
        int cur = 0;
        for(int k = 0; k < A.col; ++k) {
            int pivot = cur;
            while(pivot < A.row && A.a[pivot][k] == 0)
                ++pivot;
            if(pivot == A.row)
                continue;
            int key = A.a[pivot][k];
            swap_row(A, pivot, cur);
            first[cur] = k;

            for(int i = cur + 1; i < A.row; ++i) {
                int lead = A.a[i][k];
                int t = gcd(abs(lead), abs(key));
                add_row(A, i, cur, key / t, -lead / t);
            }
            ++cur;
        }
        // reduced
        for(int k = A.row - 1; k >= 0; --k) {
            int pivot = first[k];
            if(pivot == A.col)
                continue;
            int key = A.a[k][pivot];

            for(int i = k - 1; i >= 0; --i) {
                int lead = A.a[i][pivot];
                int t = gcd(abs(lead), abs(key));
                add_row(A, i, k, key / t, -lead / t);
            }
        }
        // normalize
        for(int k = 0; k < A.row; ++k) {
            int pivot = first[k];
            if(pivot == A.col)
                continue;

            int t = 0;
            for(int i = 0; i < A.acol; ++i)
                t = gcd(t, abs(A.a[k][i]));
            t *= sgn(A.a[k][pivot]);

            for(int i = 0; i < A.acol; ++i)
                A.a[k][i] /= t;
        }
    }
}

// equation     -> list '=' list
// list         -> component '+' list | component
// component    -> 'e' | flist charge
// flist        -> factor | factor flist
// factor       -> '(' component ')' number | element number
// charge       -> '^' number sign | \eps
// sign         -> '+' | '-'
// number       -> /[1-9][0-9]*/ | \eps
// element      -> /[A-Z][a-z]*/

namespace Alice {
    typedef map<string, int> Elements;

    struct state {
        const char *p;
    };

    void ws(state* s) {
        while(isspace(*s->p)) ++s->p;
    }

    bool match(state* s, char c) {
        ws(s);
        return *s->p == c ? (++s->p, true) : false;
    }

    int number(state* s) {
        ws(s);
        if(!isdigit(*s->p))
            return 1;

        int r = 0;
        while(isdigit(*s->p))
            r = 10 * r + (*s->p++ - '0');
        return r;
    }

    string element(state* s) {
        ws(s);

        string t;
        if(isupper(*s->p))
            t += *s->p++;
        else
            throw "Unexpected EOF";

        while(islower(*s->p))
            t += *s->p++;
        return t;
    }

    Elements component(state* s);

    Elements factor(state* s) {
        if(match(s, '(')) {
            Elements M = component(s);
            if(!match(s, ')')) throw "Expect ')'";

            int k = number(s);
            for(auto& z : M)
                z.second *= k;

            return M;
        } else {
            string t = element(s);
            int k = number(s);

            Elements M;
            M[t] = k;
            return M;
        }
    }

    Elements charge(state* s) {
        Elements M;
        if(match(s, '^')) {
            int k = number(s);
            int sign;

            if(match(s, '+'))
                sign = 1;
            else if(match(s, '-'))
                sign = -1;
            else
                throw "Input error";

            M["e"] = sign * k;
        } else
            M["e"] = 0;

        return M;
    }

    bool predict_factor(state* s) {
        ws(s);
        return isupper(*s->p) || *s->p == '(';
    }

    Elements component(state* s) {
        Elements M;
        if(match(s, 'e')) {
            M["e"] = -1;
        } else {
            do {
                Elements P = factor(s);
                for(auto& z : P)
                    M[z.first] += z.second;
            } while(predict_factor(s));

            Elements P = charge(s);
            for(auto& z : P)
                M[z.first] += z.second;
        }
        return M;
    }

    vector<Elements> list(state* s) {
        vector<Elements> l;
        do
            l.push_back(component(s));
        while(match(s, '+'));
        return l;
    }

    vector<Elements> equation(state* s) {
        vector<Elements> lhs = list(s);
        match(s, '=');
        vector<Elements> rhs = list(s);

        for(auto& M : rhs) {
            for(auto& z : M)
                z.second *= -1;
            lhs.push_back(M);
        }
        return lhs;
    }
}

namespace Marisa {
    Cirno::mat construct(const char *s) {
        Alice::state t;
        t.p = s;
        vector<Alice::Elements> l = Alice::equation(&t);

        map<string, vector<int> > T;
        for(auto& M : l)
            for(auto& z : M)
                T.insert(make_pair(z.first, vector<int>()));

        for(auto& p : T)
            for(auto& M : l)
                p.second.push_back(M[p.first]);

        Cirno::mat A;
        A.row = static_cast<int>(T.size() + 1);
        A.col = static_cast<int>(l.size());
        A.acol = A.col + 1;

        A.a.resize(A.row);
        for(int i = 0; i < A.row; ++i)
            A.a[i].resize(A.acol);

        int k = 0;
        for(auto& p : T) {
            for(int j = 0; j < A.col; ++j)
                A.a[k][j] = p.second[j];
            A.a[k][A.col] = 0;
            ++k;
        }
        A.a[k][0] = 1;
        A.a[k][A.col] = 1;

        return A;
    }

    vector<int> coef(const Cirno::mat& A) {
        {
            int i = 0;
            for(; i < A.row; ++i) {
                int nonzero = 0;
                for(int j = 0; j < A.acol; ++j)
                    if(A.a[i][j] != 0)
                        ++nonzero;

                if(nonzero == 1)
                    throw "No solution";
                if(nonzero == 0)
                    break;
            }
            if(i < A.col)
                throw "Multiple solutions";
        }

        int lcm = 1;
        // lcm(a, b, c) = lcm(lcm(a, b), c) = Xc/(X,c)
        for(int k = 0; k < A.col; ++k) {
            int t = A.a[k][k];
            lcm *= t / Cirno::gcd(lcm, t);
        }

        vector<int> x;
        for(int k = 0; k < A.col; ++k) {
            int t = A.a[k][k];
            x.push_back(lcm / t * A.a[k][A.col]);
        }
        return x;
    }

    string place(const string& t, const vector<int>& x) {
        string::size_type p = 0;
        string s;
        for(int k = 0; k < (int)x.size(); ++k) {
            string::size_type q = string::npos, r = t.find_first_of("+=", p);
            if(r != string::npos)
                q = t.find_first_not_of(" \t\n", r + 1);

            if(q != string::npos) {
                string z = t.substr(p, q - p);
                if(z.find('^') != string::npos && z.find('-') == string::npos
                        && (r = t.find_first_of("+=", q)) != string::npos)
                    q = t.find_first_not_of(" \t\n", r + 1);
            }
            if(x[k] != 1)
                s += to_string(x[k]);
            s += t.substr(p, q - p);
            p = q;
        }
        return s;
    }

    void run(const string& s) {
        try {
            if(s.empty())
                return;

            Cirno::mat A = construct(s.c_str());
            Cirno::gauss_jordan(A);
            vector<int> x = coef(A);

            string ans = place(s, x);
            cout << ans << endl;
        } catch(const char* err) {
            cerr << err << endl;
        }
    }

    void loop() {
        string s;
        while(getline(cin, s))
            run(s);
    }
}

int main() {
    Marisa::loop();
    return 0;
}

