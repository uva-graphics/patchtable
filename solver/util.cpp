
#include "util.h"
#include <algorithm>
#include "array.h"
#define VISUAL_STUDIO_WORKAROUND 1

#include "params.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

string str_replace(string original, const string &before, const string &after) {
//    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
    string ans;
    string::const_iterator end     = original.end();
    string::const_iterator current = original.begin();
    string::const_iterator next    = search(current, end, before.begin(), before.end());
    while (next != end) {
        ans.append(current, next);
        ans.append(after);
        current = next + before.size();
        next = std::search(current, end, before.begin(), before.end());
    }
    ans.append(current, next);
    return ans;
}

bool file_exists(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (f) {
        fclose(f);
        return true;
    }
    return false;
}

vector<string> str_split(const string &s, char c) {
    string sub;
    string::size_type pos = 0, last_pos = 0;
    vector<string> ans;
    bool ok = true;
     
    while (ok) {
        pos = s.find_first_of(c, pos);
        if (pos == string::npos) {
            ok = false;
            pos = s.size();
        }
        sub = s.substr(last_pos, pos - last_pos);
        ans.push_back(sub);
        last_pos = ++pos;
    }
    return ans;
}


string lstrip(string s) {
    //s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    while (s.size() >= 1 && s[0] == ' ') {
        s = s.substr(1);
    }
    return s;
}

string rstrip(string s) {
//    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    while (s.size() >= 1 && s[s.size()-1] == ' ') {
        s = s.substr(0, s.size()-1);
    }
    return s;
}

string strip(string s) {
    return lstrip(rstrip(s));
}

int count(string s, char sub) {
    return (int) std::count(s.begin(), s.end(), sub);
}

string read_binary_file(const string &filename) {
    FILE *file = fopen(filename.c_str(), "rb" );
    if (!file) {
        return string("");
    }
    fseek(file, 0, SEEK_END);
    long n = ftell(file);
    fseek(file, 0, SEEK_SET);
    string text;
    char *buf = new char[n+1];
    buf[n] = '\0';
    if (fread(buf, 1, n, file ) == (unsigned long) n) {
        text = buf;
    }
    fclose(file);
    delete[] buf;
    return text;
}

ProbabilitySampler::ProbabilitySampler(const vector<double> &prob) {
    ASSERT(prob.size(), "expected non-empty prob");
    double sum = 0.0;
    for (int i = 0; i < (int) prob.size(); i++) {
        sum += prob[i];
    }
    ASSERT(sum > 0, "expected sum > 0");
    double scale = 1.0/sum;
    C.resize(prob.size());
    sum = 0.0;
    for (int i = 0; i < (int) prob.size(); i++) {
        sum += prob[i] * scale;
        C[i] = sum;
    }
}

int ProbabilitySampler::sample() {
    double r = rand_f();
    ASSERT(r >= 0 && r < 1, "expected r in [0, 1)");
    int i = lower_bound(C.begin(), C.end(), r) - C.begin();
    if (i < 0) { i = 0; }
    else if (i >= (int) C.size()) { i = C.size()-1; }
    return i;
}

#define NBUF 256

#if !VISUAL_STUDIO_WORKAROUND
string exec_output(const char *cmd) {
    FILE *pipe = popen(cmd, "r");
    if (!pipe) { return "Error"; }
    char buffer[NBUF];
    string result;
    while (!feof(pipe)) {
    	if (fgets(buffer, NBUF, pipe) != NULL) {
    		result += buffer;
        }
    }
    pclose(pipe);
    return result;
}
#endif

double integrate_trapezoid(const vector<double> &x, const vector<double> &y) {
    double ans = 0.0;
    ASSERT(x.size() == y.size(), "x and y sizes should be equal");
    for (int i = 0; i < (int) x.size()-1; i++) {
        double h = x[i+1]-x[i];
        double yavg = (y[i]+y[i+1])*0.5;
        ans += h * yavg;
    }
    return ans;
}

double integrate_box_max(const vector<double> &x, const vector<double> &y) {
    double ans = 0.0;
    ASSERT(x.size() == y.size(), "x and y sizes should be equal");
    for (int i = 0; i < (int) x.size()-1; i++) {
        double h = x[i+1]-x[i];
        double ymax = MAX(y[i], y[i+1]);
        ans += h * ymax;
    }
    return ans;
}

void logf(const char *fmt, ...) {
    static FILE *out = NULL;
    if (out == NULL) {
        out = fopen("out.txt", "wt");
    }
    va_list args;
    va_start(args,fmt);
    char buf[2048];
    vsprintf(buf, fmt, args);
    va_end(args);

    printf("%s", buf);
    fprintf(out, "%s", buf);
    fflush(out);
}

double product(const vector<int> &L) {
    double ans = 1.0;
    for (int i = 0; i < (int) L.size(); i++) {
        ans *= L[i];
    }
    return ans;
}

double sinc(double x) {
    if (x == 0) { return 1.0; }
    double x_arg = M_PI * x;
    return sin(x_arg) / x_arg;
}

vector<int> multiply(const vector<int> &L, double factor) {
    vector<int> ans(L.size());
    for (int i = 0; i < (int) L.size(); i++) {
        ans[i] = int(L[i]*factor+0.5);
    }
    return ans;
}
