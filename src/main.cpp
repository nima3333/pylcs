#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <iostream>
#include <string.h>
#include <sstream>
#include <algorithm>

using namespace std;
namespace py = pybind11;


vector<string> utf8_split(const string &str){
    vector<string> split;
    int len = str.length();
    int left = 0;
    int right = 1;

    for (int i = 0; i < len; i++){
        if (right >= len || ((str[right] & 0xc0) != 0x80)){
            string s = str.substr(left, right - left);
            split.push_back(s);
            // printf("%s %d %d\n", s.c_str(), left, right);
            left = right;
        }
        right ++;
    }
    return split;
}


// 最长公共子序列（不连续）
int lcs_length_(const string &str1, const string &str2) {
    if (str1 == "" || str2 == "")
        return 0;
    vector<string> s1 = utf8_split(str1);
    vector<string> s2 = utf8_split(str2);
    int m = s1.size();
    int n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    int i, j;
    // printf("%d %d\n", m, n);

    for (i = 0; i <= m; i++) {
        dp[i][0] = 0;
    }
    for (j = 0; j <= n; j++) {
        dp[0][j] = 0;
    }
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                if (dp[i - 1][j] >= dp[i][j - 1])
                    dp[i][j] = dp[i - 1][j];
                else
                    dp[i][j] = dp[i][j-1];
            }
        }
    }
    return dp[m][n];
}

// 最长公共子序列（不连续）
double sw(const string &str1, const string &str2, double s, double w) {
    if (str1 == "" || str2 == "")
        return 0;
    vector<string> s1 = utf8_split(str1);
    vector<string> s2 = utf8_split(str2);
    int n1 = s1.size();
    int n2 = s2.size();

    vector<vector<double>> S_matrix(n1, vector<double>(n2));
    int i, j, t;
    
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
            if(s1[i] == s2[j])
                S_matrix[i][j] = s;
            else
                S_matrix[i][j] = -s;
        }
    }

    vector<vector<double>> H_matrix(n1+1, vector<double>(n2+1));
    vector<double> max2_matrix(n1+1);

    for (j = 1; j <= n2; j++) {
        for (t = 0; t <= n1; t++) {
            max2_matrix[t] = max(max2_matrix[t], H_matrix[t][j - 1]) - w;
        }
        double max1 = 0;
        for (i = 1; i <= n1; i++) {
            max1 = max(H_matrix[i - 1][j], max1) - w;
            double temp = H_matrix[i - 1][j - 1] + S_matrix[i - 1][j - 1];
            H_matrix[i][j] = max( max( max(temp, max1) , max2_matrix[i]), 0.0);
        }
    }

    double max_return = std::numeric_limits<double>::lowest();
    for (const auto& v : H_matrix)
    {
        double current_max = *std::max_element(v.cbegin(), v.cend());
        max_return = max_return < current_max ? current_max : max_return; // max = std::max(current_max, max);
    }

    return max_return;
}

// 最长公共子串（连续）
int lcs2_length_(const string &str1, const string &str2) {
    if (str1 == "" || str2 == "")
        return 0;
    vector<string> s1 = utf8_split(str1);
    vector<string> s2 = utf8_split(str2);
    int m = s1.size();
    int n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    int i, j;
    int max = 0;

    for (i = 0; i <= m; i++) {
        dp[i][0] = 0;
    }
    for (j = 0; j <= n; j++) {
        dp[0][j] = 0;
    }
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if (dp[i][j] > max){
                    max = dp[i][j];
                }
            }
            else {
                dp[i][j] = 0;
            }
        }
    }
    return max;
}


// TODO 返回子序列
int lcs(const string &str1, const string &str2){
    return lcs_length_(str1, str2);
}


// TODO 返回子串
int lcs2(const string &str1, const string &str2){
    return lcs2_length_(str1, str2);
}


// 编辑距离
int levenshtein_distance(const string &str1, const string &str2) {
    if (str1 == "" || str2 == "")
        return 0;
    vector<string> s1 = utf8_split(str1);
    vector<string> s2 = utf8_split(str2);
    int m = s1.size();
    int n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    int i, j;

    for (i = 0; i <= m; i++) {
        dp[i][0] = i;
    }
    for (j = 0; j <= n; j++) {
        dp[0][j] = j;
    }
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + min({dp[i][j - 1], dp[i - 1][j], dp[i - 1][j - 1]});
            }
        }
    }
    return dp[m][n];
}



double lcs_sim_score_1(py::array_t<int> input1, py::array_t<int> input2) {

    // Format of input1 and input2:
    // [ [TYPE(int), TLS_VERSION(int), SIZE(int)],
    //   [TYPE(int), TLS_VERSION(int), SIZE(int)], ...
    // ]

    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();

    int X1 = buf1.shape[0];
    int Y1 = buf1.shape[1];

    int X2 = buf2.shape[0];
    int Y2 = buf2.shape[1];
    
    assert(Y1 == Y2);

    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;

    int m = X1;
    int n = X2;
    vector<vector<double>> dp(m + 1, vector<double>(n + 1));
    int i, j;
    // printf("%d %d\n", m, n);

    // Initialize array
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            dp[i][j] = 0;
        }
    }

    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            // Check if same type
            if (ptr1[(i - 1)*Y1] == ptr2[(j - 1)*Y2]) {
                // Calculate score
                vector<double> score_list;
                score_list.push_back(1.0);
                // Fetch TLS and size
                int size1 = ptr1[(i - 1)*Y1 + 2], size2 = ptr2[(j - 1)*Y2 + 2];
                int tls1 = ptr1[(i - 1)*Y1 + 1], tls2 = ptr2[(j - 1)*Y2 + 1];
                bool is_size = (size1 != -2) and (size2 != -2);
                bool is_tls = (tls1 != -2) and (tls2 != -2);

                if(is_size and size1 != 0 and size2 != 0){
                    int size_max = max(size1, size2);
                    int size_min = min(size1, size2);
					score_list.push_back(1. - (double(size_max - size_min) / double(size_max)));
                }
                else if(is_size and ((size1 == 0) != (size2 == 0))){
                    score_list.push_back(1);
                }
                if(is_tls){
                    score_list.push_back(0.5 + 0.5 * (tls1 == tls2));
                }
                // Average
                double final_score;
                double sumTotal = 0;
                for(int k=0; k < score_list.size(); ++k){
                    // not sure what to put here basically.
                    sumTotal += score_list[k];            
                }
                final_score = sumTotal / score_list.size();
                dp[i][j] = dp[i - 1][j - 1] + final_score;
            } else {
                if (dp[i - 1][j] >= dp[i][j - 1])
                    dp[i][j] = dp[i - 1][j];
                else
                    dp[i][j] = dp[i][j-1];
            }
        }
    }
    return dp[m][n] / max(m, n);
}

double lcs_sim_score_2(py::array_t<int> input1, py::array_t<int> input2) {

    // Format of input1 and input2:
    // [ [TYPE(int), TLS_VERSION(int), SIZE(int)],
    //   [TYPE(int), TLS_VERSION(int), SIZE(int)], ...
    // ]

    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();

    int X1 = buf1.shape[0];
    int Y1 = buf1.shape[1];

    int X2 = buf2.shape[0];
    int Y2 = buf2.shape[1];
    
    assert(Y1 == Y2);

    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;

    int m = X1;
    int n = X2;
    vector<vector<double>> dp(m + 1, vector<double>(n + 1));
    int i, j;
    // printf("%d %d\n", m, n);

    // Initialize array
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            dp[i][j] = 0;
        }
    }

    double maximum;
    
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            // Check if same type
            if (ptr1[(i - 1)*Y1] == ptr2[(j - 1)*Y2]) {
                // Calculate score
                vector<double> score_list;
                score_list.push_back(1.0);
                // Fetch TLS and size
                int size1 = ptr1[(i - 1)*Y1 + 2], size2 = ptr2[(j - 1)*Y2 + 2];
                int tls1 = ptr1[(i - 1)*Y1 + 1], tls2 = ptr2[(j - 1)*Y2 + 1];
                bool is_size = (size1 != -2) and (size2 != -2);
                bool is_tls = (tls1 != -2) and (tls2 != -2);

                if(is_size and size1 != 0 and size2 != 0){
                    int size_max = max(size1, size2);
                    int size_min = min(size1, size2);
					score_list.push_back(1. - (double(size_max - size_min) / double(size_max)));
                }
                else if(is_size and ((size1 == 0) != (size2 == 0))){
                    score_list.push_back(1);
                }
                if(is_tls){
                    score_list.push_back(0.5 + 0.5 * (tls1 == tls2));
                }
                // Average
                double final_score;
                double sumTotal = 0;
                for(int k=0; k < score_list.size(); ++k){
                    // not sure what to put here basically.
                    sumTotal += score_list[k];            
                }
                final_score = sumTotal / score_list.size();
                dp[i][j] = dp[i - 1][j - 1] + final_score;
                if (dp[i][j] > maximum){
                    maximum = dp[i][j];
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }
    return maximum / max(m, n);
}


PYBIND11_MODULE(pylcs, m) {
    m.def("lcs", &lcs, R"pbdoc(
        Longest common subsequence
    )pbdoc");

    m.def("lcs2", &lcs2, R"pbdoc(
        Longest common substring
    )pbdoc");

    m.def("levenshtein_distance", &levenshtein_distance, R"pbdoc(
        Levenshtein Distance of Two Strings
    )pbdoc");

    m.def("edit_distance", &levenshtein_distance, R"pbdoc(
        Same As levenshtein_distance(): Levenshtein Distance of Two Strings
    )pbdoc");

    m.def("custom_lcs", &lcs_sim_score_1, R"pbdoc(
        Custom Longest common subsequence taking into account attributes
    )pbdoc");

    m.def("custom_lcs2", &lcs_sim_score_2, R"pbdoc(
        Custom Longest common subsequence (2) taking into account attributes
    )pbdoc");

    m.def("smith_w", &sw, R"pbdoc(
        Smith_waterman
    )pbdoc");
}
