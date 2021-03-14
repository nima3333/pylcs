#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <iostream>
#include <string.h>
#include <sstream>

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


vector<int> lcs_of_list(const string &str1, vector<string> &str_list){
    int size = str_list.size();
    vector<int> ls(size);
    for (int i = 0; i < size; i++){
        int l = lcs(str1, str_list[i]);
        ls[i] = l;
    }
    return ls;
}


vector<int> lcs2_of_list(const string &str1, vector<string> &str_list){
    int size = str_list.size();
    vector<int> ls(size);
    for (int i = 0; i < size; i++){
        int l = lcs2(str1, str_list[i]);
        ls[i] = l;
    }
    return ls;
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


vector<int> levenshtein_distance_of_list(const string &str1, vector<string> &str_list){
    int size = str_list.size();
    vector<int> ls(size);
    for (int i = 0; i < size; i++){
        int l = levenshtein_distance(str1, str_list[i]);
        ls[i] = l;
    }
    return ls;
}

py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2) {
  py::buffer_info buf1 = input1.request();
  py::buffer_info buf2 = input2.request();

  if (buf1.size != buf2.size) {
    throw std::runtime_error("Input shapes must match");
  }

  /*  allocate the buffer */
  py::array_t<double> result = py::array_t<double>(buf1.size);

  py::buffer_info buf3 = result.request();

  double *ptr1 = (double *) buf1.ptr,
         *ptr2 = (double *) buf2.ptr,
         *ptr3 = (double *) buf3.ptr;
  int X = buf1.shape[0];
  int Y = buf1.shape[1];

  for (size_t idx = 0; idx < X; idx++) {
    for (size_t idy = 0; idy < Y; idy++) {
      ptr3[idx*Y + idy] = ptr1[idx*Y+ idy] + ptr2[idx*Y+ idy];
    }
  }
 
  // reshape array to match input shape
  result.resize({X,Y});

  return result;
}


double lcs_length_2(py::array_t<int> input1, py::array_t<int> input2) {

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
    
    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;

    int m = X1;
    int n = X2;
    vector<vector<double>> dp(m + 1, vector<double>(n + 1));
    int i, j;
    // printf("%d %d\n", m, n);

    // Initialize array
    for (i = 0; i <= m; i++) {
        dp[i][0] = 0;
    }
    for (j = 0; j <= n; j++) {
        dp[0][j] = 0;
    }

    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            // Check if same type
            if (ptr1[(i - 1)*Y1] == ptr2[(j - 1)*Y2]) {
                // Calculate score
                vector<double> score_list;
                // Fetch TLS and size
                int size1 = ptr1[(i - 1)*Y1 + 2], size2 = ptr2[(j - 1)*Y2 + 2];
                int tls1 = ptr1[(i - 1)*Y1 + 1], tls2 = ptr2[(j - 1)*Y2 + 1];
                bool is_size = (size1 != -2) and (size2 != -2);
                bool is_tls = (tls1 != -2) and (tls2 != -2);

                score_list.push_back(1.0);
                if(is_size and size1 != 0 and size2 != 0){
                    size1, size2 = max(size1, size2), min(size1, size2);
					score_list.push_back(1 - (size1 - size2) / size1);
                }
                else if(is_size and (size1 == 0 ^ size2 == 0)){
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
    return dp[m][n];
}

/*
double lcs(py::array_t<int> input1, py::array_t<int> input2) {
  py::buffer_info buf1 = input1.request();
  py::buffer_info buf2 = input2.request();

  if (buf1.size != buf2.size) {
    throw std::runtime_error("Input shapes must match");
  }

  py::array_t<double> result = py::array_t<double>(buf1.size);

  py::buffer_info buf3 = result.request();

  double *ptr1 = (double *) buf1.ptr,
         *ptr2 = (double *) buf2.ptr,
         *ptr3 = (double *) buf3.ptr;
  int X = buf1.shape[0];
  int Y = buf1.shape[1];

  for (size_t idx = 0; idx < X; idx++) {
    for (size_t idy = 0; idy < Y; idy++) {
      ptr3[idx*Y + idy] = ptr1[idx*Y+ idy] + ptr2[idx*Y+ idy];
    }
  }
 
  // reshape array to match input shape
  result.resize({X,Y});

  return result;
}
*/

PYBIND11_MODULE(pylcs, m) {
    m.def("lcs", &lcs, R"pbdoc(
        Longest common subsequence
    )pbdoc");

    m.def("lcs_of_list", &lcs_of_list, R"pbdoc(
        Longest common subsequence of list
    )pbdoc");

    m.def("lcs2", &lcs2, R"pbdoc(
        Longest common substring
    )pbdoc");

    m.def("lcs2_of_list", &lcs2_of_list, R"pbdoc(
        Longest common substring of list
    )pbdoc");

    m.def("levenshtein_distance", &levenshtein_distance, R"pbdoc(
        Levenshtein Distance of Two Strings
    )pbdoc");

    m.def("edit_distance", &levenshtein_distance, R"pbdoc(
        Same As levenshtein_distance(): Levenshtein Distance of Two Strings
    )pbdoc");

    m.def("levenshtein_distance_of_list", &levenshtein_distance_of_list, R"pbdoc(
        Levenshtein Distance of one string to a list of strings
    )pbdoc");

    m.def("edit_distance_of_list", &levenshtein_distance_of_list, R"pbdoc(
        Levenshtein Distance of one string to a list of strings
    )pbdoc");

    m.def("add_arrays", &add_arrays, R"pbdoc(
        Add arrays
    )pbdoc");

    m.def("test", &lcs_length_2, R"pbdoc(
        Test
    )pbdoc");

}
