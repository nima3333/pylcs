# encoding=utf-8
import numpy as np
import pylcs


def test_lcs():
	assert pylcs.lcs("aaa", "bbb") == 0
	assert pylcs.lcs("aaa", "aabbbaa") == 3
	assert pylcs.lcs("你好", "中国") == 0
	assert pylcs.lcs("aaa你好", "你好呀") == 2


def test_lcs_of_list():
	assert pylcs.lcs_of_list("aaa", ["aabbbaa"] * 10) == [3] * 10
	assert pylcs.lcs_of_list("aaa你好", ["你好呀"] * 10) == [2] * 10


def test_lcs2():
	assert pylcs.lcs2("aaa", "bbb") == 0
	assert pylcs.lcs2("aaa", "aabbbaa") == 2
	assert pylcs.lcs2("你好", "中国") == 0
	assert pylcs.lcs2("aaa你好", "好呀你") == 1


def test_lcs2_of_list():
	assert pylcs.lcs2_of_list("aaa", ["aabbbaa"] * 10) == [2] * 10
	assert pylcs.lcs2_of_list("aaa你好", ["好呀你"] * 10) == [1] * 10


def test_edit_distance():
	assert pylcs.edit_distance("aaa", "bbb") == 3
	assert pylcs.edit_distance("aaa", "aabbbaa") == 4
	assert pylcs.edit_distance("你好", "中国") == 2
	assert pylcs.edit_distance("aaa你好", "你好呀") == 4


def test_edit_distance_of_list():
	assert pylcs.edit_distance_of_list("aaa", ["bbb"] * 10) == [3] * 10
	assert pylcs.edit_distance_of_list("aaa你好", ["你好呀"] * 10) == [4] * 10
	assert pylcs.edit_distance_of_list("aaa你好", ["bbb", "你好呀"] * 10) == [5, 4] * 10


a = np.ones((10, 4)) * 3
b = np.ones((10, 4)) * 3
b[1][0] = 4

def sw(x1, x2, w=0.03, s=0.03):
    n1 = len(x1)
    n2 = len(x2)
    
    S = np.zeros((n1, n2))
    
    for i in range(n1):
        for j in range(n2):
            if x1[i] == x2[j]:
                S[i, j] = s
            else:
                S[i, j] = -s
    
    H = np.zeros((n1 + 1, n2 + 1))
    max2 = np.zeros(n1 + 1)
    
    for j in range(1, n2 + 1):
        for t in range(n1 + 1):
            max2[t] = max(max2[t], H[t, j - 1]) - w
        max1 = 0
        for i in range(1, n1 + 1):
            max1 = max(H[i - 1, j], max1) - w
            H[i, j] = max(max(H[i - 1, j - 1] + S[i - 1, j - 1], max1), max2[i], 0)
    return np.max(H)



a = pylcs.smith_w('aaabaaabaaabaaab', 'aaabaaabaaababba', 0.03, 0.03)
print(a)
b = sw('aaabaaabaaabaaab', 'aaabaaabaaababba')
print(b)