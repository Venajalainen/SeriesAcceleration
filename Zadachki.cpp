#define _USE_MATH_DEFINES // for C++
#include <iostream>
#include <math.h>
#include <string.h>
#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cmath>

double Wynn1(std::vector<double> s) 
{
    int n = s.size();
    if (n % 2 == 0) {
        s.pop_back();
        n = n - 1;
    }
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        A[i][0] = s[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            if (j == 1) {
                A[i][j] = 1 / (A[i + 1][j - 1] - A[i][j - 1]);
            }
            else {
                A[i][j] = A[i + 1][j - 2] + 1 / (A[i + 1][j - 1] - A[i][j - 1]);
            }
        }
    }
    std::cout << A[0][n - 1] << std::endl;
    return A[0][n - 1];
}

double Wynn2(std::vector<double> sn, int n) {
    std::vector<std::vector<double>> e(n + 1, std::vector<double>(n + 1, 0.0));
    for (int i = 1; i <= n; i++) {
        e[i][1] = sn[i - 1];
    }
    
    for (int i = 3; i <= n + 1; i++) {
        for (int j = 3; j <= i; j++) {
            e[i - 1][j - 1] = e[i - 2][j - 3] + 1 / (e[i - 1][j - 2] - e[i - 2][j - 2]);
        }
    }

    std::vector<std::vector<double>> er;

    for (int i = 0; i < n + 1; i++) {
        std::vector<double> row;
        for (int j = 0; j < n + 1; j += 2) {
            row.push_back(e[i][j]);
        }
        er.push_back(row);
    }

    std::cout << er.back().back() << std::endl;
    return er.back().back();
}

int main() {
    int n = 21;
    std::vector<double> s(n, 0);
    s[0] = 4;
    for (int i = 1; i < n; i++) {
        s[i] = (s[i - 1] + 4 * pow(-1, i) * 1 / (2 * i + 1));
    }
    std::cout << "s = ";
    for (int i = 0; i < n; i++) {
        std::cout << s[i] << " ";
    }

    std::cout << "Pi = " << M_PI << std::endl;
    std::cout << "Error before Wynn = " << std::abs(s[n - 1] - M_PI) / M_PI << std::endl;

    std::cout << std::endl;
    double res = Wynn1(s);
    std::cout << "Error after Wynn 1 = " << std::abs(res - M_PI) / M_PI << std::endl;
     res = Wynn2(s, n);
    std::cout << "Error after Wynn 2 = " << std::abs(res - M_PI) / M_PI << std::endl;
    return 0;
}

