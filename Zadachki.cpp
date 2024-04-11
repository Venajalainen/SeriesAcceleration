#define _USE_MATH_DEFINES // for C++
#include <iostream>
#include <cstdio>
#include <iomanip>
#include<fstream>
#include <math.h>
#include <vector>
//#include <string.h>
//#include <string>
//#include <cmath>
//#include <stdio.h>
//#include <stdlib.h>
//#include <set>
//#include <cstdlib>
//#include <ctime>
//#include <random>
//#include <cmath>

std::ofstream fout("SUKA.txt");

long double Wynn1(std::vector<long double> s) //From Python of https://github.com/pjlohr/WynnEpsilon
{
    int n = s.size();
    if (n % 2 == 0) {
        s.pop_back();
        n = n - 1;
    }
    std::vector<std::vector<long double>> A(n, std::vector<long double>(n, 0));
    for (int i = 0; i < n; i++) {
        A[i][0] = s[i];
        fout << A[i][0] << " ";
    }
    fout << "\n";

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            if (j == 1) {
                A[i][j] = 1 / (A[i + 1][j - 1] - A[i][j - 1]);
            }
            else {
                A[i][j] = A[i + 1][j - 2] + 1 / (A[i + 1][j - 1] - A[i][j - 1]);
            }
            fout << A[i][j] << " ";
        }
        fout << "\n";
    }
    return A[0][n - 1];
}

long double Wynn2(std::vector<long double> sn) //From Python of https://mathproblems123.wordpress.com/2021/05/21/accelerating-convergence-wynns-algorithm/
    {
    int n = sn.size();
    std::vector<std::vector<long double>> e(n + 1, std::vector<long double>(n + 1, 0.0));
    for (int i = 1; i <= n; i++) {
        e[i][1] = sn[i - 1];
        fout << e[i][1] << " ";
    }
    fout << "\n";

    for (int i = 3; i <= n + 1; i++) {
        for (int j = 3; j <= i; j++) {
            e[i - 1][j - 1] = e[i - 2][j - 3] + 1 / (e[i - 1][j - 2] - e[i - 2][j - 2]);
            fout << e[i-1][j-1] << " ";
        }
        fout << "\n";
    }
    
    if (n % 2 == 0) {
        e.back().pop_back();
    }

    return e.back().back();
}

long double Wynn3(std::vector<long double> A) //ChatGPT
{
    std::vector<long double> v = A;
    std::vector<long double> epsilon;

    if (v.size() % 2 == 0)
    {
        v.pop_back();
    }

    for (int k = 1; k < v.size(); k++) {
        epsilon.push_back(v[k]);
        fout << epsilon[k-1] << " ";
    }
    fout << "\n";

    for (int m = 1; m < v.size() - 1; m++) {
        for (int k = 1; k < v.size() - m; k++) {
            const long double temp = v[k];
            v[k] = epsilon[k - 1] + 1. / (v[k] - v[k + 1]);
            epsilon[k - 1] = temp;
            fout << epsilon[k - 1] << " ";
        }
        fout << "\n";
    }

    for (int i = 0; i < v.size(); i++)
        fout << v[i] << " ";
    fout << "\n";

    return v.back();
}

long double Wynn4(std::vector<long double> A) //DardklordR
{
    //https://www.sciencedirect.com/science/article/pii/S0377042700003551#BIB2
    //https://www.adamponting.com/wynns-epsilon-method/
    int n = A.size();

    std::vector<std::vector<long double>> Eps(n+2, std::vector<long double>(n+1, 0)); //-1 -> 0; 0->1
    //T[n] = (S[n + 1] * S[n - 1] - S[n] ^ 2) / (S[n + 1] - 2S[n] + S[n]);

    fout << "Eps -1 ";
    for (int i = 1; i != n; i++)
    {
        fout << Eps[0][i] << " ";
    }
    fout << '\n';

    fout << "Eps 0  ";
    for (int i = 0; i != n; i++)
    {
        Eps[1][i]=A[i];
        
        fout << Eps[1][i] << " ";
    }
    fout << '\n';

    int k = n;

    for (int r = 1; r != n+1; r++)
    {
        fout << "Eps " << r << "  ";
        for (int i = 0; i != k; i++)
        {
            Eps[r + 1][i] = Eps[r - 1][i + 1] + 1 / (Eps[r][i + 1] - Eps[r][i]);
            fout << Eps[r + 1][i] << " ";
        }
        k--;

        fout << '\n';
    }

    if (n%2 != 0)
    { 
        n = n - 1;
    }
    return Eps[n+1][0];
}

/*
long double Wynn5(std::vector<long double> s)
{

}
*/

std::vector <long double> L4Row(int n)
{
    std::vector <long double> s(n, 0);
    s[0] = 4.0f;
    for (int i = 1; i < n; i++) {
        s[i] = (s[i - 1] + 4.0f * long double(pow(-1, i) * 1.0f / (2.0f * long double (i) + 1.0))); //Ряд Лейбница * 4
    }
    fout << "s = ";
    for (int i = 0; i < n; i++) {
        fout << s[i] << " ";
    }

    fout << "Pi = " << M_PI << '\n';
    fout << "Error before Wynn = " << std::abs(s[n - 1] - M_PI) / M_PI << '\n';
    fout << '\n';
    return s;
}

std::vector <long double> TwonN(int n)
{
    std::vector <long double> s(n, 0);

    s[0] = 0.5f;
    for (int i = 1; i < n; i++) {
        s[i] = (s[i - 1] + (long double(i)+1.0) / long double(pow(2, i+1)));
    }
    fout << "s = ";
    for (int i = 0; i < n; i++) {
        fout << s[i] << " ";
    }

    fout << "Summ = " << 2.0f << '\n';
    fout << "Error before Wynn = " << std::abs(s[n - 1] - 2.0f) / 2.0f << '\n';
    fout << '\n';
    return s;
}

std::vector <long double> Merkator2(int n)
{
    std::vector <long double> s(n, 0);

    s[0] = 1.0f;
    for (int i = 1; i < n; i++) {
        s[i] = (s[i - 1] + long double(pow(-1, i + 1 + 1)) / long double(i+1));
    }
    fout << "s = ";
    for (int i = 0; i < n; i++) {
        fout << s[i] << " ";
    }

    fout << "Summ = " << 0.69314718f << '\n';
    fout << "Error before Wynn = " << std::abs(s[n - 1] - 0.69314718f) / 0.69314718f << '\n';
    fout << '\n';
    return s;
}

void printer(std::vector <long double> s, auto lim)
{
    long double res;

    res = Wynn1(s);
    fout << res << " Error after Wynn 1 = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Wynn2(s);
    fout << res << " Error after Wynn 2 = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Wynn3(s);
    fout << res << " Error after Wynn 3 = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Wynn4(s);
    fout << res << " Error after Wynn 4 = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';
}

int main() {
    //fout << std::setprecision(2);
    printer(L4Row(21),M_PI); //Тут при больше 21 жопа
    fout << "\n" << "_______________________________________________________________________" << "\n";
    printer(TwonN(21), 2.0f); 
    fout << "\n" << "_______________________________________________________________________" << "\n";
    printer(Merkator2(21), 0.69314718f); //Тут лучше нечетно

    return 0;
}


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

