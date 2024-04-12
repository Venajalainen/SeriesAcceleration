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

std::ofstream fout("Vivod.txt");

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
            v[k] = epsilon[k - 1] + 1.0f/ (v[k] - v[k + 1]);
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
    //https://shareok.org/bitstream/handle/11244/15792/Thesis-1986-C532s.pdf?sequence=1&isAllowed=y
    //https://www.codeproject.com/Articles/5353031/Summation-of-Series-with-Convergence-Acceleration
    int n = A.size();

    std::vector<std::vector<long double>> Eps(n+2, std::vector<long double>(n + 1, 0)); //-1 -> 0; 0->1

    fout << "Eps -1 ";
    for (int i = 1; i != n; i++)
    {
        fout << Eps[0][i] << " ";
    }
    fout << '\n';

    fout << "Eps 0  ";
    for (int i = 0; i != n; i++)
    {
        
        Eps[1][i] = A[i];
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

long double Aitken(std::vector<long double> x) //Aitken's delta-squared process
{
    int n = x.size();
    std::vector<long double> A(n-2);
    for (int i = 0; i < n - 2; i++)
    {
        A[i] = (x[i] * x[i + 2] - x[i + 1] * x[i + 1]) / (x[i] + x[i + 2] - 2.0f * x[i + 1]);
        fout << A[i] << " ";
    }
    fout << '\n';
    return A.back();
}

long double Brenzinski1(std::vector<long double> s) //DarklordR
{
    int n = s.size();

    std::vector<std::vector<long double>> v(n + 4, std::vector<long double>(n + 1, 0)); //-1 -> 0; 0->1

    fout << "v 0";
    for (int i = 0; i != n; i++)
    {
        v[2][i] = s[i];
        fout << v[2][i] << " ";
    }
    fout << '\n';

    for (int k = 3; k < n + 2; k++)
    {
        for (int i = 0; i < n-2; i++)
        {
            v[2 * k][i] = v[2 * k - 1][i + 1] + 1.0f / (v[2 * k][i + 1] - v[2 * k][i]);
            v[2 * k + 1][i] = v[2 * k][i + 1] + (v[2*k][i+2]-v[2*k][i+1]) * (v[2 * k + 1][i + 2] - v[2 * k + 1][i + 1]) / (v[2 * k][i + 2] - v[2 * k][i + 1]);
        }
    }
    return v[0][0];

}

long double GetBinCoeff(int N, int K) //https://www.cyberforum.ru/post347832.html
{
    return ((N < K) ? 0 : ((K == 0) ? 1 : ((N - K + 1) / double(K) * GetBinCoeff(N, K - 1))));
}

long double Levin(std::vector<long double> s) //DarklordR C# translate from https://www.codeproject.com/Articles/5353031/Summation-of-Series-with-Convergence-Acceleration
{
    long double numerator = 0;
    long double denominator = 0;
    int z = s.size();
    int k = z / 2 - 1; //k порядок алгоритма
    int n = z / 2; //n нная точка аппроксимации
    //k+n <= z-1!!!!!
    
    for (int j = 0; j < n; j++)
    {
        long double rest = std::pow(-1.0f, j) * GetBinCoeff(k,j);

        long double C_jkn_U = std::pow((long double)(n + j + 1), k - 1);
        long double C_jkn_L = std::pow((long double)(n + k + 1), k - 1);

        long double C_njk = C_jkn_U / C_jkn_L;

        long double S_nj = s[(int)n + (int)j];

        // t transform that calculates a_n
        double g_n = s[(int)n + j] - s[(int)n + j - 1];
        // u transform that calculates (n+k) * a_n
        // g_n = (n + k)* (S_n[(int)n + j] - S_n[(int)n + j - 1]);

        numerator += rest * C_njk * S_nj / g_n;
        denominator += rest * C_njk / g_n;
    }

    return numerator / denominator;
}


//////////////////////////////////////////////////////////////

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
    fout << "Error before = " << std::abs(s[n - 1] - M_PI) / M_PI << '\n';
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
    fout << "Error before = " << std::abs(s[n - 1] - 2.0f) / 2.0f << '\n';
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
    fout << "Error before = " << std::abs(s[n - 1] - 0.69314718f) / 0.69314718f << '\n';
    fout << '\n';
    return s;
}

void printer(std::vector <long double> s, auto lim)
{
    std::cout << "Lim N " << lim << '\n';

    long double res;

    //res = Wynn1(s); //Т.к. я полностью повторил 
    //fout << res << " Error after Wynn 1 = " << std::abs(res - lim) / lim << '\n';
    //fout << '\n';

    //res = Wynn2(s); //Менее эффективен, чем Wynn 1 и 4
    //std::cout << "Result = " << res << " Error after Wynn 2 = " << std::abs(res - lim) / lim << '\n';
    //fout << '\n';

    res = Wynn3(s);
    std::cout << "Result = "  << res << " Error after GPT Wynn = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Wynn4(s);
    std::cout << "Result = " << res << " Error after DLR Wynn = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Aitken(s);
    std::cout << "Result = " << res << " Error after DLR Aitken = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    res = Levin(s);
    std::cout << "Result = " << res << " Error after DLR Levin = " << std::abs(res - lim) / lim << '\n';
    fout << '\n';

    //res = Levin2(s);
    //std::cout << "Result = " << res << " Error after GPT 2 Levin = " << std::abs(res - lim) / lim << '\n';

    //res = Brenzinski1(s);
    //std::cout << "Result = " << res << " Error after GPT Brenzinski1 = " << std::abs(res - lim) / lim << '\n';
}

int main() {
    fout << std::setprecision(10);

    int k = 0;
    int i = 6;
    while (i < 101)
    {
        std::cout << "Summ N " << i << '\n';
        printer(L4Row(i), M_PI); //Тут при больше 21 жопа
        fout << "\n" << "_______________________________________________________________________" << "\n";
        
        i += 5 ^ k;
        k++;
    }

    fout << "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";
    std::cout << "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";

    i = 6; k = 0;
    while (i < 101)
    {
        std::cout << "Summ N " << i << '\n';
        printer(TwonN(i), 2.0f);
        fout << "\n" << "_______________________________________________________________________" << "\n";
        
        i += 5 ^ k;
        k++;
    }

    fout << "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";
    std::cout << "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";

    i = 6; k = 0;
    while (i < 101)
    {
        std::cout << "Summ N " << i << '\n';
        printer(Merkator2(i), 0.69314718f); //Тут лучше нечетно
        fout << "\n" << "_______________________________________________________________________" << "\n";
        
        i += 5 ^ k;
        k++;
    }

    return 0;
}

