#pragma once
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

vector<vector<double>> operator * (double val, const vector<vector<double>>& vec)
{
    size_t n = vec.size();
    vector<vector<double>> res(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i][j] = val * vec[i][j];
    return res;
}

vector<double> operator * (double val, const vector<double>& vec)
{
    vector<double> res(vec.size());
    size_t n = vec.size();
    for (size_t i = 0; i < n; ++i)
        res[i] = val * vec[i];
    return res;
}

vector<double> operator + (const vector<double>& vec1, const vector<double>& vec2)
{
    vector<double> res(vec1.size());
    size_t n = vec1.size();
    for (size_t i = 0; i < n; ++i)
        res[i] = vec1[i] + vec2[i];
    return res;
}

vector<double> operator - (const vector<double>& vec1, const vector<double>& vec2)
{
    vector<double> res(vec1.size());
    size_t n = vec1.size();
    for (size_t i = 0; i < n; ++i)
        res[i] = vec1[i] - vec2[i];
    return res;
}

vector<double> operator *(const vector<vector<double>>& vec2, const vector<double>& vec1)
{
    size_t n = vec1.size();
    vector<double> res(n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i] += vec2[i][j] * vec1[j];
    return res;
}

double operator *(const vector<double>& vec1, const vector<double>& vec2)
{
    size_t n = vec1.size();
    double res = 0;

    for (int i = 0; i < n; i++)
        res += vec1[i] * vec2[i];

    return res;
}

double norm(const vector<double>& vec)
{
    return sqrt(vec * vec);
}

inline double det(vector<vector<double>> H1)
{
    return H1[0][0] * H1[1][1] - H1[0][1] * H1[1][0];
}




