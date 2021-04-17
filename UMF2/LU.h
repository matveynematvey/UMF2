#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "Func.h"
#include "Matrix.h"

using namespace std;

class LU
{
public:
    Matrix mat;
    double size;

    LU() {}
    LU(vector<double> q, vector<double> u, size_t time)
    {
        mat = Matrix(q, u, time);
        size = mat.N + 1;
    }

    void LY()
    {
        for (int i = 0; i < size; i++)
        {
            int i0 = mat.ai[i + 0], i1 = mat.ai[i + 1];

            double s = 0;

            for (int j = i - (i1 - i0), k = i0; j < i; j++, k++)
                s += mat.f[j] * mat.Al[k];

            mat.f[i] = (mat.f[i] - s) / mat.Di[i];
        }
    }

    void UX()
    {
        for (int i = size - 1; i >= 0; i--)
        {
            int i0 = mat.ai[i + 0], i1 = mat.ai[i + 1];

            double xi = mat.f[i];
            for (int j = i - (i1 - i0), k = i0; j < i; j++, k++)
                mat.f[j] -= xi * mat.Au[k];

            mat.f[i] = xi;
        }
    }

    void LUdec()
    {
        for (int i = 0; i < size; i++)
        {
            int i0 = mat.ai[i + 0], i1 = mat.ai[i + 1];

            double sd = 0;
            for (int j = i - (i1 - i0), k = i0; j < i; j++, k++)
            {
                double sl = 0, su = 0;
                int j0 = mat.ai[j + 0], j1 = mat.ai[j + 1];
                int kol_i = k - i0, kol_j = j1 - j0;
                int kol_r = kol_i - kol_j, ki = i0, kj = j0;

                if (kol_r > 0)
                    ki += kol_r;
                else
                    kj -= kol_r;

                for (; ki < k; ki++, kj++)
                {
                    sl += mat.Al[ki] * mat.Au[kj];
                    su += mat.Al[kj] * mat.Au[ki];
                }

                mat.Al[k] = mat.Al[k] - sl;
                mat.Au[k] = mat.Au[k] - su;
                mat.Au[k] /= mat.Di[j];

                sd += mat.Al[k] * mat.Au[k];
            }

            mat.Di[i] = mat.Di[i] - sd;
        }
    }
};
