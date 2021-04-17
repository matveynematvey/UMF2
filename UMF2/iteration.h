#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "LU.h"
#include "Func.h"

using namespace std;

class Iter : Func
{
public:
    LU lu;
    double time;
    vector<double> q, _u;
    Iter()
    {
        q.resize(11);
        _u.resize(11);
        lu = LU(q, _u, 0);
        _u.resize(11);
        for (size_t i = 0; i < lu.size; i++)
        {
            q[i] = u(lu.mat.hx[i], lu.mat.t[0]);
        }
        for (size_t j = 1; j < 4; j++)
        {
            time = lu.mat.t[j];
            for (size_t i = 0; i < 100; i++)
            {
                lu = LU(q, _u, j);
                lu.LUdec();
                lu.LY();
                lu.UX();
                _u = lu.mat.f;
            }
            q = _u;
            cout << "iter" << endl;
            for (size_t j = 0; j < lu.size; j++)
            {
                cout << scientific << q[j] << "     " << abs(q[j] - u(lu.mat.hx[j], time)) << endl;
            }
            cout << endl;
        }
    }

};



