#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "vector.h"
#include "Func.h"

using namespace std;

class Matrix : Func
{
public:
    size_t N, timeN, _time;
    double g, a, b, t1, t2, di, f1, f2, aul, dtime;
    string path = "test1/", out = "u";
    vector<double> AuM, AlM, DiM, AuMt, AlMt, DiMt, AuG, AlG, DiG, Au, Al, Di, ai, f, hx, t, q, _u;

    Matrix() {}

    Matrix(vector<double> Q, vector<double> U, size_t _TIME)
    {
        this->_time = _TIME;
        ifstream f;
        f.open(path + "area.txt");
        f >> N >> a >> b >> g;
        f >> t1 >> t2 >> timeN;
        f.close();
        resize();
        this->q = Q;
        this->_u = U;
        step();
        dtime = t[1] - t[0];
        globalMatrixMandG();
        globalMatrix();
        firstBorder();
    } 

    void matrixG(double x1, double x2, double q1, double q2)
    {
        double coefficient = (lambda(q1) + lambda(q2)) / (2 * (x2 - x1));
        di = coefficient;
        aul = -coefficient;
    }

    void matrixM(double x1, double x2)
    {
        double coefficient = (x2 - x1) * g / 6;
        di = 2*coefficient;
        aul = coefficient;
    }

    void vectorF(double x1, double x2)
    {
        double coefficient = ((x2 - x1) / 6);
        f1 = coefficient * (2 * func(x1) + func(x2));
        f2 = coefficient * (func(x1) + 2 * func(x2));
    }

    void globalMatrixMandG()
    {
        for (size_t i = 0; i < N; i++)
        {
            matrixG(hx[i], hx[i + 1], _u[i], _u[i+1]);
            DiG[i] += di; DiG[i + 1] += di; AlG[i] += aul; AuG[i] += aul;
            matrixM(hx[i], hx[i + 1]);
            DiM[i] += di; DiM[i + 1] += di; AlM[i] += aul; AuM[i] += aul;
            DiMt[i] += di / dtime; DiMt[i + 1] += di / dtime; AlMt[i] += aul / dtime; AuMt[i] += aul / dtime;
            vectorF(hx[i], hx[i + 1]);
            f[i] += f1; f[i + 1] += f2;
        }
    }

    void globalMatrix()
    {
        for (size_t i = 0; i < N; i++)
        {
            Di[i] = DiMt[i] + DiG[i];
            Au[i] = AuMt[i] + AuG[i];
            Al[i] = AlMt[i] + AlG[i];
        }
        f[0] += DiMt[0] * q[0];
        for (size_t i = 1; i < N + 1; i++)
        {
            f[i] += DiMt[i] * q[i];
            f[i] += AlMt[i - 1] * q[i];
            f[i - 1] += AuMt[i - 1] * q[i - 1];
        }
    }

    void firstBorder()
    {
        Di[0] = 1; Au[0] = 0; f[0] = u(hx[0], t[_time]);
        Di[N] = 1; Al[N-1] = 0; f[N] = u(hx[N], t[_time]);
    }

    void resize()
    {
        _u.resize(N + 1);
        t.resize(timeN + 1);
        hx.resize(N + 1);
        q.resize(N + 1);
        ai.resize(N + 2);
        f.resize(N+1);
        Au.resize(N);
        Al.resize(N);
        Di.resize(N + 1);
        AuM.resize(N);
        AlM.resize(N);
        DiM.resize(N + 1);
        AuMt.resize(N);
        AlMt.resize(N);
        DiMt.resize(N + 1);
        AuG.resize(N);
        AlG.resize(N);
        DiG.resize(N + 1);
    }

    void step()
    {
        double h = (b - a) / N; double ht = (t2 - t1) / timeN;
        hx[0] = a; t[0] = t1;
        for (size_t i = 1; i < N+1; i++)
            hx[i] = hx[i - 1] + h;
        for (size_t i = 1; i < timeN + 1; i++)
            t[i] = t[i - 1] + ht;
        ai[0] = 0; ai[1] = 0;
        for (size_t i = 2; i < N + 2; i++)
        {
            ai[i] = i - 1;
        }

    }

};