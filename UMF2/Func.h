#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Func
{
public:
    Func() {}

    double lambda(double u)
    {
        return u + 1;
    }

    double u(double x, double t)
    {
        return x + 2*t;
    }

    double func(double x)
    {
        return 1;
    }
};


