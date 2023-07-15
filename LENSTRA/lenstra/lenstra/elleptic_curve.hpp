#pragma once
#pragma comment (lib, "gmp-10.dll")
#pragma comment (lib, "gmpxx-4.dll")
#pragma warning( disable : 4127 )
#pragma warning( disable : 4146 )

#include <iostream>
#include <random>
#include <limits>
#include <chrono>
#include <mutex>
#include <map>
#include <iomanip>
#include <sstream>
#include <string>
#include <iterator>

#include "gmpxx.h"


#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/independent_bits.hpp>

#define MULT_THS 5
#define INF_POINT POINT{-1,-1}

using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::mutex;
using std::pair;
using std::lock_guard;
using std::multimap;
using std::back_inserter;

using boost::integer::mod_inverse;
using namespace boost::multiprecision;

typedef cpp_int_backend<8192> my_backend;
typedef number<my_backend> mcpp_int;

class timer {
public:
    std::chrono::steady_clock::time_point start, end;
    int unit;
    // 1 - seconds
    // 2 - milliseconds
    // 3 - microseconds
    timer() {
        start = std::chrono::steady_clock::now();
        unit = 2;
    }
    timer(int _unit) {
        start = std::chrono::steady_clock::now();
        unit = _unit;
    }
    ~timer()
    {
        string unit_s;
        end = std::chrono::steady_clock::now();

        long long exec_time;
        if (unit == 1)
        {
            exec_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
            unit_s = unit_s = " seconds.";
        }
        else if (unit == 2)
        {
            exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            unit_s = " millisec.";
        }
        else if (unit == 3)
        {
            exec_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            unit_s = " microsec.";
        }
        else
        {
            exec_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
            unit_s = unit_s = " seconds.";
        }


        cout << "\n+------------------------------------------------+\n";
        string out_str = "| Execution time = " + std::to_string(exec_time) + unit_s;
        for (int i = out_str.size(); i < 49; i++)
            out_str += " ";
        out_str += "|\n";
        cout << out_str;
        cout << "+------------------------------------------------+\n";
    }

    //Возвращает кол-во времени в миллисекундах, которое прошло с момента создания таймера
    unsigned long long get_current_time()
    {
        end = std::chrono::steady_clock::now();
        unsigned long long exec_time;
        string unit_s;
        if (unit == 1)
        {
            exec_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
            unit_s = unit_s = " seconds.";
        }
        else if (unit == 2)
        {
            exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            unit_s = " millisec.";
        }
        else if (unit == 3)
        {
            exec_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            unit_s = " microsec.";
        }
        else
        {
            exec_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
            unit_s = unit_s = " seconds.";
        }

        return exec_time;
    }
};

class POINT
{
public:
    mcpp_int x, y;

    POINT(mcpp_int x, mcpp_int y) : x(x), y(y) {};
    POINT() : x(0), y(0) {};

    //Конструктор копирования
    POINT(const POINT& other)
    {
        x = other.x;
        y = other.y;
    }


    /*void operator*=(const mcpp_int scalar)
    {
        POINT temp;
        temp = *this;

        for (mcpp_int i = 0; i < scalar; i++)
        {
            if (i == 0)
                *this = *this + *this;
            else
                *this = *this + temp;
        }

        return;
    }*/

    //Приведение точки к заданному полю
    POINT operator%(const mcpp_int mod)
    {
        POINT result = *this;
        if (result.x != 0)
        {
            result.x = result.x % mod;
            if (result.x < 0)
                result.x = mod + result.x;
        }
        if (result.y != 0)
        {
            result.y = result.y % mod;
            if (result.y < 0)
                result.y = mod + result.y;
        }
        return result;
    }

    void operator=(const POINT& b)
    {
        this->x = b.x;
        this->y = b.y;
        return;
    }

    bool operator==(const POINT& b)
    {
        return this->x == b.x and this->y == b.y;
    }

    bool operator!=(const POINT& b)
    {
        return !(*this == b);
    }
};

static bool bad_point = false;

class ELLEPTIC_CURVE
{
public:
    ELLEPTIC_CURVE() { a = 0; b = 0; x = 0; y = 0; limit = 0; orderCurve = 0; orderBP = 0;};
    ELLEPTIC_CURVE(mcpp_int a, mcpp_int b, mcpp_int x, mcpp_int y, mcpp_int limit, mcpp_int orderC, mcpp_int orderBp) : a(a), b(b), x(x), y(y), limit(limit), orderCurve(orderC), orderBP(orderBp) {};
    ELLEPTIC_CURVE(const ELLEPTIC_CURVE& other)
    {
        basePoint = other.basePoint;
        allPoint = other.allPoint;
        limit = other.limit;
        a = other.a;
        b = other.b;
        x = other.x;
        y = other.y;
        orderCurve = other.orderCurve;
        orderBP = other.orderBP;
    }

    ELLEPTIC_CURVE& operator= (const ELLEPTIC_CURVE& other)
    {
        if (this == &other)
            return *this;
        a = other.a;
        b = other.b;
        x = other.x;
        y = other.y;
        limit = other.limit;
        orderCurve = other.orderCurve;
        orderBP = other.orderBP;
        allPoint = other.allPoint;
        basePoint = other.basePoint;
    }

    POINT basePoint;
    vector<POINT> allPoint;
    mcpp_int limit, a, b, x, y, orderCurve, orderBP, cofactor;
    mutex mt;

    vector<int> empty_numbers{ 10001521, 10000019, 10000789, 14843, 29207, 27947, 33311, 34031, 38149, 41959, 102397, 99377, 99409, 99131, 99793, 94723, 88951, 82487, 78167,
        77747, 77509, 77849, 77863, 77867, 66271, 66293, 66301, 66337, 66343, 66347, 66359, 66361, 66373, 66377, 65777, 65789, 65809, 65827, 65831, 65837, 65839, 65843, 65851, 65867,
        60821, 60859, 60869, 60887, 60889, 60899, 60901, 60913, 60917, 60919 };


    void find_order_BP()
    {
        if (orderBP != 0)
            return;
        mcpp_int order = 1, counter = 1;
        POINT tmp = basePoint;
        while (true)
        {
            tmp = plus(tmp, basePoint);
            counter++;
            if (tmp == POINT{ -1,-1 })
            {
                orderBP = counter;
                break;
            }
        }
        return;
    }

    // Генерирует случайную базовую точку определенную на ЭК
    POINT gen_base_point()
    {
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        boost::random::mt19937 gen(seed);
        boost::random::uniform_int_distribution<mcpp_int> dist(1, limit);
genBP:
        basePoint.x = dist(gen);
        mcpp_int _y;
        _y = pow(basePoint.x, 3) + basePoint.x * a + b;
        _y = mod_inverse(_y, limit);

        if (_y == 0)
            goto genBP;
        else
            basePoint.y = _y;

        return basePoint;
    }

    //Находит порядок ЭК и все точки заданные на ней
    unsigned long long findAllPoint()
    {
        if (!allPoint.empty())
            return 0;
        timer t;
        //allPoint.reserve(UINT_MAX/2);

        POINT temp;

        mcpp_int _y, tmp;
        pair<mcpp_int, mcpp_int> _pair;

        multimap<mcpp_int, mcpp_int> deductions;

        //Находим все y^2
        cout << "\nLIMIT = " << limit;
        for (mcpp_int i = 1; i < limit; i++)
        {
            tmp = pow(i, 2);
            tmp %= limit;
            _pair = make_pair(tmp, i);
            deductions.insert(_pair);
        }

        //Находим все x и формируем точки
        for (mcpp_int i = 0; i < limit; i++)
        {
            temp.x = i;
            _y = pow(i, 3) + a * i + b;
            _y %= limit;
            //Условие нужно проверить
            if (gcd(_y, limit) == 1)
            {
                auto it = deductions.find(_y);
                if (it != deductions.end())
                {
                    while (it != deductions.end() and it->first == _y)
                    {
                        temp.y = _y;
                        if (basePoint == POINT{ 0,0 } /*and temp.x != 0 and temp.y != 0*/)
                            basePoint = temp;
                        //temp.x = i;
                        //temp.y = it->second;
                        //temp.ptrEC = this;
                        //allPoint.push_back(temp);

                        it++;
                        orderCurve++;
                    }
                }
            }
            else
                continue;
        }
        return t.get_current_time();
    }

    // Задает параметры ЭК согласно рекомендациям NIST 
    void set_default_ec()
    {
        mcpp_int _limit("6277101735386680763835789423207666416083908700390324961279");
        //mcpp_int _limit("6277101");
        limit = _limit;

        a = -3;
        a = a % limit;
        mcpp_int _b("2455155546008943817740293915197451784769108058161191238065");
        //mcpp_int _b("24551");
        b = _b;

        mcpp_int _x("602046282375688656758213480587526111916698976636884684818"),
            _y("174050332293622031404857552280219410364023488927386650641"),
            _orderBP("6277101735386680763835789423176059013767194773182842284081");
        /*mcpp_int _x("60204"),
            _y("17405"),
            _orderBP("6277");*/

        basePoint.x = _x;
        basePoint.y = _y;
        orderBP = _orderBP;
    }


    //Проверяет, определена ли точка на ЭК
    bool indeterminate_point_in_addition(POINT point)
    {
        mcpp_int _x = point.x, _y = point.y, _a = a, _b = b, tmp1, tmp2;
        tmp2 = (pow(_x, 3) + _a * _x + _b) % limit;
        tmp1 = mod_inverse(tmp2, limit);
        if (tmp1 != 0)
        {
            if (tmp1 != tmp2)
                return false;
            else
                return true;
        }
        else
            return false;
    }

    //Сложение двух точек на кривой
    POINT plus(POINT aP, POINT bP)
    {
        if (aP == INF_POINT and bP != INF_POINT)
            return bP;
        else if (aP != INF_POINT and bP == INF_POINT)
            return aP;
        else if (aP == INF_POINT and bP == INF_POINT)
            return INF_POINT;

        int switch_on;
        if (aP == bP)
            switch_on = 1;
        else
            switch_on = 0;

        POINT res;
        mcpp_int tmp, tmp2, lambda;

        switch (switch_on)
        {
        case 0:
            tmp = bP.x - aP.x;
            tmp = tmp % limit;
            if (tmp < 0)
                tmp = limit + tmp;

            tmp2 = mod_inverse(tmp, limit);
            if (tmp2 == 0)
                return INF_POINT;
            if (((aP.x == 0) and bP.x == 0) or (bP.x == aP.x))
            {
                POINT answer(-1, -1);
                return answer;
            }
            else {
                lambda = ((bP.y - aP.y) * tmp2)%limit;
                res.x = (pow(lambda, 2) - aP.x - bP.x) % limit;
                if (res.x < 0)
                    res.x = limit + res.x;
                res.y = (lambda * (aP.x - res.x) - aP.y) % limit;
                if (res.y < 0)
                    res.y = limit + res.y;
            }
            break;
        case 1:
            if (aP.y == 0)
                return INF_POINT;
            else {
                tmp = 2 * aP.y;

                tmp = tmp % limit;
                if (tmp < 0)
                    tmp = limit + tmp;
                tmp2 = mod_inverse(tmp, limit);
                if (tmp2 == 0)
                    return INF_POINT;
                lambda = (((pow(aP.x, 2)) * 3 + a) * tmp2) % limit;
                res.x = (pow(lambda, 2) - aP.x * 2) % limit;
                if (res.x < 0)
                    res.x = limit + res.x;
                res.y = (lambda * (aP.x - res.x) - aP.y) % limit;
                if (res.y < 0)
                    res.y = limit + res.y;
            }
            break;
        }

        return res;
    }

    POINT stupidmult(POINT p, mcpp_int scalar)
    {
        POINT res = p;

        for (mcpp_int i = 0; i < scalar - 1; i++)
        {
            res = plus(res, p);
            if (res == INF_POINT)
                return res;
        }
        return res;
    }

    //Умножение точки на число (методом удвоения точки)
    POINT mult(POINT p, mcpp_int scalar)
    {
        //timer t;
        POINT result(0, 0), tmp = p;
        result = p;

        mcpp_int new_scalar = 1,
            counter = 0;

        while (new_scalar < scalar)
        {
            result = plus(result, result);
            if (result == INF_POINT)
                return result;

            new_scalar *= 2;
            if (new_scalar*2 > scalar)
                break;
        }

        for (mcpp_int i = new_scalar; i < scalar; i++)
        {
            result = plus(result, p);
            if (result == INF_POINT)
                return result;
        }

        return result;
    }

    void _mult(POINT res, POINT p, mcpp_int scalar)
    {
        POINT result(0, 0), tmp = p;
        result = p;

        mcpp_int new_scalar = 1,
            counter = 0;

        while (new_scalar < scalar)
        {
            result = plus(result, result);
            if (result == INF_POINT)
            {
                res = result;
                return;
            }

            new_scalar *= 2;
            if (new_scalar * 2 > scalar)
                break;
        }

        for (mcpp_int i = new_scalar; i < scalar; i++)
        {
            result = plus(result, p);
            if (result == INF_POINT)
            {
                res = result;
                return;
            }
        }

        res = result;
        return;
    }

    POINT mt_mult(POINT p, mcpp_int scalar)
    {
        vector<std::thread> ths;
        vector<mcpp_int> step_box;
        vector<POINT> subRes;
        subRes.resize(10);
        POINT res, resP;

        // 12 zero
        if (scalar > 100000000)
        {
            cout << "\nMultithread mult\n";
            //bad_point = false;
            mcpp_int step = scalar / MULT_THS, remains = scalar % MULT_THS;
            //cout << "\nSource scalar = " << scalar << "\nNew scalar in thread = " << step << endl;
            for (int i = 0; i < MULT_THS; i++)
            {
                ths.push_back(std::thread([&]() {
                    _mult(std::ref(subRes[i]), p, step);
                    }));
            }
            for (int i = 0; i < MULT_THS; i++)
                ths[i].join();

            
            resP = subRes[0];
            if (resP == INF_POINT)
                return INF_POINT;
            for (int i = 1; i < MULT_THS; i++)
            {
                if (subRes[i] == INF_POINT)
                    return INF_POINT;
                resP = plus(resP, subRes[i]);
            }

            if (remains != 0)
            {
                for (mcpp_int i = 0; i < remains; i++)
                {
                    resP = plus(resP, p);
                }
            }

        }
        else
        {
            cout << "\nSimple mult\n";
            resP = mult(p, scalar);
        }

        return resP;
    }

    POINT fastmult(POINT p, mcpp_int scalar)
    {
        POINT P = INF_POINT, R = p;
        unsigned long long bit_size;

        // Переводим число в двоичный формат
        mpz_class x(scalar.str().c_str());
        size_t sz = mpz_sizeinbase(x.get_mpz_t(), 2);
        char* buffer = new char[sz];
        mpz_get_str(buffer, 2, x.get_mpz_t());

        string bin_scalar(buffer);

        bool fl = false;
        for (long long i = sz-1; i >= 0; i--)
        {
            if (bin_scalar[i] == '0')
            {
                R = plus(R, R);
                if (R == INF_POINT)
                    return R;
            }
            else {
                R = plus(P, R);
                if (R == INF_POINT)
                    return R;
            }
        }

        if (!indeterminate_point_in_addition(R))
            return INF_POINT;
        else
            return R;
    }
};

static std::vector<ELLEPTIC_CURVE> badEC;