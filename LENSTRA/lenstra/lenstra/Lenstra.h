#pragma once

#include <fstream>
#include <thread>
#include <mutex>
#include <map>
#include <future>
#include <condition_variable>

#include "elleptic_curve.hpp"


#define MAIN_THREADS_MAX 5
#define DIVISOR_THREADS_MAX 2

#define MT 1
#define SMT 2
#define ST 3

static bool KEEP_GOING = true;
static int th_id_counter = 0;

static bool FINDED_DIV = false;

static std::mutex mt;


using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::multimap;
using std::pair;
using std::uniform_int_distribution;
using std::random_device;
using std::mt19937;
using std::mutex;
using std::thread;


class Lenstra
{
private:
    unsigned th_id;
    mcpp_int input_N, B, iteration;
    ELLEPTIC_CURVE ec;
    POINT P;
    vector<POINT> allPoint;

    void main_func(int FLAG);
    mcpp_int intermediate_divisor(mcpp_int input, int FLAG);

    void init_b_value(mcpp_int input);
    mcpp_int generate_ec_coefficient(mcpp_int input);
    mcpp_int get_divisor(mcpp_int input);
    mcpp_int get_divisor_multithread(mcpp_int input, mcpp_int down_limit, mcpp_int up_limit, mcpp_int &res);

    bool check_point(POINT p);
    bool check_ec();

    void getAllPoint();
    void generateStartPoint();
    POINT getRandomPoint();

public:
    vector<mcpp_int> divisors;
    /*Lenstra(mcpp_int number, mcpp_int limit) {
        input_N = number;
        ec.limit = limit;
    };*/
    Lenstra(mcpp_int number)
    {
        input_N = number;
    };

    vector<mcpp_int> run_sub_multi_thread();
    void run_same_thread(vector<mcpp_int> &box);
};

class _Lenstra {
public:
    mcpp_int input;
    vector<mcpp_int> res_box;
    mutex _mt;
    bool mute;
    _Lenstra(mcpp_int _input) {
        input = _input;
        mute = false;
    };

    void mute_out(bool _mute)
    {
        mute = _mute;
    }

    void print_divisors()
    {
        cout << endl;
        for (int i = 0; i < res_box.size(); i++)
        {
            cout << "[" << i << "] divisor = " << res_box[i] << endl;
        }
    }

    //Specify explicitly in which mode to factorize the number
    // ST - run in same thread
    // SMT - run in multithreading for a single elliptic curve
    // MT - run in multithread
    void run(int flag)
    {
        if (flag == ST)
        {
            run_same_thread();
            mute == false ? cout << "\nRun in same_thread mode.\n" << endl : cout;
        }
        else if (flag == SMT)
        {
            run_sub_multithread();
            mute == false ? cout << "\nRun in sub_multithread mode.\n" << endl : cout;
        }
        else if (flag == MT)
        {
            run_multi_thread();
            mute == false ? cout << "\nRun in multithread mode.\n" << endl : cout;
        }
    }

    //Automatically determines the operating mode depending
    // on the size of the factorized number
    void run()
    {
        if (input < 100000)
        {
            run_same_thread();
            mute == false ? cout << "\nRun in same_thread mode.\n" << endl : cout;
        }
        else if (input > 100000 and input < 10000000000)
        {
            run_sub_multithread();
            mute == false ? cout << "\nRun in sub_multithread mode.\n" << endl : cout;
        }
        else
        {
            run_multi_thread();
            mute == false ? cout << "\nRun in multithread mode.\n" << endl : cout;
        }
    }

    //—оздает несколько экземпл€ров класса Ћенстры в отдельных потоках
    vector<thread> ths;
    void run_multi_thread()
    {
        int counter = 0;
        while (true)
        {
            if (input % 2 == mcpp_int(0))
            {
                res_box.push_back(2);
                input = input / 2;
                counter++;
            }
            else if (input % 3 == mcpp_int(0))
            {
                res_box.push_back(3);
                input = input / 3;
                counter++;
            }
            else
                break;
        }
        if (counter!= 0)
            mute == false ? cout << "\nThe input number does not meet the requirements. The input of the algorithm will be submitted: " << input << endl : cout;

        vector<mcpp_int> res_box;
        for(int i=0; i<MAIN_THREADS_MAX; i++)
        {
            ths.push_back(std::thread([&]() {
                Lenstra{ input }.run_same_thread(std::ref(res_box));
                }));
        }

        for(int i=0; i< MAIN_THREADS_MAX; i++)
            ths[i].join();

        if (!mute)
        {
            cout << endl;
            for (int i = 0; i < res_box.size(); i++)
            {
                cout << "[" << i << "] divisor = " << res_box[i] << endl;
            }
        }
        ths.clear();
    }

    void run_same_thread()
    {
        int counter = 0;
        while (true)
        {
            if (input % 2 == mcpp_int(0))
            {
                res_box.push_back(2);
                input = input / 2;
                counter++;
            }
            else if (input % 3 == mcpp_int(0))
            {
                res_box.push_back(3);
                input = input / 3;
                counter++;
            }
            else
                break;
        }
        if (counter != 0)
            mute == false ? cout << "\nThe input number does not meet the requirements. The input of the algorithm will be submitted: " << input << endl : cout;

        Lenstra l(input);
        l.run_same_thread(std::ref(res_box));

        if (!mute)
        {
            for (int i = 0; i < res_box.size(); i++)
            {
                std::cout << "[" << i << "] divisor = " << res_box[i] << endl;
            }
            cout << "\nBAD EC SIZE = " << badEC.size();
        }
    };

    void run_sub_multithread()
    {
        int counter = 0;
        while (true)
        {
            if (input % 2 == mcpp_int(0))
            {
                res_box.push_back(2);
                input = input / 2;
                counter++;
            }
            else if (input % 3 == mcpp_int(0))
            {
                res_box.push_back(3);
                input = input / 3;
                counter++;
            }
            else
                break;
        }
        if (counter != 0)
            mute == false ? cout << "\nThe input number does not meet the requirements. The input of the algorithm will be submitted: " << input << endl : cout;

        Lenstra l(input);
        res_box = l.run_sub_multi_thread();

        if (!mute)
        {
            for (int i = 0; i < res_box.size(); i++)
            {
                std::cout << "[" << i << "] divisor = " << res_box[i] << endl;
            }
            cout << "\nBAD EC SIZE = " << badEC.size();
        }
    }
};