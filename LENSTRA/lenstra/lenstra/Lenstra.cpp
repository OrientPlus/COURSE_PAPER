#include "Lenstra.h"

vector<mcpp_int> Lenstra::run_sub_multi_thread()
{
    main_func(SMT);

    return divisors;
}

void Lenstra::run_same_thread(vector<mcpp_int> &box)
{
    main_func(ST);
    if (box.empty())
        box = divisors;
    return;
}

void Lenstra::main_func(int FLAG)
{
    mcpp_int input = input_N;
    vector<mcpp_int> res;

    mcpp_int int_divisor = input;
    while (true)
    {
        //Проверяем перед каждой итерацией не разложил ли другой поток уже число
        //Если KEEP_GOING == false - прекращаем работу до того как разложим число
        {
            std::lock_guard<mutex> lock(mt);
            if (!KEEP_GOING)
                break;
            iteration++;
        }

        if (miller_rabin_test(input, 25))
        {
            std::lock_guard<mutex> lock(mt);
            res.push_back(input);
            if (KEEP_GOING == true)
            {
                KEEP_GOING = false;
                divisors = res;
            }
            break;
        }

        int_divisor = intermediate_divisor(input, FLAG);

        if (int_divisor > 0)
        {
            res.push_back(int_divisor);
            input = input / int_divisor;
        }
        else if (int_divisor == -1)
            break;
    };
}

mcpp_int Lenstra::intermediate_divisor(mcpp_int input, int FLAG)
{
    mcpp_int newInput, divisor;
    init_b_value(input);

    divisor = generate_ec_coefficient(input);

    if (divisor != 1)
        return divisor;

    P.x = ec.x;
    P.y = ec.y;

    FINDED_DIV = false;
    vector<std::thread> threads;
    vector<mcpp_int> box;
    box.resize(DIVISOR_THREADS_MAX);

    mcpp_int step = input / DIVISOR_THREADS_MAX,
        divis;

    if (FLAG == SMT)
    {
        for (int i = 0; i < DIVISOR_THREADS_MAX; i++)
        {
            mcpp_int u, d;
            if (input % DIVISOR_THREADS_MAX != 0 and i == DIVISOR_THREADS_MAX - 1)
            {
                u = step * i;
                d = input;
            }
            else
            {
                u = step * i;
                d = u + step;
            }
            threads.push_back(std::thread([&]() {
                get_divisor_multithread(input, u, d, std::ref(divis));
                }));
        }

        for (auto& t : threads)
            t.join();

        divisor = divis;
    }
    else if (FLAG == ST) {
        divisor = get_divisor(input);
    }


    if (divisor == -1)
        return -1;

    if (divisor > 0)
    {
        return divisor;
    }
    else
    {
        lock_guard<mutex> lock(mt);
        //Сохраняем неудачные ЭК
        badEC.push_back(ec);
        return 0;
    }
}

// B = N/4
void Lenstra::init_b_value(mcpp_int input)
{
    B = input / 4;
} 

// Generates a random curve in Weierstrass normal form
mcpp_int Lenstra::generate_ec_coefficient(mcpp_int input)
{
    {
        std::lock_guard<mutex> lock(mt);
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        seed %= 1000;


        boost::random::mt19937 gen(seed);
        boost::random::uniform_int_distribution<mcpp_int> dist(1, input);

        ec.limit = input;

        ec.x = dist(gen);
        ec.y = dist(gen);
        ec.a = dist(gen);

        ec.b = pow(ec.y, 2) - pow(ec.x, 3) - ec.a * ec.x;
        if (ec.b < 0 or ec.b > ec.limit)
        {
            ec.b = ec.b % ec.limit;
            if (ec.b < 0)
                ec.b = ec.limit + ec.b;
        }
    }

    mcpp_int check = pow(ec.a, 3) * 4 + pow(ec.b, 2) * 27;
    check %= ec.limit;
    if (check == 0)
    {   
        std::lock_guard<mutex> lock(mt);
        badEC.push_back(ec);
        generate_ec_coefficient(input);
    }

    if (check_ec())
        generate_ec_coefficient(input);

    mcpp_int g = gcd(input, (pow(ec.a, 3) * 4 + pow(ec.b, 2) * 27));

    if (g == input)
        generate_ec_coefficient(input);
    else if (g > 1 and g < input)
        return g;
    else if (g == 1)
        return 1;
}

// The main function of the Lenstra algorithm. Finds the divisor of a number
mcpp_int Lenstra::get_divisor(mcpp_int input)
{
    bool undefined_point_flag = false, fl = false;
    POINT Q1, Q2, Q, oldPoint = P;
    mcpp_int oldCoef = 0;

    Q = P;
    Q1 = P;

    //Для каждого простого i<B найдем наибольший коэф. r такой что i^r < B
    for (mcpp_int i = 2; i < B; i++)
    {
        if (!miller_rabin_test(i, 15))
            continue;
        mcpp_int tmp = i;
        unsigned long long r = 0, degree = 2;
        while (true)
        {
            if (!KEEP_GOING)
                return -1;

            tmp = pow(i, degree);
            if (tmp > B)
            {
                r = degree - 1;
                break;
            }
            degree++;
        }

        mcpp_int k = pow(i, r);

        //Умножаем методом удвоения точки, чтобы сократить время затрачиваемое на умножение
        //Если получим точку, которая определена, переходим к следующей итерации
        //Иначе идем дальше и последовательно складываем точку Q + Q + Q + Q +...+ Q и ищем итерацию на которой получается неопределенная точка
        // СИЛЬНО УВЕЛИЧИВАЕТ ВРЕМЯ ВЫПОЛНЕНИЯ АЛГОРИТМА
        /*Q = P * k;
        if (Q != POINT{ -1, -1 })
            continue;*/

        //Q = P;
        //Q *= k;
        //Q = Q % ec.limit;
        //auto it = find(allPoint.begin(), allPoint.end(), Q);
        //if (it != allPoint.end())
            //return 0;

        if (k > oldCoef)
            Q = oldPoint;
        else
        {   
            oldCoef = 0;
            Q = P;
        }
        Q1 = P;
        Q2 = P;


        //Последовательно умножаем точку Q на k = i^r
        for (mcpp_int j = 0; j < k - oldCoef; j++)
        {
            if (!KEEP_GOING)
                return -1;
            Q = ec.plus(Q1, Q2);
            if (Q == POINT{ -1,-1 })
                fl = true;
            else
            {
                Q = Q % ec.limit;
                Q1 = Q;
                continue;
            }

            // Флаг неопределенной точки на кривой
            if (fl == true)
            {
                mcpp_int _x1, _x2, _y1, _y2;
                _x1 = Q1.x % ec.limit;
                _x2 = Q2.x % ec.limit;
                _y1 = Q1.y % ec.limit;
                _y2 = Q2.y % ec.limit;
                _y2 = ec.limit - _y2;

                if (Q == Q1)
                {
                    mcpp_int d = _y1 + _y2;
                    d %= ec.limit;
                    d = gcd(d, input);
                    
                    return d;
                }
                else if (_x1 != _x2)
                {
                    mcpp_int d = _x1 - _x2;
                    d %= ec.limit;
                    d = gcd(d, input);
                    if (d > 1 and d < input)
                    {
                        return d;
                    }
                }
                else if (_y1 != _y2)
                {
                    mcpp_int d = _y1 + _y2;
                    d %= ec.limit;
                    d = gcd(Q.y + Q1.y, input);
                    if (d > 1 and d < input)
                    {
                        return d;
                    }
                }
            }
        }
        if (i == B - 1)
        {
            std::lock_guard<mutex> lock(mt);
            B = input;
            //cout << "\nNo suitable point was found! The border has been increased to the maximum possible!";
        }
        oldPoint = Q;
        oldCoef = k;
    }
    return 0;
}

// Finds the divisor by implementing a multithreaded version of the main loop
mcpp_int Lenstra::get_divisor_multithread(mcpp_int input, mcpp_int down_limit, mcpp_int up_limit, mcpp_int &res)
{
    bool undefined_point_flag = false, fl = false;
    POINT Q1, Q2, Q, oldPoint = P;
    mcpp_int oldCoef = 0;
    //Для каждого i<B найдем наибольший коэф. r такой что i^r < B
    for (mcpp_int i = down_limit; i < up_limit; i++)
    {
        if (!miller_rabin_test(i, 15))
        {
            continue;
        }
        mcpp_int tmp = i;
        unsigned long long r, degree = 2;
        while (true)
        {
            if (!KEEP_GOING)
                return 0;
            if (FINDED_DIV)
                return 0;

            tmp = pow(i, degree);
            if (tmp > B)
            {
                r = degree - 1;
                break;
            }
            degree++;
        }

        mcpp_int k = pow(i, r);
        POINT Q1, Q2, Q;

        if (k > oldCoef)
            Q = oldPoint;
        else
        {
            oldCoef = 0;
            Q = P;
        }
        Q1 = P;
        Q2 = P;

        //Последовательно умножаем точку Q на p^r
        for (mcpp_int j = 0; j < k - oldCoef; j++)
        {
            if (!KEEP_GOING)
                return -1;
            {
                if (FINDED_DIV)
                    return 0;
            }
            Q = ec.plus(Q1, Q2);
            if (Q == INF_POINT)
                fl = true;
            else
            {
                Q = Q % ec.limit;
                Q1 = Q;
                continue;
            }
            if (fl == true)
            {
                mcpp_int _x1, _x2, _y1, _y2;
                _x1 = Q1.x % ec.limit;
                _x2 = Q2.x % ec.limit;
                _y1 = Q1.y % ec.limit;
                _y2 = Q2.y % ec.limit;
                _y2 = ec.limit - _y2;

                if (Q == Q1)
                {
                    mcpp_int d = _y1 + _y2;
                    d %= ec.limit;
                    d = gcd(d, input);
                    if (d > 1 and d < input)
                    {
                        std::lock_guard<std::mutex> lock(mt);
                        if (FINDED_DIV == false)
                        {
                            FINDED_DIV = true;
                            res = d;
                            return d;
                        }
                        else
                            return d;
                    }
                }
                else if (_x1 != _x2)
                {
                    mcpp_int d = _x1 - _x2;
                    d %= ec.limit;
                    d = gcd(d, input);
                    if (d > 1 and d < input)
                    {
                        std::lock_guard<std::mutex> lock(mt);
                        if (FINDED_DIV == false)
                        {
                            FINDED_DIV = true;
                            res = d;
                            return d;
                        }
                        else
                            return d;
                    }
                }
                else if (_y1 != _y2)
                {
                    mcpp_int d = _y1 + _y2;
                    d %= ec.limit;
                    d = gcd(Q.y + Q1.y, input);
                    if (d > 1 and d < input)
                    {
                        std::lock_guard<std::mutex> lock(mt);
                        if (FINDED_DIV == false)
                        {
                            FINDED_DIV = true;
                            res = d;
                            return d;
                        }
                        else
                            return d;
                    }
                }
            }
        }
        if (i == B - 1)
        {
            std::lock_guard<mutex> lock(mt);
            B = input;
            //cout << "\nNo suitable point was found! The border has been increased to the maximum possible!";
        }
        oldPoint = Q;
        oldCoef = k;
    }
    return 0;
}

// Checks whether the specified curve has already been used
bool Lenstra::check_ec()
{
    std::lock_guard<mutex> lock(mt);

    for (auto _ec : badEC)
    {
        if (ec.a == _ec.a and ec.x == _ec.x and ec.y == _ec.y)
        {
            return true;
        }
    }

    return false;
}
