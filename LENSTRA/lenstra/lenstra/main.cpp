#include "Lenstra.h"
#include "ECDH.hpp"


// Runs the Lenstra algorithm 'iterations' times, with the specified flag.
// На входе алгоритма всегда input_number
// Returns the average amount of time in milliseconds spent on a single factorization
mcpp_int lenstra_test(mcpp_int input_number, int iterations, int mode)
{
	timer t(2);

	for (int i = 0; i < iterations; i++)
	{
		_Lenstra l(input_number);
		l.mute_out(true);
		l.run(mode);
	}

	return t.get_current_time() / iterations;
}


int main(int argc, char* argv[])
{
	//376416230105040466653667338727015651682345823271261		= 1521216602642502411343 * 2768810951 * 4469581842277 * 19994801
	//1237542715128736451287354231757							= 433 * 1878803 * 1521216602642502411343
	//2858066316694541457938462429								= 1878803 * 1521216602642502411343
	//12375427151287312375427									= 2768810951 * 4469581842277
	//10905898707717952193										= 10003559 * 108967 * 10004881 
	//128547635353453											= 59 * 108967 * 19994801
	//12854763													= 3 * 3 * 71 * 20117
	//12854765													= 5 * 7 * 11 * 173 * 193
	//21403														= 17 * 1259
	/*if (argc > 1)
	{
		std::stringstream ss(argv[1]);
		string str;
		ss >> str;


		mcpp_int n(str);
		_Lenstra l(n);
		
		cout << "N = " << n << endl;

		l.run();
	}
	else {
		string str;
		cout << "Enter the input number: ";
		cin >> str;
		mcpp_int n(str);
		_Lenstra l(n);

		l.run();
	}*/

	mcpp_int in("10905898707717952193");
	mcpp_int avg1 = 0, avg2 = 0, avg3 = 0;

	//avg1 = lenstra_test(in, 100, ST);
	//avg2 = lenstra_test(in, 100, SMT);
	avg3 = lenstra_test(in, 100, MT);


	cout << "\n------------------------------------------\nAverage time for single-threaded mode: " << avg1  << " millisec." << "\n------------------------------------------";
	cout << "\n------------------------------------------\nAverage time for sub-multithreaded mode: " << avg2 << " millisec." << "\n------------------------------------------";
	cout << "\n------------------------------------------\nAverage time for multi-threaded mode: " << avg3 << " millisec." << "\n------------------------------------------";

 	cout << endl;
	system("pause");
	return 0;
}