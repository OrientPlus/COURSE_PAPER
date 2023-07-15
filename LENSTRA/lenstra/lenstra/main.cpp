#include "Lenstra.h"
#include "ECDH.hpp"


// Запускает алогритм Ленстры itearations раз, с указанным флагом.
// На входе алгоритма всегда input_number
// Возвращает среднее кол-во времени в миллисекундах затраченное на факторизацию
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

	/*ECDH ecdhAlice(1842016301), ecdhBob;

	publicParameter ppA = ecdhAlice.gen_main_parameters();
	ecdhBob.set_main_parameters(ppA);

	if (ecdhAlice.gen_shared_secret(ecdhBob.public_key) == 1 or ecdhBob.gen_shared_secret(ecdhAlice.public_key) == 1)
	{
		cout << "\nIt is necessary to regenerate the main parameters of the algorithm. \n\
The resulting shared key is not defined on this elliptic curve.";
	}
	else
		cout << "\nShared key successfully generated!";*/

	/*ELLEPTIC_CURVE ec;
	ec.limit = mcpp_int{ "37654765467534762343413417"};
	ec.a = 7;
	ec.b = 9;
	
	ec.gen_base_point();
	//ec.basePoint.x = 142;
	//ec.basePoint.y = 213;
	
	cout << "\nbp.x = " << ec.basePoint.x << "\nbp.y = " << ec.basePoint.y;


	POINT res, res2, res3;
	//mcpp_int c("4534534534534534777");
	for(int i=2; i<200; i++)
	{
		res = ec.stupidmult(ec.basePoint, i);
		res2 = ec.fastmult(ec.basePoint, i);
		res3 = ec.mult(ec.basePoint, i);

		if (res != res2 or res != res3)
		{
			cout << "\n-------------------------\nI = " << i;
			cout << "\nStupid mult:\nres.x = " << res.x << "\nres.y = " << res.y;
			cout << "\nFastmult:\nres2.x = " << res2.x << "\nres2.y = " << res2.y;
			cout << "\nSimple mult:\nres3.x = " << res3.x << "\nres3.y = " << res3.y;
			cout << "\n-------------------------\n";
			continue;
		}
		cout << "\nI->" << i;
	}*/
 	cout << endl;
	system("pause");
	return 0;
}