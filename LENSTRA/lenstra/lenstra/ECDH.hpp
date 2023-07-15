#pragma once

#include "openssl/sha.h"
#include "openssl/evp.h"
#include "openssl/bn.h"

#include "elleptic_curve.hpp"


using std::pair;

class publicParameter
{
public:
	publicParameter() {};
	publicParameter(const publicParameter& other)
	{
		ec = other.ec;
		public_key = other.public_key;
	}
	ELLEPTIC_CURVE ec;
	POINT public_key;
};


class ECDH 
{
public:
	ECDH() { 
		mcpp_int def_limit("52064339259162311797");
		ec.limit = def_limit;
		ec.orderBP = 0; 
		ec.orderCurve = 0; 
	};
	ECDH(mcpp_int limit) {
		if (limit % 2 == 0 or limit % 3 == 0)
		{
			cout << "\nIncorrect field size is set! The field size will be set by default: 52064339259162311797.\n";
			mcpp_int def_limit("52064339259162311797");
			ec.limit = def_limit;
		}
		else
			ec.limit = limit;
	};
	ECDH(publicParameter _pp)
	{
		if (_pp.ec.limit % 2 == 0 or _pp.ec.limit % 3 == 0)
		{
			cout << "\nIncorrect field size is set! The field size will be set by default: 52064339259162311797.\n\
Incorrect field size is set! The field size will be set by default. \n\
ATTENTION, public parameters will be changed! \n\
If you received these parameters from another user, \n\
send him the updated public parameters to get a private key shared with him!";
			mcpp_int def_limit("52064339259162311797");
			_pp.ec.limit = def_limit;
		}
		pp.ec = _pp.ec;
		pp.public_key = _pp.public_key;
	}

	POINT public_key;
	publicParameter pp;

	publicParameter gen_main_parameters()
	{
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		boost::random::mt19937 gen(seed);
		boost::random::uniform_int_distribution<mcpp_int> dist(0, ec.limit);

		// Генерируем случайные коэф-ты a,b ЭК
		ec.a = dist(gen);
		ec.b = dist(gen);

		// Генерируем случайную базовую точку
		ec.gen_base_point();
		cout << "\nBASE POINT:\nX = " << ec.basePoint.x << "\nY = " << ec.basePoint.y << endl;
		ec.find_order_BP();
		cout << "\bBASE POINT ORDER = " << ec.orderBP;

		pp.ec = ec;

		gen_key();
		pp.public_key = public_key;

		return pp;
	}
	int set_main_parameters(publicParameter _pp)
	{
		pp = _pp;
		ec = pp.ec;
		gen_key();
		pp.public_key = public_key;
		return 0;
	}
	
	int gen_shared_secret(POINT other_public_key)
	{
		POINT tmp = ec.mt_mult(other_public_key, secret_key);
		if (tmp == INF_POINT)
		{
			cout << "\nbad shared secret(\n";
			return 1;
		}
		shared_secret = tmp.x;
		return 0;
	} 


private:
	mcpp_int secret_key;
	mcpp_int shared_secret;
	ELLEPTIC_CURVE ec;


	mcpp_int get_hash(string data)
	{
		unsigned char hashBuf[SHA256_DIGEST_LENGTH];

		SHA256(reinterpret_cast<const unsigned char*>(data.c_str()), data.length(), hashBuf);
		std::stringstream ss;
		for (int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
			ss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(hashBuf[i]);
		}

		string s_hash = ss.str();
		s_hash.insert(0, "0x");
		mcpp_int hash(s_hash);

		return hash;
	}
	void gen_key()
	{
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		boost::random::mt19937 gen(seed);
		boost::random::uniform_int_distribution<mcpp_int> dist(1, ec.orderBP);
	genSK:
		cout << "\nNew iteration";
		secret_key = dist(gen);
		public_key = ec.mt_mult(ec.basePoint, secret_key);
		if (public_key == INF_POINT)
			goto genSK;
	}
	mcpp_int mpow(mcpp_int a, mcpp_int b)
	{
		mcpp_int res = a, tmp = a;
		while (b > 10000)
		{
			res = pow(res, 10000);
			b = b / 10000;
		}
		int last_pow = stoi(b.str());
		res = pow(res, last_pow);

		return res;
	}
};