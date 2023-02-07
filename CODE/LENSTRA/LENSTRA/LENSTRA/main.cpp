#include <climits>
#include <openssl/bn.h>
#include <stdio.h>

int main() {
	BN_CTX* ctx = BN_CTX_new(); // �������� ��� ����������
	BIGNUM* mybignum = nullptr;   // �����, ������� ����� ��������
	BIGNUM* mul = BN_new();     // ��������� ���������

	BN_dec2bn(&mybignum, "18446744073709551615"); // 2^64-1
	BN_mul(mul, mybignum, mybignum, ctx);

	// ����������� ���������
	char* dec = BN_bn2dec(mul);
	if (dec) {
		printf("%s\n", dec);
		OPENSSL_free(dec);
	}

	// ��������� ���������� ������
	BN_free(mybignum);
	BN_free(mul);
}