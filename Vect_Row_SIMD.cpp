//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_Row_SIMD.h"
#include "Matrix_Row_SIMD.h"

Vect_Row_SIMD::Vect_Row_SIMD(const float tx, const float ty, const float tz, const float tw)
	: x(tx), y(ty), z(tz), w(tw)
{
}

Vect_Row_SIMD Vect_Row_SIMD::operator * ( const Matrix_Row_SIMD &ma)
{
	//create 2 temps to load with x and y
	Vect_Row_SIMD temp1;
	Vect_Row_SIMD temp2;
	//cerate 2 temps to hold my ending multiplications
	Vect_Row_SIMD vOut1;
	Vect_Row_SIMD vOut2;

	//create temp to hold the first row of matrix
	Vect_Row_SIMD vectFromMat = Vect_Row_SIMD(ma.m0, ma.m1, ma.m2, ma.m3);
	//load temps with x and y, multiply with first row of matrix
	temp1._m = _mm_set1_ps(this->x);
	temp2._m = _mm_set1_ps(this->y);
	temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
	//change temp to hold the second row of matrix, multiply with y's
	vectFromMat = Vect_Row_SIMD(ma.m4, ma.m5, ma.m6, ma.m7);
	temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

	vOut1._m = _mm_add_ps(temp1._m, temp2._m);
	// load temps with z and w, multiply with next rows of matrix
	temp1._m = _mm_set1_ps(this->z);
	temp2._m = _mm_set1_ps(this->w);
	vectFromMat = Vect_Row_SIMD(ma.m8, ma.m9, ma.m10, ma.m11);
	temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
	vectFromMat = Vect_Row_SIMD(ma.m12, ma.m13, ma.m14, ma.m15);
	temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

	vOut2._m = _mm_add_ps(temp1._m, temp2._m);

	vOut1._m = _mm_add_ps(vOut1._m, vOut2._m);
	return vOut1;
};

//start of proxies
// vector * matrix
struct vM
{
	const Vect_Row_SIMD v;
	const Matrix_Row_SIMD m;

	vM(const Vect_Row_SIMD& _v, const Matrix_Row_SIMD& _m)
		:v(_v), m(_m) {};

	operator Vect_Row_SIMD()
	{
		//create 2 temps to load with x and y
		Vect_Row_SIMD temp1;
		Vect_Row_SIMD temp2;
		//cerate 2 temps to hold my ending multiplications
		Vect_Row_SIMD vOut1;
		Vect_Row_SIMD vOut2;

		//create temp to hold the first row of matrix
		Vect_Row_SIMD vectFromMat = Vect_Row_SIMD(m.m0, m.m1, m.m2, m.m3);
		//load temps with x and y, multiply with first row of matrix
		temp1._m = _mm_set1_ps(v.x);
		temp2._m = _mm_set1_ps(v.y);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		//change temp to hold the second row of matrix, multiply with y's
		vectFromMat = Vect_Row_SIMD(m.m4, m.m5, m.m6, m.m7);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vOut1._m = _mm_add_ps(temp1._m, temp2._m);
		// load temps with z and w, multiply with next rows of matrix
		temp1._m = _mm_set1_ps(v.z);
		temp2._m = _mm_set1_ps(v.w);
		vectFromMat = Vect_Row_SIMD(m.m8, m.m9, m.m10, m.m11);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		vectFromMat = Vect_Row_SIMD(m.m12, m.m13, m.m14, m.m15);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vOut2._m = _mm_add_ps(temp1._m, temp2._m);

		vOut1._m = _mm_add_ps(vOut1._m, vOut2._m);
		return vOut1;
	}
};

inline vM operator * (const Vect_Row_SIMD& v, const Matrix_Row_SIMD& m)
{
	return vM(v, m);
}

struct vTimesMatTimesMat
{
	const Vect_Row_SIMD v;
	const Matrix_Row_SIMD m; 
	const Matrix_Row_SIMD m1;

	vTimesMatTimesMat(const vM& _vM, const Matrix_Row_SIMD& _m1)
		:v(_vM.v), m(_vM.m), m1(_m1) {};

	operator Vect_Row_SIMD()
	{
		//create 2 temps to load with x and y
		Vect_Row_SIMD temp1;
		Vect_Row_SIMD temp2;
		//cerate 2 temps to hold my ending multiplications
		Vect_Row_SIMD vOut1;
		Vect_Row_SIMD vectResult;

		//create temp to hold the first row of matrix
		Vect_Row_SIMD vectFromMat = Vect_Row_SIMD(m.m0, m.m1, m.m2, m.m3);
		//load temps with x and y, multiply with first row of matrix
		temp1._m = _mm_set1_ps(v.x);
		temp2._m = _mm_set1_ps(v.y);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		//change temp to hold the second row of matrix, multiply with y's
		vectFromMat = Vect_Row_SIMD(m.m4, m.m5, m.m6, m.m7);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vOut1._m = _mm_add_ps(temp1._m, temp2._m);
		// load temps with z and w, multiply with next rows of matrix
		temp1._m = _mm_set1_ps(v.z);
		temp2._m = _mm_set1_ps(v.w);
		vectFromMat = Vect_Row_SIMD(m.m8, m.m9, m.m10, m.m11);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		vectFromMat = Vect_Row_SIMD(m.m12, m.m13, m.m14, m.m15);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vectResult._m = _mm_add_ps(vOut1._m, _mm_add_ps(temp1._m, temp2._m));
		
		//now new matrix
		//create temp to hold the first row of matrix
		vectFromMat = Vect_Row_SIMD(m1.m0, m1.m1, m1.m2, m1.m3);
		//load temps with x and y, multiply with first row of matrix
		temp1._m = _mm_set1_ps(vectResult.x);
		temp2._m = _mm_set1_ps(vectResult.y);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		//change temp to hold the second row of matrix, multiply with y's
		vectFromMat = Vect_Row_SIMD(m1.m4, m1.m5, m1.m6, m1.m7);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vOut1._m = _mm_add_ps(temp1._m, temp2._m);
		// load temps with z and w, multiply with next rows of matrix
		temp1._m = _mm_set1_ps(vectResult.z);
		temp2._m = _mm_set1_ps(vectResult.w);
		vectFromMat = Vect_Row_SIMD(m1.m8, m1.m9, m1.m10, m1.m11);
		temp1._m = _mm_mul_ps(temp1._m, vectFromMat._m);
		vectFromMat = Vect_Row_SIMD(m1.m12, m1.m13, m1.m14, m1.m15);
		temp2._m = _mm_mul_ps(temp2._m, vectFromMat._m);

		vectResult._m = _mm_add_ps(vOut1._m, _mm_add_ps(temp1._m, temp2._m));
		
		return vectResult;
	}
};

inline vTimesMatTimesMat operator *(const vM& _vM, const Matrix_Row_SIMD& _m)
{
	return vTimesMatTimesMat(_vM, _m);
}

// ---  End of File ---------------


