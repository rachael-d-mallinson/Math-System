//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_Col_SIMD.h"
#include "Matrix_Col_SIMD.h"

Matrix_Col_SIMD::Matrix_Col_SIMD(const Vect_Col_SIMD &tV0,
								 const Vect_Col_SIMD &tV1,
								 const Vect_Col_SIMD &tV2,
								 const Vect_Col_SIMD &tV3)
	: v0(tV0), v1(tV1), v2(tV2), v3(tV3)
{
}

Vect_Col_SIMD Matrix_Col_SIMD::operator * (const Vect_Col_SIMD &vb)
{
	Vect_Col_SIMD a;
	//a._m = _mm_mul_ps(this->v0._m, vb._m);
	//a._m = _mm_mul_ps(this->v1._m, a._m);
	//a._m = _mm_mul_ps(this->v2._m, a._m);
	//a._m = _mm_mul_ps(this->v3._m, a._m);
	a = vb;
	//nested everything for efficiency
	a._m = _mm_mul_ps(this->v3._m, _mm_mul_ps(this->v3._m, _mm_mul_ps(this->v1._m, _mm_mul_ps(this->v0._m, vb._m))));
	return a;
}

Matrix_Col_SIMD Matrix_Col_SIMD::operator * (const Matrix_Col_SIMD &mb)
{
	Vect_Col_SIMD colA;
	Vect_Col_SIMD colHelp;
	
	Matrix_Col_SIMD newMat;
	newMat.m0 = this->m0;
	newMat.m1 = this->m4;
	newMat.m2 = this->m8;
	newMat.m3 = this->m12;

	newMat.m4 = this->m1;
	newMat.m5 = this->m5;
	newMat.m6 = this->m9;
	newMat.m7 = this->m13;

	newMat.m8 = this->m2;
	newMat.m9 = this->m6;
	newMat.m10 = this->m10;
	newMat.m11 = this->m14;

	newMat.m12 = this->m3;
	newMat.m13 = this->m7;
	newMat.m14 = this->m11;
	newMat.m15 = this->m15;

	Matrix_Col_SIMD newMat2;
	newMat2.m0 = mb.m0;
	newMat2.m1 = mb.m4;
	newMat2.m2 = mb.m8;
	newMat2.m3 = mb.m12;
			
	newMat2.m4 = mb.m1;
	newMat2.m5 = mb.m5;
	newMat2.m6 = mb.m9;
	newMat2.m7 = mb.m13;

	 newMat2.m8 = mb.m2;
	 newMat2.m9 = mb.m6;
	newMat2.m10 = mb.m10;
	newMat2.m11 = mb.m14;
				 
	newMat2.m12 = mb.m3;
	newMat2.m13 = mb.m7;
	newMat2.m14 = mb.m11;
	newMat2.m15 = mb.m15;


	colA._m = _mm_add_ps(_mm_mul_ps(newMat.v0._m, newMat2.v0._m), _mm_mul_ps(newMat.v0._m, newMat2.v1._m));
	colHelp._m = _mm_add_ps(_mm_mul_ps(newMat.v0._m, newMat2.v2._m), _mm_mul_ps(newMat.v0._m, newMat2.v3._m));
	colA._m = _mm_add_ps(colA._m, colHelp._m);

	Vect_Col_SIMD colB;
	colB._m = _mm_add_ps(_mm_mul_ps(newMat.v1._m, newMat2.v0._m), _mm_mul_ps(newMat.v1._m, newMat2.v1._m));
	colHelp._m = _mm_add_ps(_mm_mul_ps(newMat.v1._m, newMat2.v2._m), _mm_mul_ps(newMat.v1._m, newMat2.v3._m));
	colB._m = _mm_add_ps(colB._m, colHelp._m);

	Vect_Col_SIMD colC;
	colC._m = _mm_add_ps(_mm_mul_ps(newMat.v2._m, newMat2.v0._m), _mm_mul_ps(newMat.v2._m, newMat2.v1._m));
	colHelp._m = _mm_add_ps(_mm_mul_ps(newMat.v2._m, newMat2.v2._m), _mm_mul_ps(newMat.v2._m, newMat2.v3._m));
	colC._m = _mm_add_ps(colC._m, colHelp._m);

	Vect_Col_SIMD colD;
	colD._m = _mm_add_ps(_mm_mul_ps(newMat.v3._m, newMat2.v0._m), _mm_mul_ps(newMat.v3._m, newMat2.v1._m));
	colHelp._m = _mm_add_ps(_mm_mul_ps(newMat.v3._m, newMat2.v2._m), _mm_mul_ps(newMat.v3._m, newMat2.v3._m));
	colD._m = _mm_add_ps(colD._m, colHelp._m);

	//Trace::out("%d, %d, %d, %d\n", colA, colB, colC, colD);
	return Matrix_Col_SIMD(colA, colB, colC, colD);
}



// ---  End of File ---------------
