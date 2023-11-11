//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_M_SIMD.h"
#include "Matrix_M_SIMD.h"

Matrix_M_SIMD::Matrix_M_SIMD(const Vect_M_SIMD &tV0,
							 const Vect_M_SIMD &tV1,
							 const Vect_M_SIMD &tV2,
							 const Vect_M_SIMD &tV3)
	: v0(tV0), v1(tV1), v2(tV2), v3(tV3)
{
}

Matrix_M_SIMD Matrix_M_SIMD::operator * (const Matrix_M_SIMD &mb) const
{		
	//temps to hold the ending new matrix lines
		Matrix_M_SIMD MatOut;
	
		//first row (or vector v0)
		MatOut.v0._m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m0), mb.v0._m), _mm_mul_ps(_mm_set1_ps(this->m1), mb.v1._m)),
								  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m2), mb.v2._m), _mm_mul_ps(_mm_set1_ps(this->m3), mb.v3._m)));

		//second row (or vector v1)
		MatOut.v1._m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m4), mb.v0._m), _mm_mul_ps(_mm_set1_ps(this->m5), mb.v1._m)),
								  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m6), mb.v2._m), _mm_mul_ps(_mm_set1_ps(this->m7), mb.v3._m)));

		//third row (or vector v2)
		MatOut.v2._m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m8), mb.v0._m), _mm_mul_ps(_mm_set1_ps(this->m9), mb.v1._m)),
								  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m10), mb.v2._m), _mm_mul_ps(_mm_set1_ps(this->m11), mb.v3._m)));

		//fourth row (or vector v3)
		MatOut.v3._m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m12), mb.v0._m), _mm_mul_ps(_mm_set1_ps(this->m13), mb.v1._m)),
								  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->m14), mb.v2._m), _mm_mul_ps(_mm_set1_ps(this->m15), mb.v3._m)));

	return MatOut;
}


// ---  End of File ---------------
