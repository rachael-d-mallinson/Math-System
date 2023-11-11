//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_Mv_SIMD.h"
#include "Matrix_Mv_SIMD.h"

Matrix_Mv_SIMD::Matrix_Mv_SIMD(const Vect_Mv_SIMD &tV0,
							   const Vect_Mv_SIMD &tV1,
							   const Vect_Mv_SIMD &tV2,
							   const Vect_Mv_SIMD &tV3)
	: v0(tV0), v1(tV1), v2(tV2), v3(tV3)
{
}

Vect_Mv_SIMD Matrix_Mv_SIMD::operator * (const Vect_Mv_SIMD &v) const
{
	//create temps 
	Vect_Mv_SIMD A;
	
	//nested for efficiency
	A._m = _mm_hadd_ps(_mm_hadd_ps(_mm_mul_ps(this->v0._m, v._m), _mm_mul_ps(this->v1._m, v._m)), _mm_hadd_ps(_mm_mul_ps(this->v2._m, v._m), _mm_mul_ps(this->v3._m, v._m)));
	return A;
}
// ---  End of File ---------------
