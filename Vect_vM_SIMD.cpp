//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_vM_SIMD.h"
#include "Matrix_vM_SIMD.h"

Vect_vM_SIMD::Vect_vM_SIMD(const float tx, const float ty, const float tz, const float tw)
	: x(tx), y(ty), z(tz), w(tw)
{
}

//------------------------------
// For No proxy:
//------------------------------

Vect_vM_SIMD Vect_vM_SIMD::operator * (const Matrix_vM_SIMD &ma) const
{
	//cerate 2 temps to hold my ending multiplications
	Vect_vM_SIMD vOut1;
	
	// load temps with z and w, multiply with next rows of matrix, and so on
	vOut1._m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->x), ma.v0._m), _mm_mul_ps(_mm_set1_ps(this->y), ma.v1._m)),
						  _mm_add_ps(_mm_mul_ps(_mm_set1_ps(this->z), ma.v2._m), _mm_mul_ps(_mm_set1_ps(this->w), ma.v3._m)));
	return vOut1;
};


// ---  End of File ---------------


