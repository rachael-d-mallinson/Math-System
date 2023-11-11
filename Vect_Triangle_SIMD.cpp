//----------------------------------------------------------------------------
// Copyright 2021, Ed Keenan, all rights reserved.
//----------------------------------------------------------------------------

#include "Vect_TRIANGLE_SIMD.h"


Vect_TRIANGLE_SIMD::Vect_TRIANGLE_SIMD(const float tx, const float ty, const float tz, const float tw)
	:x(tx),y(ty),z(tz),w(tw)
{
}

//sqrtf((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
float Vect_TRIANGLE_SIMD::Length() const
{
	Vect_TRIANGLE_SIMD a;
	//a._m = _mm_mul_ps(this->_m,this->_m);
	//a._m = _mm_hadd_ps(a._m, a._m);
	//a._m = _mm_hadd_ps(a._m, a._m);
	//a._m = _mm_sqrt_ss(a._m);

	//nested the equations for efficiency
	a._m = _mm_sqrt_ss(_mm_hadd_ps(_mm_hadd_ps(_mm_mul_ps(this->_m, this->_m), _mm_mul_ps(this->_m, this->_m)), _mm_hadd_ps(_mm_mul_ps(this->_m, this->_m), _mm_mul_ps(this->_m, this->_m))));
	return a.x;
}

//Vect_TRIANGLE_SIMD(   ((this->y* R.z) - (R.y * this->z)),
//						((this->z* R.x) - (R.z * this->x)),
//						((this->x* R.y) - (R.x * this->y)));
Vect_TRIANGLE_SIMD Vect_TRIANGLE_SIMD::Cross(const Vect_TRIANGLE_SIMD &R) const
{
	//got down to 2 temp variables
	Vect_TRIANGLE_SIMD a;

	a._m = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(this->_m, this->_m, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(R._m, R._m, _MM_SHUFFLE(3, 1, 0, 2))), _mm_mul_ps(_mm_shuffle_ps(R._m, R._m, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(this->_m, this->_m, _MM_SHUFFLE(3, 1, 0, 2))));
	return a;
}

float Vect_TRIANGLE_SIMD::Area(const Vect_TRIANGLE_SIMD &T, const Vect_TRIANGLE_SIMD &R) 
{
	//cross then find the new area 
	return (0.5f * (T.Cross(R)).Length());

}


// ---  End of File ---

