#pragma once
#ifndef H_VECTOR2
#define H_VECTOR2

#include "numeric.h"

#include <iostream>
#include <cmath>
#include <type_traits>

namespace dt {

	template<typename T>
	struct Vector2
	{
		using Type = T;
		Vector2() = default;
		//Vector2(const Vector2<T> &v) = default;
		//Vector2(Vector2<T>&&) = default;
		Vector2(const T vx, const T vy, const int idx);
		Vector2(const T vx, const T vy);
		Vector2(const Vector2<T> & pt);


		T dist2(const Vector2<T> &v) const;
		T dist(const Vector2<T> &v) const;
		T norm2() const;

		Vector2 &operator=(const Vector2<T>&) = default;
		Vector2 &operator=(Vector2&&) = default;
		//bool operator ==(const Vector2<T> &v) const;
		template<typename U>
		friend std::ostream &operator <<(std::ostream &str, const Vector2<U> &v);
		void setPoint(const T vx, const T vy, const int index);
		void setPoint(const T vx, const T vy);

		void setIndex(const int index);
		
		T x;
		T y;
		int index;
		int sub_idx;//用于subdivision的新增索引，工具索引

		//static_assert(std::is_floating_point<Vector2<T>::Type>::value,
		//	"Type must be floating-point");
	};

	template<typename T>
	bool almost_equal(const Vector2<T> &v1, const Vector2<T> &v2)
	{
		return almost_equal(v1.x, v2.x) && almost_equal(v1.y, v2.y);
	}

} // namespace dt

#endif
