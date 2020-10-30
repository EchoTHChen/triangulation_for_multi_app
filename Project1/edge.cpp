#include "edge.h"

namespace dt {

	template<typename T>
	Edge<T>::Edge(const VertexType &v1, const VertexType &v2) 
		//:v(&v1), w(&v2)
	{
		//深拷贝构造函数；
		v = new VertexType(v1.x, v1.y, v1.index);
		w = new VertexType(v2.x, v2.y, v2.index);
	}
	template<typename T>
	void Edge<T>::setEdge(const VertexType &v1, const VertexType &v2) {
		v = new VertexType(v1.x, v1.y, v1.index);
		w = new VertexType(v2.x, v2.y, v2.index);
	}
	//template<typename T>
	//bool
	//	Edge<T>::operator ==(const Edge<T> &e) const
	//{
	//	return (*(this->v) == *e.v && *(this->w) == *e.w) ||
	//		(*(this->v) == *e.w && *(this->w) == *e.v);
	//}

	template<typename U>
	std::ostream&
		operator <<(std::ostream &str, const Edge<U> &e)
	{
		return str << "Edge " << *e.v << ", " << *e.w;
	}

	template struct Edge<float>;
	template struct Edge<double>;

} // namespace dt
