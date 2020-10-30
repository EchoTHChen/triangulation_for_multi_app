#pragma once
#ifndef H_INTERSECTION
#define H_INTERSECTION

namespace dt{

	/*template<typename T>
	bool pt_cmp(const dt::Vector2<T> &p1, const dt::Vector2<T>& p2);
*/
	bool pt_cmp(const dt::Vector2<double> &p1, const dt::Vector2<double>& p2);


	/*template<typename T>
	struct Intersection{
		using Type = T;
		Intersection() = default;
*/
	template<typename T>
	T min_value(T a, T b);
	template double min_value(double a, double b);

	template<typename T>
	T max_value(T a, T b);
	template double max_value(double a, double b);
	
	/*template<typename T>
	bool almost_equal(T a, T b, T eps = 1e-5);
	template bool almost_equal(double a, double b, double eps = 1e-5);

*/
	template<typename T>
	bool PointInPolygon(dt::Vector2<T> p, std::vector<dt::Vector2<T>> &  ptPolygon);
	template bool PointInPolygon(dt::Vector2<double> p, std::vector<dt::Vector2<double>> &  ptPolygon);


		//void print_edges(const std::vector<dt::Edge<double>> edges);
		//void print_vec(std::vector<double> vec);
	template<typename T>
	bool isRectCross(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> q1, dt::Vector2<T> q2);
	template bool isRectCross(dt::Vector2<double> p1, dt::Vector2<double> p2, dt::Vector2<double> q1, dt::Vector2<double> q2);
	
	template<typename T>
	T cross_multi(dt::Vector2<T> pp0, dt::Vector2<T>pp1, dt::Vector2<T>pp2);//pp0为公用点，叉积计算
	template double cross_multi(dt::Vector2<double> pp0, dt::Vector2<double>pp1, dt::Vector2<double>pp2);//pp0为公用点，叉积计算
	
	template<typename T>
	bool isLineSegmentCross(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> q1, dt::Vector2<T> q2);
	template bool isLineSegmentCross(dt::Vector2<double> p1, dt::Vector2<double> p2, dt::Vector2<double> q1, dt::Vector2<double> q2);

	template<typename T>
	T cross_multi(dt::Vector2<T> vec1, dt::Vector2<T>vec2);
	template double cross_multi(dt::Vector2<double> vec1, dt::Vector2<double>vec2);

	template<typename T>
	bool GetCrossPoint(dt::Edge<T> cur_edge1, dt::Edge<T>cur_edge2, dt::Vector2<T>& cross_pt);
	template bool GetCrossPoint(dt::Edge<double> cur_edge1, dt::Edge<double>cur_edge2, dt::Vector2<double>& cross_pt);

	////bool pt_cmp(const dt::Vector2<T> &p1, const dt::Vector2<T>& p2);

	//3.去重复
	template<typename T>
	void del_repeat(std::vector<dt::Vector2<T>> &ans_tmp, std::vector<dt::Vector2<T>> &ans);
	template void del_repeat(std::vector<dt::Vector2<double>> &ans_tmp, std::vector<dt::Vector2<double>> &ans);

	////4.逆时针或顺时针排序
	////若点A大于B,即点A在点B的顺时针方向，返回true,否则返回false;//
	template<typename T>
	bool CSortCmp(const dt::Vector2<T> &a, const dt::Vector2<T> &b, const dt::Vector2<T> & center);
	template bool CSortCmp(const dt::Vector2<double> &a, const dt::Vector2<double> &b, const dt::Vector2<double> & center);

	template<typename T>
	void ClockwiseSortPoints(std::vector<dt::Vector2<T>> &ans_tmp);
	template void ClockwiseSortPoints(std::vector<dt::Vector2<double>> &ans_tmp);
	
	template<typename T>
	void intersect_triangles(dt::Triangle<T> t1, dt::Triangle<T> t2, std::vector<dt::Vector2<T>> &ans_tmp, int total);
	template void intersect_triangles(dt::Triangle<double> t1, dt::Triangle<double> t2, std::vector<dt::Vector2<double>> &ans_tmp, int total);
	
	template<typename T>
	void PolygonClip(const dt::Delaunay<T> &fg, const dt::Delaunay<T>& bg, std::vector<dt::Vector2<T>> & ans_tmp, int total);
	template void PolygonClip(const dt::Delaunay<double> &fg, const dt::Delaunay<double>& bg, std::vector<dt::Vector2<double>> & ans_tmp, int total);

	template<typename T>
	T g(std::vector<dt::Vector2<T>> pts, int a, int b, int c);
	template double g(std::vector<dt::Vector2<double>> pts, int a, int b, int c);

	template<typename T>
	bool judge_convex_concave(std::vector<dt::Vector2<T>> pts);
	template bool judge_convex_concave(std::vector<dt::Vector2<double>> pts);
	
	template<typename T>
	std::vector<dt::Vector2<T>> regulation_sort_pt_map(std::vector<dt::Vector2<T>> tmp_pts, std::map<int, int> &pt_map);
	template std::vector<dt::Vector2<double>> regulation_sort_pt_map(std::vector<dt::Vector2<double>> tmp_pts, std::map<int, int> &pt_map);

	//TODO:三角面片对点进行规整；？？？把多个离散 的多边形分开。
	template<typename T>
	void regulation_triangles(std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>> &new_pts);//修改index；points;
	template void regulation_triangles(std::vector<dt::Triangle<double>> &inter_triangles_ans, std::vector<dt::Vector2<double>> &new_pts);//修改index；points;

	template<typename T>
	void index_intersection(dt::Delaunay<T>& fg, dt::Delaunay<T>& bg, std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>>& new_pts);
	template void index_intersection(dt::Delaunay<double>& fg, dt::Delaunay<double>& bg, std::vector<dt::Triangle<double>> &inter_triangles_ans, std::vector<dt::Vector2<double>>& new_pts);

	template<typename T>
	void intersection(dt::Delaunay<T>& fg, dt::Delaunay<T>& bg, std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>>& new_pts);
	template void intersection(dt::Delaunay<double>& fg, dt::Delaunay<double>& bg, std::vector<dt::Triangle<double>> &inter_triangles_ans, std::vector<dt::Vector2<double>>& new_pts);

	/*};*/
}
#endif
