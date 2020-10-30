#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>
#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"
//#include"numeric.h"
#include "intersection.h"

namespace dt{
	//template<typename T>

	bool pt_cmp(const dt::Vector2<double>& p1, const dt::Vector2<double>& p2) {
		if (dt::almost_equal(p1.x, p2.x)) {
			if (p1.y < p2.y)//>
				return true;
			else//>
				return false;
		}
		else if (p1.x < p2.x) //<
			return true;
		else //>
			return false;
	}

	template<typename T>
	T min_value(T a, T b) {
		if (a < b)
			return a;
		else
			return b;
	}
	template<typename T>
	T max_value(T a, T b) {
		if (a < b)
			return b;
		else
			return a;
	}
	
	/*template<typename T>
	bool almost_equal(T a, T b, T eps = 1e-5) {
		if (abs(a - b) <= eps)
			return true;
		else
			return false;

	}*/


	template<typename T>
	bool PointInPolygon(dt::Vector2<T> p, std::vector<dt::Vector2<T>> &  ptPolygon) {
		int nCount = ptPolygon.size();
		int nCross = 0;
		for (int i = 0; i < nCount; i++) {
			auto p1 = ptPolygon[i];//当前节点
			auto p2 = ptPolygon[(i + 1) % nCount];//下一个节点

												  // 求解 y=p.y 与 p1p2 的交点

			if (almost_equal(p1.y, p2.y)) // p1p2 与 y=p0.y平行
				continue;

			if (p.y < min_value(p1.y, p2.y)) // 交点在p1p2延长线上
				continue;
			if (p.y >= max_value(p1.y, p2.y)) // 交点在p1p2延长线上
				continue;

			// 从P发射一条水平射线 求交点的 X 坐标 ------原理: ((p2.y-p1.y)/(p2.x-p1.x))=((y-p1.y)/(x-p1.x))
			//直线k值相等 交点y=p.y
			T x = (T)(p.y - p1.y) * (T)(p2.x - p1.x) / (T)(p2.y - p1.y) + p1.x;

			if (x >= p.x)//if包括边界点=
				nCross++; // 只统计单边交点
		}

		// 单边交点为偶数，点在多边形之外 ---
		return (nCross % 2 == 1);

	}

	//void Intersection::print_edges(const std::vector<dt::Edge<double>> edges) {
	//	for (auto iter = edges.begin(); iter != edges.end(); iter++) {
	//		std::cout << "iter v:(" << iter->v->x << ", " << iter->v->y << ")" << std::endl;
	//		std::cout << "iter w:(" << iter->w->x << ", " << iter->w->y << ")" << std::endl;

	//	}

	//}

	//void Intersection::print_vec(std::vector<double> vec) {
	//	for (auto iter = vec.begin(); iter != vec.end(); iter++) {
	//		std::cout << *iter << std::endl;
	//	}

	//}
	template<typename T>
	bool isRectCross(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> q1, dt::Vector2<T> q2) {
		bool ret = std::min(p1.x, p2.x) <= std::max(q1.x, q2.x) && std::min(p1.y, p2.y) <= std::max(q1.y, q2.y)
			&& std::min(q1.x, q2.x) <= std::max(p1.x, p2.x) && std::min(q1.y, q2.y) <= std::max(p1.y, p2.y);
		return ret;
	}
	template<typename T>
	T cross_multi(dt::Vector2<T> pp0, dt::Vector2<T>pp1, dt::Vector2<T>pp2)//pp0为公用点，叉积计算
	{
		return (pp0.x - pp1.x)*(pp0.y - pp2.y) - (pp0.x - pp2.x)*(pp0.y - pp1.y);

	}
	template<typename T>
	bool isLineSegmentCross(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> q1, dt::Vector2<T> q2) {
		//2.计算边与边的交点；（如果边重叠/（统计，q1/q2刚好在p1p2上，统计；

		if (cross_multi(p1, q1, p2)*cross_multi(p1, p2, q2) < 0) return false;//=0边重叠或者某个端点在另一线段上。
		if (cross_multi(q1, p1, q2)*cross_multi(q1, q2, p2) < 0) return false;

		return true;

	}
	template<typename T>
	T cross_multi(dt::Vector2<T> vec1, dt::Vector2<T>vec2) {
		return vec1.x * vec2.y - vec2.x*vec1.y;
	}

	template<typename T>
	bool GetCrossPoint(dt::Edge<T> cur_edge1, dt::Edge<T>cur_edge2, dt::Vector2<T>& cross_pt) {
		//1.快速排斥实验(对角线对应的矩形是否相交)
		dt::Vector2<T> p1 = *(cur_edge1.v), p2 = *(cur_edge1.w), q1 = *(cur_edge2.v), q2 = *(cur_edge2.w);

		if (isRectCross(p1, p2, q1, q2)) {
			//2.跨立实验
			if (isLineSegmentCross(p1, p2, q1, q2)) {
				//s1: p1p2, s2:q1q2
				//https://blog.csdn.net/zhouzi2018/article/details/80549528
				dt::Vector2<T> base_vec{ q2.x - q1.x, q2.y - q1.y };
				dt::Vector2<T> hypo_vec1{ p1.x - q1.x, p1.y - q1.y };
				dt::Vector2<T> hypo_vec2{ p2.x - q1.x, p2.y - q1.y };

				T d1 = 1.0*abs(cross_multi(base_vec, hypo_vec1)) / base_vec.norm2();
				T d2 = 1.0*abs(cross_multi(base_vec, hypo_vec2)) / base_vec.norm2();
				T t = d1 / (d1 + d2);

				T xx = p1.x + (p2.x - p1.x)*t;
				T yy = p1.y + (p2.y - p1.y)*t;

				cross_pt.setPoint(xx, yy);
				return true;
			}
		}
		return false;



	}
	//


	//3.去重复
	template<typename T>

	void del_repeat(std::vector<dt::Vector2<T>> &ans_tmp, std::vector<dt::Vector2<T>> &ans) {
		//比较相邻项，如果差不多，即重复，去除前一项。
		if (ans_tmp.size() <= 1)//不用去重
		{
			ans = ans_tmp;
			return;
		}
		sort(ans_tmp.begin(), ans_tmp.end(), pt_cmp);


		//ans.push_back(ans_tmp[0]);
		dt::Vector2<T> last = ans_tmp[0];
		for (int i = 1; i < ans_tmp.size(); i++) {
			//比较第i, i+1项;
			if (dt::almost_equal(ans_tmp[i].x, last.x) && dt::almost_equal(ans_tmp[i].y, last.y)) {
				;//空操作
				//if (ans_tmp[i].index < last.index)//比较index;//取较小的index;
				//{
				//	last = ans_tmp[i];
				//}
			}
			else {//
				ans.push_back(last);
				last = ans_tmp[i];
			}
		}
		//
		ans.push_back(last);
		//return ans;
	}

	//4.逆时针或顺时针排序
	//若点A大于B,即点A在点B的顺时针方向，返回true,否则返回false;//
	template<typename T>
	bool CSortCmp(const dt::Vector2<T> &a, const dt::Vector2<T> &b, const dt::Vector2<T> & center) {
		//
		/*if (a.x >= 0 && b.x < 0)
		return true;
		if (a.x == 0 && b.x == 0)
		return a.y > b.y;*/
		//向量OA和OB的叉积


		T det = (a.x - center.x)*(b.y - center.y) - (b.x - center.x)*(a.y - center.y);
		if (det < 0)//
			return true;

		//det>0, A在B的顺时针方向，fanhui false,buyong xiugai weizhi
		if (det > 0)
			return false;
		//向量OA和向量OB共线，以距离判断大小。
		T d1 = (a.x - center.x)*(a.x - center.x) + (a.y - center.y)*(a.y - center.y);
		T d2 = (b.x - center.x)*(b.x - center.y) + (b.y - center.y) *(b.y - center.y);
		return d1 > d2;
	}

	template<typename T>
	void ClockwiseSortPoints(std::vector<dt::Vector2<T>> &ans_tmp) {
		if (ans_tmp.size() <= 1)return;//不用排序

									   //求重心
		dt::Vector2<T> center;
		T xx = 0, yy = 0;
		int cnt = ans_tmp.size();
		for (int i = 0; i < cnt; i++) {
			xx = xx + ans_tmp[i].x*1.0 / cnt;
			yy = yy + ans_tmp[i].y*1.0 / cnt;
		}
		center.x = xx;
		center.y = yy;

		//冒泡排序
		for (int i = 0; i<ans_tmp.size() - 1; i++)
			for (int j = 0; j < ans_tmp.size() - i - 1; ++j) {
				if (CSortCmp(ans_tmp[j], ans_tmp[j + 1], center)) {
					dt::Vector2<T> tmp = ans_tmp[j];
					ans_tmp[j] = ans_tmp[j + 1];
					ans_tmp[j + 1] = tmp;
				}
			}
		//
		return;
	}
	template<typename T>
	void intersect_triangles(dt::Triangle<T> t1, dt::Triangle<T> t2, std::vector<dt::Vector2<T>> &ans_tmp) {

		//三角形的相交区域和凸多边形的相交区域计算方式类似。
		//1.计算多边形内部点；
		//2.计算边与边的交点；（如果边重叠/（统计，q1/q2刚好在p1p2上，统计）
		//3.如果没有内部点或者交点，则重叠区域为空，返回；
		//4.对点去重；
		//5.对点进行逆时针排序

		ans_tmp.clear();
		//1.求内部点
		dt::Vector2<T> A1{ t1.a->x, t1.a->y, t1.a->index };
		dt::Vector2<T> B1{ t1.b->x, t1.b->y, t1.b->index };
		dt::Vector2<T> C1{ t1.c->x, t1.c->y, t1.c->index };

		dt::Vector2<T> A2{ t2.a->x, t2.a->y, t2.a->index };
		dt::Vector2<T> B2{ t2.b->x, t2.b->y, t2.b->index };
		dt::Vector2<T> C2{ t2.c->x, t2.c->y, t2.c->index };

		std::vector<dt::Vector2<T>> t1_pts;
		t1_pts.push_back(A1);
		t1_pts.push_back(B1);
		t1_pts.push_back(C1);
		std::vector<dt::Vector2<T>> t2_pts;
		t2_pts.push_back(A2);
		t2_pts.push_back(B2);
		t2_pts.push_back(C2);


		if (PointInPolygon(A1, t2_pts))
			ans_tmp.push_back(A1);

		if (PointInPolygon(B1, t2_pts))
			ans_tmp.push_back(B1);

		if (PointInPolygon(C1, t2_pts))
			ans_tmp.push_back(C1);

		if (PointInPolygon(A2, t1_pts))
			ans_tmp.push_back(A2);

		if (PointInPolygon(B2, t1_pts))
			ans_tmp.push_back(B2);

		if (PointInPolygon(C2, t1_pts))
			ans_tmp.push_back(C2);
		//2.计算边与边的交点；（如果边重叠/（统计，q1/q2刚好在p1p2上，统计；
		std::vector<dt::Edge<T>> edges1;
		edges1.push_back(dt::Edge<T>(A1, B1));
		edges1.push_back(dt::Edge<T>(B1, C1));
		edges1.push_back(dt::Edge<T>(C1, A1));


		std::vector<dt::Edge<T>> edges2;
		edges2.push_back(dt::Edge<T>(A2, B2));
		edges2.push_back(dt::Edge<T>(B2, C2));
		edges2.push_back(dt::Edge<T>(C2, A2));

		
		for (int i = 0; i<edges1.size(); i++)
			for (int j = 0; j < edges2.size(); j++) {
				//
				auto cur_edge1 = edges1[i];
				auto cur_edge2 = edges2[j];
				dt::Vector2<T> cross_pt;
				if (GetCrossPoint(cur_edge1, cur_edge2, cross_pt)) {
									
					ans_tmp.push_back(cross_pt);
				}
			}
		//3.去重复
		std::vector<dt::Vector2<T>> ans;
		del_repeat(ans_tmp, ans);
		//4.逆时针或顺时针排序
		ClockwiseSortPoints(ans);
		for (int i = 0; i < ans.size(); i++) {
			ans[i].index = i;
		}

		ans_tmp = ans;
	}

	template<typename T>
	void PolygonClip(const dt::Delaunay<T> &fg, const dt::Delaunay<T>& bg, std::vector<dt::Vector2<T>> & ans_tmp) {

		//凸多边形的相交区域和凸多边形的相交区域计算：
		//1.计算多边形内部点；
		//2.计算边与边的交点；（如果边重叠/（统计，q1/q2刚好在p1p2上，统计）
		//3.如果没有内部点或者交点，则重叠区域为空，返回；
		//4.对点去重；
		//5.对点进行逆时针排序

		//1.求内部点

		std::vector<dt::Vector2<T>> o_pt1 = fg._orig_vertices;
		std::vector<dt::Vector2<T>> o_pt2 = bg._orig_vertices;
		for (int i = 0; i<o_pt1.size(); i++)
			if (PointInPolygon(o_pt1[i], o_pt2))
				ans_tmp.push_back(o_pt1[i]);
		for (int i = 0; i<o_pt2.size(); i++)
			if (PointInPolygon(o_pt2[i], o_pt1))
				ans_tmp.push_back(o_pt1[i]);

		//2.计算边与边的交点；（如果边重叠/（统计，q1/q2刚好在p1p2上，统计；

		std::vector<dt::Edge<T>> edges1 = fg.getEdges();
		std::vector<dt::Edge<T>> edges2 = bg.getEdges();
		for (int i = 0; i<edges1.size(); i++)
			for (int j = 0; j < edges2.size(); j++) {
				auto cur_edge1 = edges1[i];
				auto cur_edge2 = edges2[j];
				dt::Vector2<T> cross_pt;
				if (GetCrossPoint(cur_edge1, cur_edge2, cross_pt)) {
					
					ans_tmp.push_back(cross_pt);
				}
			}

		//3.去重复
		std::vector<dt::Vector2<T>> ans;
		del_repeat(ans_tmp, ans);//
		//4.逆时针或顺时针排序
		ClockwiseSortPoints(ans);
		for (int i = 0; i < ans.size(); i++) {
			ans[i].index = i;
		}

		ans_tmp = ans;
		return;
	}

	template<typename T>

	T g(std::vector<dt::Vector2<T>> pts, int a, int b, int c) {
		double t;
		//公式: s = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);//叉积
		t = (pts[a].x - pts[c].x)*(pts[b].y - pts[c].y) - (pts[b].x - pts[c].x)*(pts[a].y - pts[c].y);
		return t;
	}
	template<typename T>
	bool judge_convex_concave(std::vector<dt::Vector2<T>> pts) {
		int p_n = pts.size();
		T t;
		int cnt_lt = 0;//t<0 less than 0的次数
		int cnt_gt = 0;//t>0 greater than 0的次数
		for (int i = 0; i < p_n; i++) {
			//
			t = g(pts, i%p_n, (i + 1) % p_n, (i + 2) % p_n);
			//默认没有重叠点
			if (t < 0)//
				cnt_lt++;
			else if (t > 0)
				cnt_gt++;
		}
		if (cnt_gt>0 && cnt_lt >0)//一定是凹多边形
			return false;
		else
			return true;
	}
	template<typename T>
	std::vector<dt::Vector2<T>> regulation_sort_pt_map(std::vector<dt::Vector2<T>> tmp_pts, std::map<int, int> &pt_map) {
		if (tmp_pts.size() <= 1)//不用去重
		{

			return tmp_pts;
		}

		sort(tmp_pts.begin(), tmp_pts.end(), pt_cmp);

		//比较相邻项，如果差不多，即重复。

		int new_cnt = 0;
		std::vector<dt::Vector2<T>> new_pts;
		new_pts.push_back(tmp_pts[0]);
		pt_map[tmp_pts[0].index] = new_cnt;
		new_cnt++;

		for (int i = 0; i < tmp_pts.size() - 1; i++) {
			//比较第i, i+1项;
			if (dt::almost_equal(tmp_pts[i].x, tmp_pts[i + 1].x) && dt::almost_equal(tmp_pts[i].y, tmp_pts[i + 1].y)) {
				pt_map[tmp_pts[i + 1].index] = pt_map[tmp_pts[i].index];//空操作
			}
			else {//
				new_pts.push_back(tmp_pts[i + 1]);
				pt_map[tmp_pts[i + 1].index] = new_cnt;
				new_cnt++;
			}
		}
		//修改索引
		for (int i = 0; i < new_pts.size(); i++)
		{
			new_pts[i].index = i;
		}
		return new_pts;
	}

	//TODO:三角面片对点进行规整；？？？把多个离散 的多边形分开。
	template<typename T>
	void regulation_triangles(std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>> &new_pts) {//修改index；points;
																																  //
																																  //去重复点；重新标索引。
		std::vector<dt::Vector2<T>> tmp_pts;


		int cur_pt_idx = 0;

		for (int i = 0; i < inter_triangles_ans.size(); i++) {


			auto p1 = *(inter_triangles_ans[i].a);
			auto p2 = *(inter_triangles_ans[i].b);
			auto p3 = *(inter_triangles_ans[i].c);
			p1.index = cur_pt_idx;
			cur_pt_idx++;

			p2.index = cur_pt_idx;
			cur_pt_idx++;

			p3.index = cur_pt_idx;
			cur_pt_idx++;

			tmp_pts.push_back(p1);
			tmp_pts.push_back(p2);
			tmp_pts.push_back(p3);
			//构建点到三角形的双向map; pt_index/3 = triangles;

		}
		//pts->new_pts排序，且构建map
		std::map<int, int> pt_map;
		new_pts = regulation_sort_pt_map(tmp_pts, pt_map);

		//3.修改三角形的点的索引值，并修改所有点。（三角形采用深拷贝，所以都得修改）
		for (int i = 0; i < inter_triangles_ans.size(); i++) {
			inter_triangles_ans[i].a->index = pt_map[i * 3 + 0];
			inter_triangles_ans[i].b->index = pt_map[i * 3 + 1];
			inter_triangles_ans[i].c->index = pt_map[i * 3 + 2];
		}
	}
	template<typename T>
	void index_intersection(dt::Delaunay<T>& fg, dt::Delaunay<T>& bg, std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>>& new_pts) {
		//判断顶点是否全部在另一个内部。
		//TODO?delaunay 算法处理多边形是有问题的，查查github。
		////TODO？
		int cnt_fg = 0, cnt_bg = 0;
		std::vector<dt::Vector2<T>> o_pt1 = fg._orig_vertices;
		//
		std::vector<dt::Vector2<T>> o_pt2 = bg._orig_vertices;





		for (int i = 0; i < o_pt1.size(); i++) {
			if (PointInPolygon(o_pt1[i], o_pt2))//TODO???修改第二个参数用edges比points好一些。可以处理多个离散多边形的情况。
				cnt_fg++;
		}
		for (int i = 0; i < o_pt2.size(); i++) {
			if (PointInPolygon(o_pt2[i], o_pt1))
				cnt_bg++;
		}
		//1. 一个在另一个的内部
		if (cnt_fg == o_pt1.size()) {//前者在后者的内部
			inter_triangles_ans = fg._orig_triangles;
			new_pts = fg._orig_vertices;
			return;

		}
		if (cnt_bg == o_pt2.size()) {//后者在前者的内部
			inter_triangles_ans = bg._orig_triangles;
			new_pts = bg._orig_vertices;
			return;
		}

		//2. if (完全没有交点) {
		std::vector<dt::Triangle<T>> ans;
		ans.clear();
		if (cnt_fg == 0 && cnt_bg == 0)
			return;
		//return ans;
		//3.如果都是凸多边形，直接进行求交
		if (judge_convex_concave(o_pt1) && judge_convex_concave(o_pt2)) {//不能判断多个离散的情况下的凹凸判断，点集是顺序输入（最好是逆序）；

			std::vector<dt::Vector2<T>> polygon_pts;
			PolygonClip(fg, bg, polygon_pts);
			dt::Delaunay<T> new_triangulation;
			inter_triangles_ans = new_triangulation.triangulate(polygon_pts);
			//规整索引
			//std::vector<dt::Vector2<double>> new_pts;
			regulation_triangles(inter_triangles_ans, new_pts);

			return;
		}
		//3. 如果只有一方是凸多边形，将另一方化成三角形再进行求交。TODO？

		std::vector<dt::Triangle<T>> ot1 = fg._orig_triangles;
		std::vector<dt::Triangle<T>> ot2 = bg._orig_triangles;

		inter_triangles_ans.clear();
		for (auto iter_fg = ot1.begin(); iter_fg != ot1.end(); ++iter_fg) {
			for (auto iter_bg = ot2.begin(); iter_bg != ot2.end(); ++iter_bg) {
				std::vector<dt::Vector2<T>> polygon_pts;
				polygon_pts.clear();
				try {
					intersect_triangles(*iter_fg, *iter_bg, polygon_pts);// 逆时针/顺时针？；
				}
				catch (std::bad_alloc) {
					//输出polygon_pts情况
					std::cout << "bad_alloc" << std::endl;
				}
				//对polygon_pts的index重排序。？？
				std::cout << "1" << std::endl;
				std::cout << "2" << std::endl;
				dt::Delaunay<T> triangulation_bg;
				const std::vector<dt::Triangle<T>> triangles_inter = triangulation_bg.triangulate(polygon_pts);//因为是凸多边形，所以不需要最后一步操作（排除异常情况true条件）

																													//... vec1,vec2赋值
				inter_triangles_ans.insert(inter_triangles_ans.end(), triangles_inter.begin(), triangles_inter.end());
			}
		}
		//规整索引
		regulation_triangles(inter_triangles_ans, new_pts);
	}
	template<typename T>
	void intersection(dt::Delaunay<T>& fg, dt::Delaunay<T>& bg, std::vector<dt::Triangle<T>> &inter_triangles_ans, std::vector<dt::Vector2<T>>& new_pts) {
		//判断顶点是否全部在另一个内部。
		//TODO?delaunay 算法处理多边形是有问题的，查查github。
		////TODO？
		int cnt_fg = 0, cnt_bg = 0;
		std::vector<dt::Vector2<T>> o_pt1 = fg._orig_vertices;
		std::vector<dt::Vector2<T>> o_pt2 = bg._orig_vertices;





		for (int i = 0; i < o_pt1.size(); i++) {
			if (PointInPolygon(o_pt1[i], o_pt2))//TODO???修改第二个参数用edges比points好一些。可以处理多个离散多边形的情况。
				cnt_fg++;
		}
		for (int i = 0; i < o_pt2.size(); i++) {
			if (PointInPolygon(o_pt2[i], o_pt1))
				cnt_bg++;
		}
		//1. 一个在另一个的内部
		if (cnt_fg == o_pt1.size()) {//前者在后者的内部
			inter_triangles_ans = fg._orig_triangles;
			new_pts = fg._orig_vertices;
			return;
		}
		if (cnt_bg == o_pt2.size()) {//后者在前者的内部
			inter_triangles_ans = bg._orig_triangles;
			new_pts = bg._orig_vertices;
			return;
		}

		//2. if (完全没有交点) {
		std::vector<dt::Triangle<T>> ans;
		ans.clear();
		if (cnt_fg == 0 && cnt_bg == 0)
			return;
		//return ans;
		//3.如果都是凸多边形，直接进行求交
		if (judge_convex_concave(o_pt1) && judge_convex_concave(o_pt2)) {//不能判断多个离散的情况下的凹凸判断，点集是顺序输入（最好是逆序）；
			std::vector<dt::Vector2<T>> polygon_pts;
			PolygonClip(fg, bg, polygon_pts);
			dt::Delaunay<T> new_triangulation;
			inter_triangles_ans = new_triangulation.triangulate(polygon_pts);
			//规整索引
			//std::vector<dt::Vector2<double>> new_pts;
			regulation_triangles(inter_triangles_ans, new_pts);
			return;
		}
		//3. 如果只有一方是凸多边形，将另一方化成三角形再进行求交。TODO？

		std::vector<dt::Triangle<T>> ot1 = fg._orig_triangles;
		std::vector<dt::Triangle<T>> ot2 = bg._orig_triangles;

		inter_triangles_ans.clear();
		for (auto iter_fg = ot1.begin(); iter_fg != ot1.end(); ++iter_fg) {
			for (auto iter_bg = ot2.begin(); iter_bg != ot2.end(); ++iter_bg) {
				std::vector<dt::Vector2<T>> polygon_pts;
				polygon_pts.clear();
				try {
					intersect_triangles(*iter_fg, *iter_bg, polygon_pts);// 逆时针/顺时针？；
				}
				catch (std::bad_alloc) {
					//输出polygon_pts情况
					std::cout << "bad_alloc" << std::endl;
				}
				//对polygon_pts的index重排序。？？
				std::cout << "1" << std::endl;
				std::cout << "2" << std::endl;
				dt::Delaunay<T> triangulation_tmp;
				const std::vector<dt::Triangle<T>> triangles_inter = triangulation_tmp.triangulate(polygon_pts);//因为是凸多边形，所以不需要最后一步操作（排除异常情况true条件）

																													//... vec1,vec2赋值
				inter_triangles_ans.insert(inter_triangles_ans.end(), triangles_inter.begin(), triangles_inter.end());
			}
		}
		//规整索引
		regulation_triangles(inter_triangles_ans, new_pts);
		//最大最小内角，//进行delaunay准则优化（最大-最小内角）
		fg.maxMinInnerAngle(inter_triangles_ans, new_pts);

		
	}
	/*template struct Intersection<double>;*/

}