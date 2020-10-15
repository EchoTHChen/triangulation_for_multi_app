#pragma once
#ifndef H_DELAUNAY
#define H_DELAUNAY

#include "vector2.h"
#include "edge.h"
#include "triangle.h"

#include <vector>
#include <algorithm>
#include <map>
#include <set>
namespace dt {

	//template<typename T>
	struct Edge_Node {
		Edge_Node();
		Edge_Node(const Edge_Node&en_other);
		//public:
		Edge_Node(int n1, int n2);
		int first_idx;
		int second_idx;
		/*bool operator == ( const Edge_Node & en1) const {
		return first_idx == en1.first_idx && second_idx == en1.second_idx;
		}
		*/
		//      bool operator () (const Edge_Node & en1, const Edge_Node & en2) const {
		//	return en1.first_idx == en2.first_idx && en1.second_idx == en2.second_idx;
		//}

		/*bool operator != (const Edge_Node & en) const {
		return first_idx != en.first_idx || second_idx != en.second_idx;
		}
		*/
		bool operator < (const Edge_Node & en) const {
			//return (first_idx + second_idx) < (en.first_idx + en.second_idx);
			//
			//if (first_idx != en.first_idx)
			//	return first_idx < en.first_idx;
			//else if (second_idx != en.second_idx)
			//	return second_idx < en.second_idx;
			//else
			//	return true;//???
			//
			//return (first_idx < en.first_idx) || (first_idx == en.first_idx&&second_idx < en.second_idx);
			if (first_idx < en.first_idx) {
				return true;
			}
			else if (first_idx > en.first_idx) {
				return false;
			}
			else {
				if (second_idx < en.second_idx)
					return true;
				else
					return false;

			}
		}

	};

	struct Delaunay_Vertex {
		//
		
		Delaunay_Vertex(const dt::Vector2<double> pt);
		double x, y;
		int index;//顶点序号
		int mark; //凹凸顶点标志
		double w; //相关三角形的权值
		
		struct Delaunay_Vertex * last;//指向前一个顶点的指针
		struct Delaunay_Vertex * next;//指向下一个顶点的指针

	};
	
	struct EL_Node {

		EL_Node(const dt::Edge_Node& en0, std::vector<int> tri_vec0);
		dt::Edge_Node en;
		std::vector<int> tri_vec;
		
	};

	struct pt_coincide {
		pt_coincide(dt::Vector2<double> p0, int cat) {
			p = p0;
			category = cat;
		}
		dt::Vector2<double> p;

		int category;
		bool operator <(const pt_coincide&pt_other)const {
			if (almost_equal(p.x, pt_other.p.x)) {
				if (p.y < pt_other.p.y) {
					return true;
				}
				else {
					return false;
				}
			}
			else if (p.x < pt_other.p.x) {
				return true;
			}
			else {
				return false;
			}

		}
	};

	/*Edge_Node{

	}*/
	
	//bool PointInPolygon(dt::Vector2<double> p, std::vector<dt::Vector2<double>> &  ptPolygon);
	//template<typename T>
	class EdgeLinkedNode{
		public:
			EdgeLinkedNode(double xx, double yy) :x(xx), y(yy) {
				idx_set.clear();

			}
			std::set<int> idx_set;
			double x;
			double y;
	};
	

	template<typename T>
	class Delaunay
	{

		using Type = T;
		using VertexType = Vector2<Type>;
		using EdgeType = Edge<Type>;
		using TriangleType = Triangle<Type>;

		/*static_assert(std::is_floating_point<Delaunay<T>::Type>::value,
			"Type must be floating-point");
*/
		
	public:
		//append new elements;
		std::vector<TriangleType> _triangles;
		std::vector<TriangleType> _orig_triangles;
		std::vector<TriangleType> _tear_triangles;

		std::vector<EdgeType> _edges;
		std::vector<VertexType> _vertices;
		std::vector<VertexType> _orig_vertices;

		std::map<Edge_Node, std::vector<TriangleType>> _edge_map;
		std::map<Edge_Node, VertexType> _edge_mid_pt;
		std::vector<EdgeLinkedNode>_edge_linked_list;

		const double _pi = 3.14;
		int _total;

	public:

		Delaunay() = default;
		Delaunay(const Delaunay&) = delete;
		Delaunay(Delaunay&&) = delete;
		//const int  remove_outlier_edges_and_triangle();
		//bool Delaunay<T>::PointInPolygon(dt::Vector2<double> p, std::vector<dt::Vector2<double>> &  ptPolygon);
		//const std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices, bool if_polygon);
		//subdivision
		T cross_multi(dt::Vector2<T> pt1, dt::Vector2<T> pt2);
		T min_value(T a, T b);
		T max_value(T a, T b);
		bool almost_equal(T a, T b, T eps = 1e-5);
		bool PointInPolygon(dt::Vector2<T> p, std::vector<dt::Vector2<T>> &  ptPolygon);
		//const std::vector<TriangleType>& loop_subdivision(int degree);
		//std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices);
		//polygon_delaunay_triangulate:
		//
		T dist2(dt::Delaunay_Vertex pt1, dt::Delaunay_Vertex pt2);
		T  getWeight(Delaunay_Vertex*cur_pt);
		
		T getWeight(int pt_idx1, int pt_idx2, int pt_idx3, std::vector<dt::Vector2<T>>points);
		dt::Delaunay_Vertex* getMaxWeight(Delaunay_Vertex*dv_head);
		void addTriangleTML(Delaunay_Vertex*max_pt, std::vector<dt::Triangle<T>>&tml);
		int CountDelaunayVertex(Delaunay_Vertex* dv_head);
		void edge_map_triangles(std::map<dt::Edge_Node, std::vector<int>> &el_map, int tri_idx, dt::Vector2<double> p1, dt::Vector2<double> p2);

		int find_other_pt(int tri_idx, std::vector<dt::Triangle<T>> &tml, int pt_idx1, int pt_idx2);
		T g(std::vector<dt::Vector2<T>> pts, int a, int b, int c);

		bool judge_convex_concave(std::vector<dt::Vector2<T>> pts);

		void el_revise(std::map<dt::Edge_Node, std::vector<int>>&el_map, dt::Edge_Node tmp_en, int tri_idx1, int tri_idx2);
		void maxMinInnerAngle(std::vector<dt::Triangle<T>>& tml, std::vector<dt::Vector2<T>>points);
		std::vector<dt::Triangle<T> > & polygon_delaunay_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles);


		//ear_clipping_triangulate
		bool judge_clockwise(std::vector<dt::Vector2<T>> points);
		T getAreaClock(dt::Vector2<T>p1, dt::Vector2<T> p2, dt::Vector2<T>p3);
		bool collisionTrianglePoint(dt::Vector2<T> a, dt::Vector2<T> b, dt::Vector2<T> c, dt::Vector2<T> point);
		T Distance(dt::Vector2<T> const& v1, dt::Vector2<T> const& v2);
		T determinant(dt::Vector2<T> u, dt::Vector2<T> v);
		std::vector<dt::Triangle<T>>& ear_clipping_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles);
		std::vector<dt::Triangle<T>>& tri_segment_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles);
		bool judge_legal_tri(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> p3);

		std::vector<dt::Triangle<T> > & triangulate(std::vector<dt::Vector2<T>> &points);

		
		//subdivision
		void getNewIndices(int total, std::vector<dt::Triangle<T>> &triangles, std::vector<dt::Vector2<T>>& new_pts);
		T max_vec(std::vector<T> vec);
		void get_area(std::vector<dt::Triangle<T>>triangles, std::vector<T> &area_vec);
		T calc_area(dt::Triangle<T> triangle);
		T calc_edge(dt::Vector2<T> *P1, dt::Vector2<T> *P2);
		void merge_triangles(std::vector<dt::Triangle<T>>&sub_triangles, std::vector<dt::Triangle<T>> triangles);
		void loop_subdivision_quarter(dt::Triangle<T> triangle, std::vector<dt::Triangle<T>> & tri_vec);
		void loop_subdivision_half(dt::Triangle<T> triangle, std::vector<dt::Triangle<T>> & tri_vec);
		std::vector<dt::Triangle<T>>& subdivision(int k);
		
		//tear
		bool pt_in_tri(dt::Vector2<T> p, dt::Triangle<T> tri);
		bool Collinear(dt::Vector2<T> pt, dt::Edge<T> e);//共线
		T dist2(dt::Vector2<T> pt1, dt::Vector2<T>pt2);
		bool almost_equal_pt(dt::Vector2<T>pt1, dt::Vector2<T>pt2);

		void triangulation_restrict(dt::Edge<T> result_segment, dt::Triangle<T>triangle, std::vector<dt::Triangle<T>>&result_triangles);
		void edge_in_tri(dt::Edge<T> edge, dt::Triangle<T> tri, bool &return_flag, dt::Edge<T> &result_segment);

		void line_in_tri(std::vector<dt::Edge<T>>line, dt::Triangle<T> tri, std::vector<dt::Edge<T>>& in_line);	
		void segments_sub_triangles(dt::Triangle<T> tri, std::vector<dt::Edge<T>>& segments, std::vector<dt::Triangle<T> >  &sub_triangles);
		bool intersect_edge_tri(dt::Edge<T>  e, dt::Triangle<T>tri, int & ans1, int &ans2, int &ans3, dt::Vector2<T> &cpt1, dt::Vector2<T> & cpt2, dt::Vector2<T> &cpt3, dt::Edge<T> &ce1, dt::Edge<T> &ce2, dt::Edge<T> &ce3);
		void process_tris_line(dt::Edge<T>segment, std::vector<dt::Triangle<T>>triangles, std::vector<dt::Triangle<T>>&sub_triangles);
		bool outside_all(std::vector<dt::Edge<T>>line, dt::Triangle<T>tri);


		int coincide_intersect_segments(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T>q1, dt::Vector2<T>q2, dt::Vector2<T> &cross_pt, dt::Edge<T> &coincide_edge);
		
		std::vector<dt::Triangle<T>>& polygon_lines(std::vector<std::vector<dt::Edge<T>>>lines);

		
		//bool add_edge_list_map(TriangleType t, dt::Vector2<double> p1, dt::Vector2<double> p2);
		//
		
		//
		//
		
		//void merge_tri(std::vector<dt::Triangle<T>>& triangles, std::vector<dt::Triangle<T>>tmp_tri_vec)

		void print_edge_map(std::map<Edge_Node, std::vector<TriangleType>> em);
		void print_edge_mid_pt(std::map<Edge_Node, VertexType> edge_mid_pt);
		void print_new_triangles(const std::vector<TriangleType> triangles) const;

		void print_triangle(TriangleType triangle);


		const std::vector<TriangleType>& getTriangles() const;
		const std::vector<EdgeType>& getEdges() const;
		const std::vector<VertexType>& getVertices() const;

		void setTriangles(std::vector<TriangleType> triangles, std::vector<VertexType> & new_pts);

		Delaunay& operator=(const Delaunay&) = delete;
		Delaunay& operator=(Delaunay&&) = delete;
	};

} // namespace dt

#endif
