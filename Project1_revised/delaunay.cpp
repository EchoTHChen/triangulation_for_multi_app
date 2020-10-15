#include "delaunay.h"
//#include "Polygon_cal.h"

#include<cmath>
#include <typeinfo>
#include<algorithm>
namespace dt {
	bool pt_cmp_index(const dt::Vector2<double>& p1, const dt::Vector2<double>& p2) {
		if (dt::almost_equal(p1.x, p2.x)) {
			if (dt::almost_equal(p1.y, p2.y)) {
				if (p1.index > p2.index)
					return true;
				else
					return false;
			}
			else if (p1.y < p2.y)//>
				return true;
			else//>
				return false;
		}
		else if (p1.x < p2.x) //<
			return true;
		else //>
			return false;
	}
	bool index_cmp(const dt::Vector2<double>& p1, const dt::Vector2<double>& p2) {
		if (p1.index < p2.index)
			return true;
		else
			return false;
		
	}

	Edge_Node::Edge_Node() {
		first_idx = -1;
		second_idx = -1;

	}
	Edge_Node::Edge_Node(const Edge_Node&en_other) {
		first_idx = en_other.first_idx;
		second_idx = en_other.second_idx;

	}
	//public:
	Edge_Node::Edge_Node(int n1, int n2) {
		if (n1 < n2) {
			first_idx = n1;
			second_idx = n2;
		}
		else {
			first_idx = n2;
			second_idx = n1;
		}
	}


	EL_Node::EL_Node(const dt::Edge_Node &en0, std::vector<int> tri_vec0)
	{
		en = dt::Edge_Node{ en0.first_idx, en0.second_idx };
		tri_vec = tri_vec0;

	}
	Delaunay_Vertex::Delaunay_Vertex(const dt::Vector2<double> pt):
		index(pt.index), x(pt.x), y(pt.y)
	{


	}
	template<typename T>
	T Delaunay<T>::min_value(T a, T b) {
		if (a < b)
			return a;
		else
			return b;
	}
	template<typename T>
	T Delaunay<T>::max_value(T a, T b) {
		if (a < b)
			return b;
		else
			return a;
	}

	template<typename T>
	bool Delaunay<T>::almost_equal(T a, T b, T eps = 1e-5) {
		if (abs(a - b) <= eps)
			return true;
		else
			return false;

	}
	template<typename T>
	T Delaunay<T>::cross_multi(dt::Vector2<T> pt1, dt::Vector2<T> pt2) {
		return pt1.x*pt2.y - pt2.x*pt1.y;

	}


	template<typename T>
	bool Delaunay<T>::PointInPolygon(dt::Vector2<T> p, std::vector<dt::Vector2<T>> &  ptPolygon) {
		int nCount = ptPolygon.size();
		int nCross = 0;
		for (int i = 0; i < nCount; i++) {
			auto p1 = ptPolygon[i];//��ǰ�ڵ�
			auto p2 = ptPolygon[(i + 1) % nCount];//��һ���ڵ�

												  // ��� y=p.y �� p1p2 �Ľ���

			if (almost_equal(p1.y, p2.y)) // p1p2 �� y=p0.yƽ��
				continue;

			if (p.y < min_value(p1.y, p2.y)) // ������p1p2�ӳ�����
				continue;
			if (p.y >= max_value(p1.y, p2.y)) // ������p1p2�ӳ�����
				continue;

			// ��P����һ��ˮƽ���� �󽻵�� X ���� ------ԭ��: ((p2.y-p1.y)/(p2.x-p1.x))=((y-p1.y)/(x-p1.x))
			//ֱ��kֵ��� ����y=p.y
			double x = (double)(p.y - p1.y) * (double)(p2.x - p1.x) / (double)(p2.y - p1.y) + p1.x;

			if (x >= p.x)//if�����߽��=
				nCross++; // ֻͳ�Ƶ��߽���
		}

		// ���߽���Ϊż�������ڶ����֮�� ---
		return (nCross % 2 == 1);

	}


	//template<typename T>
	//void Delaunay<T>::print_edge_map(std::map<Edge_Node, std::vector<TriangleType>> em) {
	//	std::cout << "���edge_map���" << std::endl;
	//	int cnt = 0;
	//	for (auto it = em.begin(); it != em.end(); it++) {
	//		//std::cout << "cnt:" << cnt << std::endl;

	//		std::cout << it->first.first_idx << std::endl;
	//		std::cout << it->first.second_idx << std::endl;
	//		cnt += 1;


	//	}
	//}
	template<typename T>
	void Delaunay<T>::print_edge_mid_pt(std::map<Edge_Node, VertexType> emt) {
		std::cout << "���edge_mid_pt���" << std::endl;
		int cnt = 0;
		for (auto it = emt.begin(); it != emt.end(); it++) {
			//std::cout << "cnt:" << cnt << std::endl;
			//std::cout << "end point indices:" << std::endl;
			std::cout << it->first.first_idx << std::endl;
			std::cout << it->first.second_idx << std::endl;
			//std::cout << "point position:" << std::endl;
			std::cout << it->second.x << std::endl;
			std::cout << it->second.y << std::endl;

			cnt += 1;


		}
	}

	template<typename T>
	void Delaunay<T>::print_new_triangles(const std::vector<TriangleType> triangles) const {
		auto iter = triangles.begin();
		std::cout << "***********start**********" << std::endl;
		std::cout << "���triangles:" << std::endl;
		int cnt = 0;
		for (int i=0; i< triangles.size(); i++) {
			//std::cout << "cnt:" << cnt << std::endl;
			auto p1 = triangles[i].a;
			auto p2 = triangles[i].b;
			auto p3 = triangles[i].c;
			std::cout << "p1:(" << p1->x << "," << p1->y <<")"<< std::endl;
			std::cout << "p2:(" << p2->x << "," << p2->y << ")" << std::endl;
			std::cout << "p3:(" << p3->x << "," << p3->y << ")" << std::endl;
			cnt += 1;
		}

		std::cout << "***********end**********" << std::endl;

	}
	template<typename T>
	void Delaunay<T>::print_triangle(TriangleType triangle) {
		auto p1 = triangle.a;
		auto p2 = triangle.b;
		auto p3 = triangle.c;
		std::cout << "triangle:(" << p1->x << "," << p1->y << ")" << ",";
		std::cout << "(" << p2->x << "," << p2->y << ")" << ",";
		std::cout << "(" << p3->x << "," << p3->y << ")" << "," << std::endl;



	}

	//template<typename T>
	//bool Delaunay<T>::PointInPolygon(dt::Vector2<double> p, std::vector<dt::Vector2<double>> &  ptPolygon) {
	//	int nCount = ptPolygon.size();
	//	int nCross = 0;
	//	for (int i = 0; i < nCount; i++) {
	//		auto p1 = ptPolygon[i];//��ǰ�ڵ�
	//		auto p2 = ptPolygon[(i + 1) % nCount];//��һ���ڵ�

	//										 // ��� y=p.y �� p1p2 �Ľ���

	//		if (almost_equal(p1.y, p2.y)) // p1p2 �� y=p0.yƽ��
	//			continue;

	//		if (p.y < min_value(p1.y, p2.y)) // ������p1p2�ӳ�����
	//			continue;
	//		if (p.y >= max_value(p1.y, p2.y)) // ������p1p2�ӳ�����
	//			continue;

	//		// ��P����һ��ˮƽ���� �󽻵�� X ���� ------ԭ��: ((p2.y-p1.y)/(p2.x-p1.x))=((y-p1.y)/(x-p1.x))
	//		//ֱ��kֵ��� ����y=p.y
	//		double x = (double)(p.y - p1.y) * (double)(p2.x - p1.x) / (double)(p2.y - p1.y) + p1.x;

	//		if (x > p.x)
	//			nCross++; // ֻͳ�Ƶ��߽���
	//	}

	//	// ���߽���Ϊż�������ڶ����֮�� ---
	//	return (nCross % 2 == 1);

	//}

	
//	template<typename T>
//	const std::vector<typename Delaunay<T>::TriangleType>&
//		Delaunay<T>::triangulate(std::vector<VertexType> &vertices, bool if_polygon)
//
//	{//std::vector<dt::Vector2<double>> 
//		// Store the vertices locally
//		if (vertices.size() <= 0)
//			return _triangles;
//		_vertices = vertices;
//		_total = vertices.size();
//
//		//initialize _edge_linked_list and _edge_map;
//		for (std::size_t i = 0; i < vertices.size(); ++i) {
//			_edge_linked_list.push_back(EdgeLinkedNode{ vertices[i].x , vertices[i].y });
//			
//		}
//
//		// Determinate the super triangle
//		T minX = vertices[0].x;
//		T minY = vertices[0].y;
//		T maxX = minX;
//		T maxY = minY;
//
//		for (std::size_t i = 0; i < vertices.size(); ++i)
//		{
//			if (vertices[i].x < minX) minX = vertices[i].x;
//			if (vertices[i].y < minY) minY = vertices[i].y;
//			if (vertices[i].x > maxX) maxX = vertices[i].x;
//			if (vertices[i].y > maxY) maxY = vertices[i].y;
//		}
//
//		const T dx = maxX - minX;
//		const T dy = maxY - minY;
//		const T deltaMax = std::max(dx, dy);
//		const T midx = (minX + maxX) / 2;
//		const T midy = (minY + maxY) / 2;
//
//		const VertexType p1(midx - 20 * deltaMax, midy - deltaMax);
//		const VertexType p2(midx, midy + 20 * deltaMax);
//		const VertexType p3(midx + 20 * deltaMax, midy - deltaMax);
//
//		// Create a list of triangles, and add the supertriangle in it
//		_triangles.push_back(TriangleType(p1, p2, p3));
//
//		for (auto p = begin(vertices); p != end(vertices); p++)
//		{
//			std::vector<EdgeType> polygon;
//
//			for (auto & t : _triangles)
//			{
//				if (t.circumCircleContains(*p))
//				{
//					t.isBad = true;
//					polygon.push_back(Edge<T>{*t.a, *t.b});
//					polygon.push_back(Edge<T>{*t.b, *t.c});
//					polygon.push_back(Edge<T>{*t.c, *t.a});
//				}
//			}
//
//			_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](TriangleType &t) {
//				return t.isBad;
//			}), end(_triangles));
//
//			for (auto e1 = begin(polygon); e1 != end(polygon); ++e1)
//			{
//				for (auto e2 = e1 + 1; e2 != end(polygon); ++e2)
//				{
//					if (almost_equal(*e1, *e2))
//					{
//						e1->isBad = true;
//						e2->isBad = true;
//					}
//				}
//			}
//
//			polygon.erase(std::remove_if(begin(polygon), end(polygon), [](EdgeType &e) {
//				return e.isBad;
//			}), end(polygon));
//
//			for (const auto e : polygon)
//				_triangles.push_back(TriangleType(*e.v, *e.w, *p));
//
//		}
//
//		_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3](TriangleType &t) {
//			return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
//		}), end(_triangles));
//
//	
//		if(if_polygon){
//			for (auto t = _triangles.begin(); t != _triangles.end(); ) {
//				auto p1 = t->a;
//				auto p2 = t->b;
//				auto p3 = t->c;
//				
//				dt::Vector2<double> Center = dt::Vector2<double>{ (p1->x + p2->x + p3->x) / 3, (p1->y + p2->y + p3->y) / 3 };
//				if (!PointInPolygon(Center, vertices)) {
//					t=_triangles.erase(t);
//				}
//				else
//					t++;
//				//std::cout << "In iteration" << std::endl;
//			
//			}
//		}
//		
//		for (const auto t : _triangles)
//		{
//			_edges.push_back(Edge<T>{*t.a, *t.b});
//			_edges.push_back(Edge<T>{*t.b, *t.c});
//			_edges.push_back(Edge<T>{*t.c, *t.a});
//		}
//		_orig_triangles = _triangles;
//		_orig_vertices = _vertices;
//
//		return _triangles;
//	}
//	//����ϸ��
//	template<typename T>
//	const std::vector<typename Delaunay<T>::TriangleType>& Delaunay<T>::loop_subdivision(int k) {
//		
//		for (int i = 0; i < k; i++) {
//			//ϸ�ִ���
//			
//			//�������ݽṹ_edge_linked_list��_edge_map
//
//			for (auto t : _triangles) {
//				//�����µ�������
//				auto p1 = *t.a;
//				auto p2 = *t.b;
//				auto p3 = *t.c;
//				
//
//				bool add_result=add_edge_list_map(t, p1, p2);
//				add_result=add_edge_list_map(t, p2, p3);
//				
//				add_result=add_edge_list_map(t, p1, p3);
//				
//
//				
//
//			}
//			//�������ж���
//			std::vector<dt::Vector2<double>> new_vertices;
//			new_vertices.clear();
//			for(int j = 0; j < _vertices.size(); j++) {
//				dt::Vector2<double> cur_p = _vertices[j];
//				dt::Vector2<double> new_p;
//				int degree = _edge_linked_list[j].idx_set.size();
//				double tmp_number = 3 / 8 + 1 / 4 * cos(2 * _pi / degree);
//				double beta = 1 / degree*(5 / 8 - tmp_number*tmp_number);
//				
//				//��������j�����ڵ�
//				new_p.x = (1 - degree*beta)*(cur_p.x);
//				new_p.y = (1 - degree*beta)*(cur_p.y);
//				new_p.index = j;
//
//				for (auto edge_it = _edge_linked_list[j].idx_set.begin(); edge_it != _edge_linked_list[j].idx_set.end(); edge_it++) {
//					new_p.x += beta* (_vertices[*edge_it].x);
//					new_p.y += beta* (_vertices[*edge_it].y);
//
//				}
//				new_vertices.push_back(new_p);
//			}
//			//new_vertices = _vertices;
//			//	dt::Vector2<double> new_p;
//
//
//			auto em = _edge_map;
//			int cnt = 0;
//			//error!
//			for (auto it = em.begin(); it != em.end(); it++) {
//
//				cnt += 1;
//				int idx0 = it->first.first_idx;
//				int idx1 = it->first.second_idx;
//				std::vector<TriangleType> edge_tri = it->second;
//				int tri_number = edge_tri.size();
//				//���ڱߵĶ˵�Ϊv0,v1;
//				dt::Vector2<double> v0= _vertices[idx0],v1= _vertices[idx1];
//				dt::Vector2<double> new_add_p;
//				//std::cout << "v0.x, y:" << v0.x << "," << v1.x << std::endl;
//				if (tri_number == 1) {//�߽綥��
//					//std::cout << "inner:" << std::endl;
//					new_add_p.x = 1.0 / 2 * (v0.x + v1.x);
//					new_add_p.y = 1.0 / 2 * (v0.y + v1.y);
//					}
//				else {
//					//�ڲ�����
//					//��Ե���������v2,v3;
//					dt::Vector2<double> v2;
//					dt::Vector2<double> v3;
//					//find_vertices_for_edge_point;
//					TriangleType tri1 = edge_tri[0];
//					TriangleType tri2 = edge_tri[1];
//
//					if (tri1.a->index != idx0 && tri1.a->index != idx1) {
//						v2 = _vertices[tri1.a->index];
//					}
//					//*tri1.b.index;
//					if (tri1.b->index != idx0 && tri1.b->index != idx1) {
//						//
//						v2 = _vertices[tri1.b->index];
//					}
//
//					//*tri1.c.index;
//					if (tri1.c->index != idx0 && tri1.c->index != idx1) {
//						//
//						v2 = _vertices[tri1.c->index];
//					}
//
//					if (tri2.a->index != idx0 && tri2.a->index != idx1) {
//						//
//						v3 = _vertices[tri2.a->index];
//					}
//					//*tri1.b.index;
//					if (tri2.b->index != idx0 && tri2.b->index != idx1) {
//						//
//						v3 = _vertices[tri2.b->index];
//					}
//
//					//*tri1.c.index;
//					if (tri2.c->index != idx0 && tri2.c->index != idx1) {
//						//
//						v3 = _vertices[tri2.c->index];
//					}
//
//					new_add_p.x = 3.0 / 8 * (v0.x) + 3.0 / 8 * (v1.x) + 1.0 / 8 * (v2.x) + 1.0 / 8 * (v3.x);
//					new_add_p.y = 3.0 / 8 * (v0.y) + 3.0 / 8 * (v1.y) + 1.0 / 8 * (v2.y) + 1.0 / 8 * (v3.y);
//						
//				}
//				new_add_p.index = _total;
//				_total += 1;
//
//				new_vertices.push_back(new_add_p);
//
//				//add_vertices.push_back(new_add_p);//index=_total;_total+=1
//				//��¼�ߵ��µ��м�ӳ���
//				_edge_mid_pt[Edge_Node{idx0, idx1}] = new_add_p ;
//		
//			}
//			////���µ�
//			//error!
//			//_vertices = new_vertices;
//			for (int j = 0; j < _vertices.size(); j++) {
//				//ֱ�Ӹ�ֵ��push_back��������,vertices�ĵ�ַ�仯�ˣ�ָ�����
//				_vertices[j].x=new_vertices[j].x;
//				_vertices[j].y = new_vertices[j].y;
//				_vertices[j].index = new_vertices[j].index;
//
//			}
//
//			
//			for (int j = _vertices.size(); j < new_vertices.size(); j++) {
//				_vertices.push_back(new_vertices[j]);
//			}
//
//			std::vector<TriangleType> new_triangles;
//			int cnt_old=0;
//			for (auto t = _triangles.begin(); t != _triangles.end(); t++ ) {
//				double x1= t->a->x, y1= t->a->y;
//				int index1=t->a->index;
//				double x2 = t->b->x, y2 = t->b->y;
//				int index2 = t->b->index;
//				double x3 = t->c->x, y3 = t->c->y;
//				int index3 = t->c->index;
//
//				cnt_old += 1;
//
//				
//				int q12_idx;
//				if (_edge_mid_pt.find(Edge_Node{ index1, index2 })!= _edge_mid_pt.end())
//					q12_idx = _edge_mid_pt.find(Edge_Node{ index1, index2 })->second.index;
//				else {
//					q12_idx = -1;
//					std::cout << "not found in _edge_mid_pt" << std::endl;
//					return _triangles;
//
//				}
//				int q13_idx = _edge_mid_pt.find(Edge_Node{ index1, index3 })->second.index;
//				int q23_idx = _edge_mid_pt.find(Edge_Node{ index2, index3 })->second.index;
//				new_triangles.push_back(TriangleType(new_vertices[index1], new_vertices[q12_idx], new_vertices[q13_idx]));
//
//				new_triangles.push_back(TriangleType(new_vertices[index2], new_vertices[q12_idx], new_vertices[q23_idx]));
//				new_triangles.push_back(TriangleType(new_vertices[index3], new_vertices[q13_idx], new_vertices[q23_idx]));
//				new_triangles.push_back(TriangleType(new_vertices[q12_idx], new_vertices[q13_idx], new_vertices[q23_idx]));
//
//			}
//
//
//			_triangles = new_triangles;
//			
////
//			//���±�
//			_edges.clear();//һ��������գ���Ϊ����ָ���ԭ�е�_vertices�����ѱ��Զ�����ˣ�_edges�����ݳ�Ϊ��Ұָ�룬���Բ���յĻ����ᱣ��Ұָ�������
//			for (const auto t : _triangles)
//			{
//				_edges.push_back(Edge<T>{*t.a, *t.b});
//				_edges.push_back(Edge<T>{*t.b, *t.c});
//				_edges.push_back(Edge<T>{*t.c, *t.a});
//			}
//			
//			
//			_edge_mid_pt.erase(_edge_mid_pt.begin(), _edge_mid_pt.end());
//			
//			_edge_map.erase(_edge_map.begin(), _edge_map.end());
//
//		
//			_edge_linked_list.erase(_edge_linked_list.begin(), _edge_linked_list.end());
//			
//			for (std::size_t i = 0; i < _vertices.size(); ++i) {
//				_edge_linked_list.push_back(EdgeLinkedNode{ _vertices[i].x , _vertices[i].y });
//
//			}
//			
//			
//		}
//		
//		return _triangles;
//	}

	//template<typename T>
	//bool Delaunay<T>::add_edge_list_map(TriangleType t, dt::Vector2<double> p1, dt::Vector2<double> p2) {
	//	//���p1p2
	//	
	//	int idx1 = p1.index;
	//	int idx2 = p2.index;
	//	std::set<int> tmp_set = _edge_linked_list[idx1].idx_set;



	//	if (tmp_set.find(idx2)!= tmp_set.end()) {//֮ǰ��ӹ�
	//												  //˵���Ѿ�ͳ�ƹ���p1p2;
	//		//std::cout << "have added this edge:"<<idx1<<","<<idx2<<std::endl;
	//		return false;
	//	}
	//	else {
	//		//std::cout << "adding this edge::" << idx1 << ","<<idx2 << std::endl;
	//		//ͳ��p1p2;

	//		//˫��ͼ
	//		_edge_linked_list[idx1].idx_set.insert(idx2);
	//		_edge_linked_list[idx1].x = p2.x;
	//		_edge_linked_list[idx1].y = p2.y;

	//		_edge_linked_list[idx2].idx_set.insert(idx1);
	//		//_edge_linked_list[idx2].pt = p1;
	//		_edge_linked_list[idx2].x = p1.x;
	//		_edge_linked_list[idx2].y = p1.y;

	//		//idx1, idx2��˳����Ҫ
	//		dt::Edge_Node tmp_en = dt::Edge_Node{ idx1, idx2 };

	//		if (_edge_map.find(tmp_en) == _edge_map.end()) {//û�ҵ�
	//														//std::cout << "not found this edge in edge_map!" << std::endl;

	//			std::vector<TriangleType> map_tmp_triangles;
	//			map_tmp_triangles.clear();
	//			map_tmp_triangles.push_back(t);
	//			_edge_map[Edge_Node{ idx1, idx2 }] = map_tmp_triangles;
	//		}
	//		else {
	//			//std::cout << "found this edge in edge_map!" << std::endl;

	//			_edge_map[Edge_Node{ idx1, idx2 }].push_back(t);
	//		}
	//		
	//		return true;
	//	}
	//}
	template<typename T>
	T Delaunay<T>::Distance(dt::Vector2<T> const& v1, dt::Vector2<T> const& v2) // basic distance formula : dist(u(x,y) , v(x',y') ) = sqrt( (x'-x)^2 + (y'-y)^2 )
	{
		T distance = sqrt(pow((v2.x - v1.x), 2) + pow((v2.y - v1.y), 2));
		return distance;
	}

	template<typename T>
	bool Delaunay<T>::collisionTrianglePoint(dt::Vector2<T> a, dt::Vector2<T> b, dt::Vector2<T> c, dt::Vector2<T> point) // checks if a point is within the triangle ABC : the point must be at the left of each edge -> be careful to the winding direction
	{
		dt::Vector2<T> ab{ b.x - a.x , b.y - a.y };
		dt::Vector2<T> bc{ c.x - b.x , c.y - b.y };
		dt::Vector2<T> ca{ a.x - c.x , a.y - c.y };
		dt::Vector2<T> ap{ point.x - a.x , point.y - a.y };
		dt::Vector2<T> bp{ point.x - b.x , point.y - b.y };
		dt::Vector2<T> cp{ point.x - c.x , point.y - c.y };
		//if ((determinant(ab, ap) < 0 && determinant(bc, bp) < 0 && determinant(ca, cp) < 0)||(determinant(ab, ap) > 0 && determinant(bc, bp) > 0 && determinant(ca, cp) > 0))//�ڲ�
		//	return true;
		//
		T res1 = determinant(ab, ap), res2= determinant(bc, bp), res3= determinant(ca, cp);
		if ((res1 < 0 || almost_equal(res1, 0)) && (res2 < 0 || almost_equal(res2, 0)) && (res3 < 0 || almost_equal(res3, 0)))//���ڲ����߽߱���
			return true;
		if ((res1 > 0 || almost_equal(res1, 0)) && (res2 > 0 || almost_equal(res2, 0)) && (res3 > 0 || almost_equal(res3, 0)))//
			return true;
		//
		return false;
	}

	
	
	template<typename T>
	T Delaunay<T>::getAreaClock(dt::Vector2<T>p1, dt::Vector2<T> p2, dt::Vector2<T>p3) {
		dt::Vector2<T> p1p2{ p2.x - p1.x, p2.y - p1.y };
		dt::Vector2<T> p2p3{ p3.x - p2.x, p3.y - p2.y };
		return cross_multi(p1p2, p2p3);
	}
	template<typename T>
	bool Delaunay<T>::judge_clockwise(std::vector<dt::Vector2<T>> points) {
	
	
		T area = 0;
		int n = points.size();
		for (int i = 1; i <= n; i++) {
			std::cout << "i-1:" << (i - 1) % n << std::endl;
			std::cout << "i:" << i % n << std::endl;
			std::cout << "i+1:" << (i + 1) % n << std::endl;

			area += getAreaClock(points[(i - 1) % n], points[i%n], points[(i + 1) % n]);
		}
		if (area > 0)
			return false;
		else
			return true;

	}
	template<typename T>
	T Delaunay<T>::determinant(dt::Vector2<T> u, dt::Vector2<T> v) // basic determinant formula : det(u(x,y) , v(x',y') ) = xy'-x'y 
	{
		T result = u.x*v.y - u.y*v.x;
		return result;
	}
	template<typename T>
	std::vector<dt::Triangle<T>>& Delaunay<T>::ear_clipping_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles)
	{


		bool clockwise_order = judge_clockwise(points);
		if (!clockwise_order)
			std::reverse(points.begin(), points.end());
		std::vector<dt::Triangle<T>> triangles; /* a dynamic array that will store the points of the triangles : if the triangle n is (An Bn Cn), then the points will be stored as [A1,B1,C1,
												A2,B2,C2,
												A3,B3,C3...] */
		std::vector<dt::Vector2<T>> initialPoints = points;

		std::vector<int> indexArray;/*
		for (int i = 0; i < initialPoints.size(); i++)
			indexArray.push_back(i);*/

		if (points.size() < 3) // let's make sure that the user don't feed the function with less than 3 points !
			return result_triangles;
		else
		{
			bool impossibleToTriangulate = false;
			bool triangleFound = true;

			while (points.size() != 0) // run the algorithm until our polygon is empty
			{
				if (!triangleFound) // if we've looped once without finding any ear, the program is stuck, the polygon is not triangulable for our algorithm (likely to be a 8 shape or such self intersecting polygon)
				{
					result_triangles = triangles;
					return result_triangles; // we return the triangle array
				}
				triangleFound = false; // we want to find a new ear at each loop

				//-2
				int n = points.size();
				for (int i = 0; i < n ; i++) // for each 3 consecutive points we check if it's an ear : an ear is a triangle that wind in the right direction and that do not contain any other point of the polygon
				{
					if (triangleFound) break;
					 // if we still didn't find an ear
					
					bool result = false;
					dt::Vector2<T> p1p2{ points[(i + 1)%n].x - points[i%n].x, points[(i + 1)%n].y - points[i%n].y };
					dt::Vector2<T> p2p3{ points[(i + 2)%n].x - points[(i + 1)%n].x, points[(i + 2)%n].y - points[(i + 1)%n].y };
					T cm_ans=determinant(p1p2, p2p3);
					bool coincide = false;

					if (almost_equal(cm_ans, 0)) {
						result = true;
						coincide = true;//�غϵ����
						for (int j(0); j < initialPoints.size(); j++) // we check if there's no point inside it
						{

							if (almost_equal_pt(initialPoints[j], points[(i + 2) % n]) || almost_equal_pt(initialPoints[j], points[(i + 1) % n]) || almost_equal_pt(initialPoints[j], points[i%n])) //���������������:
							{
								continue;
							}
							/*if (j == indexArray[(i + 2) % n] || j == indexArray[(i + 1) % n]|| j == indexArray[i % n])
								continue;*/
							if (collisionTrianglePoint(points[(i + 2) % n], points[(i + 1) % n], points[i%n], initialPoints[j]))
							{
								result = false; // if I got a point in my triangle, then it's not an ear !
							}
						}
						

					}
					else if (cm_ans < 0) // if the triangle winds in the right direction
					{
						result = true;
						for (int j(0); j < initialPoints.size(); j++) // we check if there's no point inside it
						{
							if ( almost_equal_pt(initialPoints[j], points[(i + 2)%n]) || almost_equal_pt(initialPoints[j], points[(i+1)%n]) || almost_equal_pt(initialPoints[j], points[i%n])) //���������������:
							{
								continue;
							}
							/*if (j == indexArray[(i + 2) % n] || j == indexArray[(i + 1) % n] || j == indexArray[i % n])
								continue;
*/
							if (collisionTrianglePoint(points[(i + 2)%n], points[(i + 1)%n], points[i%n], initialPoints[j]))
							{
								result = false; // if I got a point in my triangle, then it's not an ear !
							}
						}
					}

					if (result) // now, we have found an ear :
					{
						triangleFound = true;

						//triangles.push_back(points[i]); // so we add our 3 vec2f to the triangle array : it's one of our triangles !
						//triangles.push_back(points[i + 1]);
						//triangles.push_back(points[i + 2]);
						if( !coincide)//
							triangles.push_back(dt::Triangle<T>{points[i%n], points[(i + 1)%n], points[(i + 2)%n]});

						std::vector<dt::Vector2<T>> bufferArray;
						//std::vector<int> bufferIndexArray;

						for (int j(0); j < points.size(); j++) // then we delete the triangle in the points array : we already know that it's an ear, we don't need it anymore
						{
							if (j != (i + 1)%n) // we copiy all the points in a buffer array except the point we don't want
							{
								bufferArray.push_back(points[j]);
								//bufferIndexArray.push_back(indexArray[j]);//?
							}
						}
						points = bufferArray;
						//indexArray = bufferIndexArray;
					}
					
				}
			}
		}
		
	}

	template<typename T>
	bool Delaunay<T>::judge_legal_tri(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T> p3) {
		//�ж������Ƿ��ߣ�
		dt::Vector2<T> p1p2{ p2.x - p1.x, p2.y - p1.y };
		dt::Vector2<T> p2p3{ p3.x - p2.x, p3.y - p2.y };
		T cm_ans = determinant(p1p2, p2p3);
		if (almost_equal(cm_ans, 0)) {
			return false;

		}
		else {
			return true;
		}
	}

	template<typename T>
	std::vector<dt::Triangle<T>>& Delaunay<T>::tri_segment_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles) {
		dt::Vector2<T> v = points[0];
		dt::Vector2<T> w = points[5];
		dt::Vector2<T> pt1 = points[1];
		dt::Vector2<T> pt2 = points[2];
		dt::Vector2<T> pt3 = points[3];
		
		bool right1 = judge_legal_tri(pt1, pt2, v);
		bool right2 = judge_legal_tri(pt1, pt3, w);
		std::vector<dt::Vector2<T>> tmp_points;
		tmp_points.push_back(v);
		tmp_points.push_back(pt2);
		tmp_points.push_back(pt3);
		tmp_points.push_back(w);
		std::vector<dt::Triangle<T>> tmp_triangles;
		ear_clipping_triangulate(tmp_points,tmp_triangles);
		merge_triangles(result_triangles, tmp_triangles);
		if (right1)
			result_triangles.push_back(dt::Triangle<T> {pt1, pt2, v});
		if (right2)
			result_triangles.push_back(dt::Triangle<T> {pt1, pt3, w});

		return result_triangles;

	}

	template<typename T>
	std::vector<dt::Triangle<T> > & Delaunay<T>::triangulate(std::vector<dt::Vector2<T>> &points) {
		_orig_vertices = points;
		_vertices = points;
		//ear_clipping_triangulate(points, _triangles);
		polygon_delaunay_triangulate(points, _triangles);
		_orig_triangles = _triangles;
		//_triangles = _triangles;
		return _triangles;

	}
	template<typename T>
	T Delaunay<T>::dist2(dt::Delaunay_Vertex pt1, dt::Delaunay_Vertex pt2) {
		return sqrt((pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y));
	}

	template<typename T>
	T  Delaunay<T>::getWeight(dt::Delaunay_Vertex*cur_pt) {//�������������С��

		//dt::Delaulay_Vextex* ;
		//PQR: last, cur, next;
		
		T PR = dist2(*(cur_pt->last), *(cur_pt->next));
		T QR = dist2(*(cur_pt), *(cur_pt->next));
		T PQ = dist2(*(cur_pt->last), *(cur_pt));
		//���Ҷ���;
		T cosQ = (QR*QR + PQ*PQ - PR*PR) / (2 * QR*PQ);
		T cosP = (PR*PR + PQ*PQ - QR*QR) / (2 * PR*PQ);
		T cosR = (PR*PR + QR*QR - PQ*PQ) / (2 * PR*QR);

		//�ҵ���С��//
		int index = -1;
		T max_v = -1e-9;
		if (cosQ > max_v)
		{
			index = 0;
			max_v = cosQ;
		}
		if (cosP > max_v)
		{
			index = 1;
			max_v = cosP;
		}
		if (cosR > max_v)
		{
			index = 2;
			max_v = cosR;
		}
		return acos((max_v > 0.99) ? 0.99 : (max_v < -0.99) ? -0.99 : max_v);//
	}
	template<typename T>
	T Delaunay<T>::getWeight(int pt_idx1, int pt_idx2, int pt_idx3, std::vector<dt::Vector2<T>>points) {
		//dt::Delaulay_Vextex* ;
		//PQR: points[pt_idx1], points[pt_idx2], points[pt_idx3];

		T PR = dist2(points[pt_idx1], points[pt_idx3]);
		T QR = dist2(points[pt_idx2], points[pt_idx3]);
		T PQ = dist2(points[pt_idx1], points[pt_idx2]);
		//���Ҷ���;
		T cosQ = (QR*QR + PQ*PQ - PR*PR) / (2 * QR*PQ);
		T cosP = (PR*PR + PQ*PQ - QR*QR) / (2 * PR*PQ);
		T cosR = (PR*PR + QR*QR - PQ*PQ) / (2 * PR*QR);
		//�ҵ���С��//
		int index = -1;
		T max_v = -1e-9;
		if (cosQ > max_v)
		{
			index = 0;
			max_v = cosQ;
		}
		if (cosP > max_v)
		{
			index = 1;
			max_v = cosP;
		}
		if (cosR > max_v)
		{
			index = 2;
			max_v = cosR;
		}
		//acos((alpha > 0.99) ? 0.99 : (alpha < -0.99) ? -0.99 : alpha);
		return acos((max_v > 0.99) ? 0.99 : (max_v < -0.99) ? -0.99 : max_v);//

	}
	template<typename T>
	dt::Delaunay_Vertex* Delaunay<T>::getMaxWeight(dt::Delaunay_Vertex*dv_head) {
		//�������ҵ����Ȩֵ��
		
		dt::Delaunay_Vertex* cur_pt = dv_head;
		T max_w = -1e-9;
		Delaunay_Vertex* max_pt=NULL;

		while (1) {
			//����

			if (cur_pt->mark == 1) //͹
				if (cur_pt->w > max_w)
				{
					max_w = cur_pt->w;
					max_pt = cur_pt;
				}

			if (cur_pt->next == dv_head)
				break;
			cur_pt = cur_pt->next;
		}
		return max_pt;

	}
	template<typename T>
	void Delaunay<T>::addTriangleTML(dt::Delaunay_Vertex*max_pt, std::vector<dt::Triangle<T>>&tml) {
		dt::Delaunay_Vertex*ls = max_pt->last;
		dt::Delaunay_Vertex*nx = max_pt->next;
		
		dt::Vector2<T> ls_pt = dt::Vector2<T>{ ls->x, ls->y, ls->index };
		dt::Vector2<T> cr_pt = dt::Vector2<T>{ max_pt->x, max_pt->y, max_pt->index };
		dt::Vector2<T> nx_pt = dt::Vector2<T>{ nx->x, nx->y, nx->index };
		dt::Triangle<T> triangle{ls_pt, cr_pt, nx_pt};
		tml.push_back(triangle);		
	}
	template<typename T>
	int Delaunay<T>::CountDelaunayVertex(Delaunay_Vertex* dv_head) {
		Delaunay_Vertex* cur_pt = dv_head;
		int cnt = 0;
		while (1) {
			cnt++;
			if (cur_pt->next == dv_head) {
				break;
			}
			cur_pt = cur_pt->next;
		}
		return cnt;
	}
	template<typename T>
	void Delaunay<T>::edge_map_triangles(std::map<dt::Edge_Node, std::vector<int>> &el_map , int tri_idx, dt::Vector2<double> p1, dt::Vector2<double> p2) {
		//���p1p2		
		int idx1 = p1.index;
		int idx2 = p2.index;
		//std::set<int> tmp_set = _edge_linked_list[idx1].idx_set;

		//	//std::cout << "adding this edge::" << idx1 << ","<<idx2 << std::endl;
		//	//ͳ��p1p2;

		//	//˫��ͼ
		//	_edge_linked_list[idx1].idx_set.insert(idx2);
		//	_edge_linked_list[idx1].x = p2.x;
		//	_edge_linked_list[idx1].y = p2.y;

		//	_edge_linked_list[idx2].idx_set.insert(idx1);
		//	//_edge_linked_list[idx2].pt = p1;
		//	_edge_linked_list[idx2].x = p1.x;
		//	_edge_linked_list[idx2].y = p1.y;

		//idx1, idx2��˳����Ҫ
		if (idx1 > idx2) {
			int tmp = idx1;
			idx1 = idx2;
			idx2 = tmp;
		}
		dt::Edge_Node tmp_en = dt::Edge_Node{ idx1, idx2 };

		if (el_map.find(tmp_en) == el_map.end()) {//û�ҵ�
			//std::vector<dt::Triangle<T>> map_tmp_triangles;
			std::vector<int> tri_vec;
			tri_vec.clear();
			tri_vec.push_back(tri_idx);
			el_map[tmp_en] = tri_vec;
		}
		else {
			//std::cout << "found this edge in edge_map!" << std::endl;
			el_map[tmp_en].push_back(tri_idx);
		}		
	}


	template<typename T>
	int Delaunay<T>::find_other_pt(int tri_idx,std::vector<dt::Triangle<T>> &tml, int pt_idx1, int pt_idx2) {
		//�ҵ���ͬ�Ķ���
		int id1 = tml[tri_idx].a->index;
		int id2 = tml[tri_idx].b->index;
		int id3 = tml[tri_idx].c->index;
		int pt_idx3 = -1;
		if (id1 != pt_idx1 && id1 != pt_idx2) {
			pt_idx3 = id1;
		}
		else if (id2 != pt_idx1 && id2 != pt_idx2) {
			pt_idx3 = id2;
		}
		else if (id3 != pt_idx1 && id3 != pt_idx2) {
			pt_idx3 = id3;
		}
		return pt_idx3;
	}


	template<typename T>
	void Delaunay<T>::el_revise(std::map<dt::Edge_Node, std::vector<int>>&el_map, dt::Edge_Node tmp_en, int tri_idx1, int tri_idx2) {
		for (int k = 0; k < el_map[tmp_en].size(); k++) {
			if (el_map[tmp_en][k] == tri_idx1) {
				el_map[tmp_en][k] = tri_idx2;
				break;
			}
		}
	}
	template<typename T>

	T Delaunay<T>::g(std::vector<dt::Vector2<T>> pts, int a, int b, int c) {
		double t;
		//��ʽ: s = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);//���
		t = (pts[a].x - pts[c].x)*(pts[b].y - pts[c].y) - (pts[b].x - pts[c].x)*(pts[a].y - pts[c].y);
		return t;
	}
	template<typename T>
	bool Delaunay<T>::judge_convex_concave(std::vector<dt::Vector2<T>> pts) {
		int p_n = pts.size();
		T t;
		int cnt_lt = 0;//t<0 less than 0�Ĵ���
		int cnt_gt = 0;//t>0 greater than 0�Ĵ���
		for (int i = 0; i < p_n; i++) {
			//
			t = g(pts, i%p_n, (i + 1) % p_n, (i + 2) % p_n);
			//Ĭ��û���ص���
			if (t < 0)//
				cnt_lt++;
			else if (t > 0)
				cnt_gt++;
		}
		if (cnt_gt>0 && cnt_lt >0)//һ���ǰ������
			return false;
		else
			return true;
	}
	template<typename T>
	void Delaunay<T>::maxMinInnerAngle(std::vector<dt::Triangle<T>>& tml, std::vector<dt::Vector2<T>>points) {
		if (points[0].index != 0) {
			std::reverse(points.begin(), points.end());

		}
		while(1){
			int interchanged_edge = 0;
			std::map<dt::Edge_Node, std::vector<int>> el_map;
			//1.��TML���γ�����������ı߱�EL;

			for (int i = 0; i < tml.size(); i++) {
				//������//
				auto p1 = *(tml[i].a);
				auto p2 = *(tml[i].b);
				auto p3 = *(tml[i].c);


				edge_map_triangles(el_map, i, p1, p2);
				edge_map_triangles(el_map, i, p2, p3);
				edge_map_triangles(el_map, i, p1, p3);
			}

			//2.�Ա߱�EL�е�ÿһ���Ǳ߽�ߣ�������3.
			//3.ȡ���Ըñ�Ϊ�ڽ�������2�����������Σ�����һ���ı��Ρ��������ı����ǰ��ģ������κδ����������ı�����͹�ġ�����2�������ε���С�ڽ�alpha���������Խ��ߣ��������������ε���С�ڽ�beta��
			//���beta>alpha,���Խ����Խ��ߺ���������������ԭ�е�������(�޸�el,�޸ı�-���������ڹ�ϵ)���������������޸�TML�������߼�������interchanged_edge +1;(��ʼֵΪ0)
	
			std::vector<dt::EL_Node> el;//
			std::map<dt::Edge_Node, int> edge2ind;//edge->����el������
			//to_el(el_map, el);
			int cnt = 0;
			for (auto it = el_map.begin(); it != el_map.end(); it++) {
				dt::Edge_Node en_tmp =it->first;

				el.push_back(dt::EL_Node{en_tmp, el_map[en_tmp]});
				edge2ind[en_tmp] = cnt;
				cnt++;
			}
			//��mapת�����б�
			for (int i = 0; i < el.size();i++) {
				if (el[i].tri_vec.size() <= 1) //�߽��
					continue;
				else {
					int pt_idx1 = el[i].en.first_idx;//it->first->first;
					int pt_idx2 = el[i].en.second_idx;
					auto tri_vec = el[i].tri_vec;

					int tri_idx1 = tri_vec[0];
					int tri_idx2 = tri_vec[1];

					int pt_idx3 = find_other_pt(tri_idx1, tml, pt_idx1, pt_idx2);
					int pt_idx4 = find_other_pt(tri_idx2, tml, pt_idx1, pt_idx2);

					//
					std::vector<int> pt_idx_vec;
					pt_idx_vec.push_back(pt_idx1);
					pt_idx_vec.push_back(pt_idx3);
					pt_idx_vec.push_back(pt_idx2);
					pt_idx_vec.push_back(pt_idx4);

					std::vector<dt::Vector2<T>> pt_vec;
					for (int j = 0; j < pt_idx_vec.size(); j++)
						pt_vec.push_back(points[pt_idx_vec[j]]);
					if (judge_convex_concave(pt_vec)) {//͹��
						//���Ž����Խ��ߣ�������С�ڽ�
						//�����������ã�////
						//����pt1pt2pt3, pt1pt2pt4����С�ڽ�alpha
						T alpha1 = getWeight(pt_idx1, pt_idx2, pt_idx3, points);
						T alpha2 = getWeight(pt_idx1, pt_idx2, pt_idx4, points);

						T alpha = alpha1 <= alpha2 ? alpha1 : alpha2;//min(alpha1, alpha2);
						//����pt1pt3pt4, pt2pt3pt4����С�ڽ�beta;
						T beta1 = getWeight(pt_idx1, pt_idx3, pt_idx4, points);
						T beta2 = getWeight(pt_idx2, pt_idx3, pt_idx4, points);
						T beta = beta1 <= beta2 ? beta1 : beta2;
						/*if (almost_equal(beta, alpha)) {
							continue;
						}*/
						if (beta > alpha) {
							//�����Խ���
							
							//�޸�el��
							//�޸�tml��
							//triangle tri1: (pt_idx1, pt_idx2, pt_idx3)->(pt_idx1, pt_idx3, pt_idx4);
							//triangle tri2: (pt_idx1, pt_idx2, pt_idx4)->(pt_idx2, pt_idx3, pt_idx4);
							//tml�������������ηֱ����ԭ������������λ��

							//0.
							tml[tri_idx1] = dt::Triangle<T>{ points[pt_idx1], points[pt_idx3], points[pt_idx4] };
							
							//1.
							tml[tri_idx2] = dt::Triangle<T>{ points[pt_idx2], points[pt_idx3], points[pt_idx4] };

							//el��ɾ��pt1pt2�Ͷ�Ӧ��������;,����pt3pt4���Ͷ�Ӧ����������;
							//el_delete(pt_idx1, pt_idx2, i);

							dt::Edge_Node tmp_en0{ pt_idx1, pt_idx2 };
							el_map.erase(el_map.find(tmp_en0));
							std::vector<int> tmp_tri_vec;
							tmp_tri_vec.push_back(tri_idx1);
							tmp_tri_vec.push_back(tri_idx2);
							el_map[tmp_en0] = tmp_tri_vec;
							el[i] = dt::EL_Node{ tmp_en0, tmp_tri_vec };
							edge2ind.erase(edge2ind.find(tmp_en0));
							//�޸�ԭ�����������εı������������ڹ�ϵ;

							//pt1pt3:tri_idx1(�����޸�)  pt2pt3:tri_idx2;�޸�
							dt::Edge_Node tmp_en1{ pt_idx2, pt_idx3 };
							el_revise(el_map, tmp_en1, tri_idx1, tri_idx2);
							el[edge2ind[tmp_en1]] = dt::EL_Node{ tmp_en1,  el_map[tmp_en1]};

							//pt1pt4:tri_idx1���޸ģ�  pt2pt4:tri_idx2;�������޸ģ�
							dt::Edge_Node tmp_en2{ pt_idx1, pt_idx4 };
							el_revise(el_map, tmp_en2, tri_idx2, tri_idx1);
							el[edge2ind[tmp_en2]] = dt::EL_Node{ tmp_en2,  el_map[tmp_en2] };						
							//�޸ı߽�������
							interchanged_edge++;
						}
					}
					else {//����
						//��������
					}
				}
			}
			//4.���interchanged_edge Ϊ0������������ת����1
			if (interchanged_edge == 0) {
				break;
			}
		}
		
		
	}
	//�����ʷ֣�https://wenku.baidu.com/view/3d92d322a5e9856a561260b7.html
	template<typename T>
	std::vector<dt::Triangle<T> > & Delaunay<T>::polygon_delaunay_triangulate(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles) {
		//_orig_vertices = points;
		////ear_clipping_triangulate(points, _triangles);

		//_orig_triangles = _triangles;
		//return _triangles;
		bool clockwise_order = judge_clockwise(points);
		//1.��ʱ�����򵥶���εĶ��㡣
		if (clockwise_order)
			std::reverse(points.begin(), points.end());
		std::vector<dt::Triangle<T>> triangles; /* a dynamic array that will store the points of the triangles : if the triangle n is (An Bn Cn), then the points will be stored as [A1,B1,C1,
												A2,B2,C2,
												A3,B3,C3...] */
		std::vector<dt::Vector2<T>> initialPoints = points;

		std::vector<int> indexArray;/*
									for (int i = 0; i < initialPoints.size(); i++)
									indexArray.push_back(i);*/
		std::vector<dt::Triangle<T>> tml;
		if (points.size() < 3) // let's make sure that the user don't feed the function with less than 3 points !
			return result_triangles;
		else if (points.size() == 3) {
			dt::Triangle<T>triangle{ points[0], points[1], points[2] };
			triangles.push_back(triangle);
			result_triangles = triangles;
			return result_triangles;
		}
		else

		{

			
			//��������˫������
			dt::Delaunay_Vertex *dv_head = NULL;// , *dv_tail = NULL;
			std::vector<dt::Delaunay_Vertex*>dv_vec;
			for (int i = 0; i < points.size(); i++)
				dv_vec.push_back(new Delaunay_Vertex(points[i]));
			

			//dv_tail = new Delaunay_Vertex(points[points.size()-1]);
			int n = points.size();
			for (int i = 0; i < n; i++)
			{
				if (i == 0) {
					dv_head = dv_vec[i];
				}
				dv_vec[i]->last = dv_vec[(i - 1 + n) % n];
				dv_vec[i]->next = dv_vec[(i + 1) % n];
			}

			bool impossibleToTriangulate = false;
			bool triangleFound = true;
			while(1){
				//2.����ÿ������İ�͹�ԣ�
				//i-1 i, i i+1;
				dt::Delaunay_Vertex* cur_pt = dv_head;

				while(1){
					//cur_pt = cur_pt->next;

					dt::Vector2<T> p1p2{ cur_pt->x - cur_pt->last->x, cur_pt->y - cur_pt->last->y };
					dt::Vector2<T> p2p3{ cur_pt->next->x - cur_pt->x, cur_pt->next->y - cur_pt->y };
					if (cross_multi(p1p2, p2p3)>=0) //͹
						cur_pt->mark = 1;
					else //��
						cur_pt->mark = 0;
					//����
					if (cur_pt->next == dv_head) 
						break;
					cur_pt = cur_pt->next;	
				}
				//3.��ÿ��͹����Q������ǰ�󶥵�P,R��ɵ�������Ϊ������PQR����������PQR������������ϵ��������㣨���㷨2�жϣ���������������ε�Ȩֵ������Щ�����������Ȩֵ���������Σ�����Ϊ������ABC,�������εĶ�����ű��浽TML���У�����������ɾ�����B��
				cur_pt = dv_head;

				while (1) {
					//cur_pt = cur_pt->next;
					//PQR;
					/*dt::Vector2<T> p1p2{ cur_pt->x - cur_pt->last->x, cur_pt->y - cur_pt->last->y };
					dt::Vector2<T> p2p3{ cur_pt->next->x - cur_pt->x, cur_pt->next->y - cur_pt->y };*/
					if (cur_pt->mark == 1) //͹
					{//����
						dt::Vector2<T> cr = { cur_pt->x, cur_pt->y };
						dt::Vector2<T> nx = { cur_pt->next->x, cur_pt->next->y };
						dt::Vector2<T> ls = { cur_pt->last->x, cur_pt->last->y };



						bool result = true;
						for (int j(0); j < initialPoints.size(); j++) // we check if there's no point inside it
						{

							if (almost_equal_pt(initialPoints[j],cr ) || almost_equal_pt(initialPoints[j], nx) || almost_equal_pt(initialPoints[j], ls)) //���������������:
								continue;
						
							/*if (j == indexArray[(i + 2) % n] || j == indexArray[(i + 1) % n]|| j == indexArray[i % n])
							continue;*/
							if (collisionTrianglePoint(ls, cr, nx, initialPoints[j]))
							{
								result = false; // if I got a point in my triangle, then it's not an ear !
								break;
							}
						}
						if (result) {
							cur_pt->w = getWeight(cur_pt);//TODO?
						}
						else {
							cur_pt->w = -1;
						}
					}
					else {
						cur_pt->w = -1;
					}
				



					if (cur_pt->next == dv_head)
						break;
					cur_pt = cur_pt->next;
				}
				
				//max_weight
				Delaunay_Vertex *max_pt = getMaxWeight(dv_head);
				//TML table;//delete max_pt;
				addTriangleTML(max_pt, tml);
				if (max_pt == dv_head) {
					dv_head = max_pt->next;
				}
				Delaunay_Vertex*ls = max_pt->last;
				Delaunay_Vertex*nx = max_pt->next;
				ls->next = nx;
				nx->last = ls;
				delete max_pt;//
				
				
				//4.�������л�����3�����ϵĽ�㣬��ת����2������ת����5
				if (CountDelaunayVertex(dv_head) <= 3)
					break;
					
			}
			//5.�����������������������Ӧ�Ķ���ζ��㹹����һ�������Σ�ɾ�����������3�����㡣
			addTriangleTML(dv_head, tml);
			delete dv_head->last;
			delete dv_head->next;

			delete dv_head;
			
			

			//6.�����-��С�ڽ�׼��ͨ���ֲ��任���õ�Delaunay�����ʷ֡����ο��㷨3��
			maxMinInnerAngle(tml, initialPoints);
		}
		result_triangles = tml;

		//result_triangles = _triangles;
		return result_triangles;
	}


	//subdivision
	template<typename T>
	T Delaunay<T>::calc_edge(dt::Vector2<T> *P1, dt::Vector2<T> *P2) {
		return sqrt((P1->x - P2->x)*(P1->x - P2->x) + (P1->y - P2->y)*(P1->y-P2->y));
	}

	template<typename T>
	T Delaunay<T>::calc_area(dt::Triangle<T> triangle) {
		dt::Vector2<T> *A = triangle.a, *B = triangle.b, *C = triangle.c;

		T a = calc_edge(B, C), b = calc_edge(A, C), c = calc_edge(A, B);//������a,b,c
		T p = (a + b + c) / 2;//p�ǰ��ܳ�
						   //���׹�ʽ
		T S = sqrt(p*(p - a)*(p - b)*(p - c));
		return S;

	}
	template<typename T>
	void Delaunay<T>::get_area(std::vector<dt::Triangle<T>>triangles, std::vector<T> &area_vec) {
		area_vec.clear();
		for (int i = 0; i<triangles.size(); i++) {
			area_vec.push_back(calc_area(triangles[i]));
		}
	}
	template<typename T>
	T Delaunay<T>::max_vec(std::vector<T> vec) {
		T max_v = -1e-9;
		for (int i = 0; i<vec.size(); i++) {
			if (vec[i]>max_v)
				max_v = vec[i];
		}
		return max_v;
	}
	template<typename T>
	void Delaunay<T>::merge_triangles(std::vector<dt::Triangle<T>>&sub_triangles, std::vector<dt::Triangle<T>> triangles) {
		for (int i = 0; i<triangles.size(); i++) {
sub_triangles.push_back(triangles[i]);
		}
	}

	template<typename T>
	void Delaunay<T>::loop_subdivision_quarter(dt::Triangle<T> triangle, std::vector<dt::Triangle<T>> & tri_vec) {
		dt::Vector2<T> *a = triangle.a, *b = triangle.b, *c = triangle.c;
		//��ab�е㡢bc�е㡢ac�е�;

		//���ж�����������TODO
		dt::Vector2<T> ab_mid{ (triangle.a->x + triangle.b->x)*1.0 / 2, (triangle.a->y + triangle.b->y)*1.0 / 2 , -1 };
		dt::Vector2<T> bc_mid{ (triangle.b->x + triangle.c->x)*1.0 / 2, (triangle.b->y + triangle.c->y)*1.0 / 2 , -1 };
		dt::Vector2<T> ac_mid{ (triangle.a->x + triangle.c->x)*1.0 / 2, (triangle.a->y + triangle.c->y)*1.0 / 2 , -1 };

		dt::Triangle<T> tri1{ *a, ab_mid, ac_mid };
		dt::Triangle<T> tri2{ *b, ab_mid, bc_mid };
		dt::Triangle<T> tri3{ *c, bc_mid, ac_mid };
		dt::Triangle<T> tri4{ ab_mid, bc_mid, ac_mid };
		tri_vec.push_back(tri1);
		tri_vec.push_back(tri2);
		tri_vec.push_back(tri3);
		tri_vec.push_back(tri4);

	}

	template<typename T>
	void Delaunay<T>::loop_subdivision_half(dt::Triangle<T> triangle, std::vector<dt::Triangle<T> > & tri_vec) {
		T ab = calc_edge(triangle.a, triangle.b);
		T bc = calc_edge(triangle.b, triangle.c);
		T ac = calc_edge(triangle.a, triangle.c);
		T max_v = -1e-9;
		int ind = -1;
		if (ab > max_v)
		{
			ind = 1;
			max_v = ab;
		}
		if (bc > max_v)
		{
			ind = 2;
			max_v = bc;
		}

		if (ac > max_v)
		{
			ind = 3;
			max_v = ac;
		}
		if (ind == 1) {
			//�з�ab
			dt::Vector2<T> ab_mid{ (triangle.a->x + triangle.b->x)*1.0 / 2, (triangle.a->y + triangle.b->y)*1.0 / 2 };
			dt::Triangle<T> tri1{ *(triangle.a), *(triangle.c), ab_mid };
			dt::Triangle<T> tri2{ *(triangle.b), *(triangle.c), ab_mid };
			tri_vec.push_back(tri1);
			tri_vec.push_back(tri2);

		}
		else if (ind == 2) {
			//bc
			dt::Vector2<T> bc_mid{ (triangle.b->x + triangle.c->x)*1.0 / 2, (triangle.b->y + triangle.c->y)*1.0 / 2 };
			dt::Triangle<T> tri1{ *(triangle.a), *(triangle.b), bc_mid };
			dt::Triangle<T> tri2{ *(triangle.a), *(triangle.c), bc_mid };
			tri_vec.push_back(tri1);
			tri_vec.push_back(tri2);



		}
		else if (ind == 3) {
			dt::Vector2<T> ac_mid{ (triangle.a->x + triangle.c->x)*1.0 / 2, (triangle.a->y + triangle.c->y)*1.0 / 2 };
			dt::Triangle<T> tri1{ *(triangle.b), *(triangle.a), ac_mid };
			dt::Triangle<T> tri2{ *(triangle.b), *(triangle.c), ac_mid };
			tri_vec.push_back(tri1);
			tri_vec.push_back(tri2);
		}

	}


	template<typename T>
	void Delaunay<T>::getNewIndices(int total, std::vector<dt::Triangle<T>> &triangles, std::vector<dt::Vector2<T>>& new_pts) {
		//�������εĵ����ӳ��
		std::vector<dt::Vector2<T>> tmp_pts_sort;

		for (int i = 0; i < triangles.size(); i++) {
			triangles[i].a->sub_idx = i * 3 + 0;
			tmp_pts_sort.push_back(*(triangles[i].a));
			triangles[i].b->sub_idx = i * 3 + 1;
			tmp_pts_sort.push_back(*(triangles[i].b));
			triangles[i].c->sub_idx = i * 3 + 2;
			tmp_pts_sort.push_back(*(triangles[i].c));
			
			/*tmp_pts_sort[i * 3 + 0].sub_idx = i * 3 + 0;
			tmp_pts_sort[i * 3 + 1].sub_idx = i * 3 + 1;
			tmp_pts_sort[i * 3 + 2].sub_idx = i * 3 + 2;*/
		}
		sort(tmp_pts_sort.begin(), tmp_pts_sort.end(), pt_cmp_index);
		//
		
		if (tmp_pts_sort[0].index == -1) {
			tmp_pts_sort[0].index = total;
			total++;
		}
		new_pts.push_back(tmp_pts_sort[0]);

		for (int i = 0; i < tmp_pts_sort.size() - 1; i++) {
			if (almost_equal(tmp_pts_sort[i + 1].x, tmp_pts_sort[i].x) && almost_equal(tmp_pts_sort[i + 1].y, tmp_pts_sort[i].y)){
				tmp_pts_sort[i + 1].index = tmp_pts_sort[i].index;
			}
			else {
				if (tmp_pts_sort[i + 1].index == -1) {
					tmp_pts_sort[i + 1].index = total;
					total++;
				}
				new_pts.push_back(tmp_pts_sort[i+1]);
			}
		}
		for (int i = 0; i < tmp_pts_sort.size(); i++) {
			//������
			int tri_idx = tmp_pts_sort[i].sub_idx / 3;
			int tri_pt_idx = tmp_pts_sort[i].sub_idx % 3;
			if (tri_pt_idx == 0) {
				triangles[tri_idx].a->index = tmp_pts_sort[i].index;
			}
			else if (tri_pt_idx == 1) {
				triangles[tri_idx].b->index = tmp_pts_sort[i].index;
			}
			else 
				triangles[tri_idx].c->index = tmp_pts_sort[i].index;
		}
	}
	template<typename T>
	std::vector<dt::Triangle<T> > & Delaunay<T>::subdivision(int k) {//TODO������������

		//_triangles;area
		//��triangles���������������
		//
		

		for(int i = 0; i<k; i++) {
			std::vector<T> area_vec;
			get_area(_triangles, area_vec);
			T max_s = max_vec(area_vec);
			std::vector<dt::Triangle<T>> new_triangles;

			//if (i == 0) {
			for (int j = 0; j<_triangles.size(); j++) {
				T cur_area = area_vec[j];
				if (cur_area >= 3 * 1.0 / 4 * max_s) {
					//����1/4����
					std::vector<dt::Triangle<T>> tmp_tri_vec;
					loop_subdivision_quarter(_triangles[j], tmp_tri_vec);

					merge_triangles(new_triangles, tmp_tri_vec);

				}
				//else if (cur_area < 3 * 1.0 / 4 * max_s && cur_area >= 3 * 1.0 / 8 * max_s) {
				//	//����1/2����
				//	std::cout << "1/2" << std::endl;
				//	std::vector<dt::Triangle<T>> tmp_tri_vec;
				//	loop_subdivision_half(_triangles[j], tmp_tri_vec);
				//	merge_triangles(new_triangles, tmp_tri_vec);

				//}
				else {

					//������
					std::vector<dt::Triangle<double>> tmp_tri_vec;
					tmp_tri_vec.push_back(_triangles[j]);
					//merge_triangles(new_triangles, tmp_tri_vec);
					new_triangles.push_back(_triangles[j]);
				}
			}
			//}
			//else {
			//	//����1/4����
			//	std::vector<dt::Triangle<double>> tmp_tri_vec;
			//	loop_subdivision_quarter(_triangles[j], tmp_tri_vec);
			//	merge_tri(new_triangles, tmp_tri_vec);
			//}
			//����triangles;
			_triangles = new_triangles;
		}

		//add new indices;
		//�ҵ� _vertices
		//n��ð��;

		int total=_vertices.size();//��ǰ��ĸ���
		//����new_triangles.a, b, c���µ������ӳ��;
		//���е�->�������Ȱ���index����ǰ��
		//����_triangles�ĵ������
		std::vector<dt::Vector2<T>> new_pts;
		getNewIndices(total, _triangles, new_pts);
		sort(new_pts.begin(), new_pts.end(), index_cmp);
		maxMinInnerAngle(_triangles, new_pts);
		

		_vertices = new_pts;
		//����_vertices;
		return _triangles;
	}



	template<typename T>
	bool Delaunay<T>::pt_in_tri(dt::Vector2<T> p, dt::Triangle<T> tri) {
		//�жϵ��Ƿ����������ڲ���
		/*
		���㷨���㷨2���ƣ����Կ����Ƕ��㷨2�ļ򻯣�Ҳ���õ������Ĳ�ˡ����������ε������㰴��˳ʱ�루������ʱ�룩˳����A,B,C������ĳһ��P�������������PA,PB,PC, Ȼ���������������ˣ�^��ʾ��˷��ţ���

		t1 = PA^PB,

		t2 = PB^PC,

		t3 = PC^PA,

		���t1��t2��t3ͬ�ţ�ͬ����ͬ��������ôP���������ڲ����������ⲿ��
		*/
		dt::Vector2<T> *a = tri.a, *b = tri.b, *c = tri.c;
		dt::Vector2<T> pa{ a->x - p.x, a->y - p.y };
		dt::Vector2<T> pb{ b->x - p.x, b->y - p.y };
		dt::Vector2<T> pc{ c->x - p.x, c->y - p.y };
		T t1 = cross_multi(pa, pb);
		T t2 = cross_multi(pb, pc);
		T t3 = cross_multi(pc, pa);
		if ((t1 >= 0 && t2 >= 0 && t3 >= 0) || (t1 <= 0 && t2 <= 0 && t3 <= 0)) {//�ڲ�
			return true;
		}
		else {
			return false;
		}

	}

	template<typename T>
	bool Delaunay<T>::Collinear(dt::Vector2<T> p, dt::Edge<T> e) {//����
						   //q1, q2;
						   //p:p;
		dt::Vector2<T> q1 = *(e.v), q2=*(e.w);


		dt::Vector2<T> p_q1{ q1.x - p.x, q1.y - p.y };
		dt::Vector2<T> q{ q2.x - q1.x, q2.y - q1.y };
		if (almost_equal(cross_multi(p_q1, q), 0))
			return true;
		else
			return false;
	}


	
	template<typename T>
	bool Delaunay<T>::almost_equal_pt(dt::Vector2<T>pt1, dt::Vector2<T>pt2) {
		if (almost_equal(pt1.x, pt2.x) && almost_equal(pt1.y, pt2.y))
			return true;
		else
			return false;
	}
	template<typename T>
	T Delaunay<T>::dist2(dt::Vector2<T> pt1, dt::Vector2<T>pt2) {
		return sqrt((pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y));
	}
	
	template<typename T>
	void Delaunay<T>::triangulation_restrict(dt::Edge<T> result_segment, dt::Triangle<T>triangle, std::vector<dt::Triangle<T>>&result_triangles) {
		//���result_segment���˵㲻�غ�
		if (almost_equal_pt(*(result_segment.v), *(result_segment.w))) {
			//TODO
			//��������
			//

		}
		else {
			////�ҵ������ߵĶ���A(������һ��)
			dt::Vector2<T> pt1, pt2, pt3;
			if (!Collinear(*(triangle.a), result_segment)) {
				pt1 = *(triangle.a);
				//pt2��pt1:vec1,  pt1��start��vec2; pt1��end:vec3;
				//vec1
				// pt2 = triangle.b;
				//����vec1=triangle.b-pt1;
				// vec2 = result_segment.u - pt1;
				// vec3 = result_segment.v - pt1;
				dt::Vector2<T> vec1{ triangle.b->x - pt1.x, triangle.b->y - pt1.y };

				dt::Vector2<T> vec2{ result_segment.v->x - pt1.x, result_segment.v->y - pt1.y };
				dt::Vector2<T> vec3{ result_segment.w->x - pt1.x, result_segment.w->y - pt1.y };

				//vec2 x vec1,
				//vec2 x vec3;
				T cm1 = cross_multi(vec2, vec1);
				T cm2 = cross_multi(vec2, vec3);

				if (almost_equal(cm1, 0)) {
					//�غ�
					pt2 = *(triangle.b);
					pt3 = *(triangle.c);

				}
				else if (cm1*cm2 < 0)//���
				{
					pt2 = *(triangle.b);
					pt3 = *(triangle.c);

				}
				else {
					pt2 = *(triangle.c);
					pt3 = *(triangle.b);
				}
			}
			else if (!Collinear(*(triangle.b), result_segment)) {
				pt1 = *(triangle.b);
				dt::Vector2<T> vec1{ triangle.a->x - pt1.x, triangle.a->y - pt1.y };

				dt::Vector2<T> vec2{ result_segment.v->x - pt1.x, result_segment.v->y - pt1.y };
				dt::Vector2<T> vec3{ result_segment.w->x - pt1.x, result_segment.w->y - pt1.y };

				//vec2 x vec1,
				//vec2 x vec3;
				//vec2 x vec1,
				//vec2 x vec3;
				T cm1 = cross_multi(vec2, vec1);
				T cm2 = cross_multi(vec2, vec3);

				if (almost_equal(cm1, 0)) {
					//�غ�
					pt2 = *(triangle.a);
					pt3 = *(triangle.c);
				}
				else if (cm1*cm2 < 0)//���
				{
					pt2 = *(triangle.a);
					pt3 = *(triangle.c);

				}
				else {
					pt2 = *(triangle.c);
					pt3 = *(triangle.a);
				}
			}
			else {
				pt1 = *(triangle.c);
				dt::Vector2<T> vec1{ triangle.a->x - pt1.x, triangle.a->y - pt1.y };
				dt::Vector2<T> vec2{ result_segment.v->x - pt1.x, result_segment.v->y - pt1.y };
				dt::Vector2<T> vec3{ result_segment.w->x - pt1.x, result_segment.w->y - pt1.y };

				//vec2 x vec1,
				//vec2 x vec3;
				//vec2 x vec1,
				//vec2 x vec3;
				T cm1 = cross_multi(vec2, vec1);
				T cm2 = cross_multi(vec2, vec3);

				if (almost_equal(cm1, 0)) {
					pt2 = *(triangle.a);
					pt3 = *(triangle.b);

				}

				else if (cross_multi(vec2, vec1)*cross_multi(vec2, vec3) < 0)//���
				{
					pt2 = *(triangle.a);
					pt3 = *(triangle.b);

				}
				else {
					pt2 = *(triangle.b);
					pt3 = *(triangle.a);
				}
			}
			//start-end-A, start-pt1-pt2-pt3-pt1-end���������ʷ�
			dt::Triangle<T>one_triangle{*( result_segment.v), *(result_segment.w), pt1 };
			std::vector<dt::Vector2<T>> pt_vec;
			pt_vec.push_back(*(result_segment.v));
			pt_vec.push_back(pt1);
			pt_vec.push_back(pt2);
			pt_vec.push_back(pt3);
			pt_vec.push_back(pt1);
			pt_vec.push_back(*(result_segment.w));


			//ear_clipping_triangulate(pt_vec,result_triangles);//TODO?
			tri_segment_triangulate(pt_vec, result_triangles);

			result_triangles.push_back(one_triangle);
			//result_triangles = tmp_triangles;
		}
	}

	template<typename T>
	void Delaunay<T>::edge_in_tri(dt::Edge<T> e, dt::Triangle<T> tri, bool &return_flag, dt::Edge<T> &result_segment) {
		//��e��tri�󽻼�����
		dt::Vector2<T> *p1 = new dt::Vector2<T>{ *(e.v) };
		dt::Vector2<T> *p2 = new dt::Vector2<T>{ *(e.w) }; //e.w;

		if (!pt_in_tri(*p1, tri) && !pt_in_tri(*p2, tri)) {//e���˵���tri��

			dt::Vector2<T> cpt1, cpt2, cpt3;
			int ans1, ans2, ans3;
			dt::Edge<T> ce1, ce2, ce3;

			bool flag = intersect_edge_tri(e, tri, ans1, ans2, ans3, cpt1, cpt2, cpt3, ce1, ce2, ce3);
			//if(coincide(e, tri)){
			//}
			if (!flag) {//if(���û�н���)continue;
						//û�н���
				return_flag = false;
				// continue;
				return;
			}
			else {////e��triĳ�����غ�;ֱ��return ��غϱ�
				if (ans1 == 1 || ans2 == 1 || ans3 == 1) {//�غ�
														  //��غϱ�;
					T max_length = -1e-9;
					int max_idx = -1;

					if (ans1 == 1) {
						//(tri.a), *(tri.b)
						T cur_length = dist2(*(tri.a), *(tri.b));
						if (cur_length > max_length) {
							max_length = cur_length;
							max_idx = 1;
						}
					}
					if (ans2 == 1) {
						//(tri.a), *(tri.b)
						T cur_length = dist2(*(tri.b), *(tri.c));
						if (cur_length > max_length) {
							max_length = cur_length;
							max_idx = 2;
						}
					}

					if (ans3 == 1) {
						//(tri.a), *(tri.b)
						T cur_length = dist2(*(tri.c), *(tri.a));
						if (cur_length > max_length) {
							max_length = cur_length;
							max_idx = 3;
						}
					}
					switch (max_idx) {
					case 1:
						result_segment.setEdge(*(tri.a), *(tri.b));
						break;
					case 2:
						result_segment.setEdge(*(tri.b), *(tri.c));

						break;
					case 3:
						result_segment.setEdge(*(tri.c), *(tri.a));

						break;
					}
					// continue;
					return_flag = true;
					return;
				}
				//��e��tri�����󽻣�//
				else {//�н���
					int cross_cnt = 0;
					if (ans1 == 2)
						cross_cnt++;
					if (ans2 == 2)
						cross_cnt++;
					if (ans3 == 2)
						cross_cnt++;

					if (cross_cnt == 3)//�������н���)//ĳ������һ�������Ƕ���
					{//ֻ��һ�����
					 //�ų���ͬ�Ľ���
						std::vector<pt_coincide> tmp_vec;
						tmp_vec.push_back(pt_coincide(cpt1, 0));
						tmp_vec.push_back(pt_coincide(cpt2, 1));
						tmp_vec.push_back(pt_coincide(cpt3, 2));
						sort(tmp_vec.begin(), tmp_vec.end());
						//result_segment=dt::Edge<T >> {tmp_vec[0].p, tmp_vec[2].p};
						result_segment.setEdge(tmp_vec[0].p, tmp_vec[2].p);
						// continue;
						return_flag = true;
						return;

					}

					else if (cross_cnt == 2)//
					{
						//
						//����������(������ͬһ�����㣬����)
						if (ans1 == 2 && ans2 == 2)
							//result_segment=dt::Edge<T >> {cpt1, cpt2};
							result_segment.setEdge(cpt1, cpt2);


						if (ans1 == 2 && ans3 == 2) 
							//result_segment=dt::Edge<T >> {cpt1, cpt3};
							result_segment.setEdge(cpt1, cpt3);

						if (ans2 == 2 && ans3 == 2) 
							//result_segment=dt::Edge<T >> {cpt2, cpt3};
							result_segment.setEdge(cpt2, cpt3);
						return_flag = true;
						return;
					}
					else if (cross_cnt == 1) {//һ������
											  //�������������
						return_flag = false;
						return;
					}
				}
			}
		}
		else if (pt_in_tri(*p1, tri) && pt_in_tri(*p2, tri)) {//e���˵���tri��){
			result_segment.setEdge(*(e.v), *(e.w));
			return_flag = true;
			return;
			//continue;
		}
		else {
			dt::Vector2<T> cpt1, cpt2, cpt3;
			int ans1, ans2, ans3;
			dt::Edge<T> ce1, ce2, ce3;

			bool flag = intersect_edge_tri(e, tri, ans1, ans2, ans3, cpt1, cpt2, cpt3, ce1, ce2, ce3);

			int cross_cnt = 0;
			if (ans1 == 2)
				cross_cnt++;
			if (ans2 == 2)
				cross_cnt++;
			if (ans3 == 2)
				cross_cnt++;

			//e��triĳ�����غ�;ֱ��return ��غϱ�
			//continue;
			if (ans1 == 1 || ans2 == 1 || ans3 == 1) {//�غ�
													  //��غϱ�;
				T max_length = -1e-9;
				int max_idx = -1;

				if (ans1 == 1) {
					//(tri.a), *(tri.b)

					T cur_length = dist2(*(ce1.v), *(ce1.w));
					if (cur_length > max_length) {
						max_length = cur_length;
						max_idx = 1;
					}
				}
				if (ans2 == 1) {
					//(tri.a), *(tri.b)
					T cur_length = dist2(*(ce2.v), *(ce2.w));
					if (cur_length > max_length) {
						max_length = cur_length;
						max_idx = 2;
					}
				}

				if (ans3 == 1) {
					//(tri.a), *(tri.b)
					T cur_length = dist2(*(ce3.v), *(ce3.w));
					if (cur_length > max_length) {
						max_length = cur_length;
						max_idx = 3;
					}
				}
				switch (max_idx) {
				case 1:
					//result_segment = ce1;
					result_segment.setEdge(*(ce1.v), *(ce1.w));
					break;
				case 2:
					//result_segment = ce2;
					result_segment.setEdge(*(ce2.v), *(ce2.w));

					break;
				case 3:
					//result_segment = ce3;
					result_segment.setEdge(*(ce3.v), *(ce3.w));

					break;
				}
				return_flag = true;
				return;
			}

			//��e��������tri�ı���
			if (cross_cnt == 1)
			{
				//1�������ֱ����
				dt::Vector2<T> cpt;
				if (ans1 == 2)
					cpt = cpt1;
				else if (ans2 == 2)
					cpt = cpt2;
				else
					cpt = cpt3;
				if (pt_in_tri(*p1, tri)) {
					//result_segment = dt::Edge<T>{ cpt, *p1 };
					result_segment.setEdge(cpt, *p1);
					return_flag = true;
					return;
				}
				if (pt_in_tri(*p2, tri)) {
					//result_segment = dt::Edge<T>{ cpt, *p2};
					result_segment.setEdge(cpt, *p2);
					return_flag = true;
					return;
				}
			}
			else if (cross_cnt == 2)
			{
				//2�����
				//(1)����������ĳ��ͬһ����
				//(2)һ�������ڱ��ϣ��Ƕ��㣩����һ������e�Ķ˵㣬�������Ǳ��ϡ�
				dt::Vector2<T> in_pt;
				if (pt_in_tri(*p1, tri))
					in_pt = *p1;
				if (pt_in_tri(*p2, tri))
					in_pt = *p2;
				if (ans1 == 2 && ans2 == 2) {
					if (almost_equal_pt(cpt1, cpt2))
						//result_segment = dt::Edge<T>{ cpt1, in_pt };
						result_segment.setEdge(cpt1, in_pt);

					else
						result_segment.setEdge(cpt1, cpt2);
						//result_segment = dt::Edge<T>{ cpt1, cpt2 };
					
				}
				else if (ans1 == 2 && ans3 == 2) {
					if (almost_equal_pt(cpt1, cpt3))
						result_segment.setEdge(cpt1, in_pt);

						//result_segment = dt::Edge<T>{ cpt1, in_pt };
					else
						result_segment.setEdge(cpt1, cpt3);

						//result_segment = dt::Edge<T>{ cpt1, cpt3 };
					
				}
				else {
					if (almost_equal_pt(cpt2, cpt3))
						result_segment.setEdge(cpt2, in_pt);
						//result_segment = dt::Edge<T>{ cpt2, in_pt };
					else
						result_segment.setEdge(cpt2, cpt3);
						//result_segment = dt::Edge<T>{ cpt2, cpt3 };
					
				}
				return_flag = true;
				return;
			}
			else if (cross_cnt == 3) {
				//ĳ����������ͬһ���㣬��һ������e�Ķ˵㣬�������Ǳ��ϡ�
				//�ų���ͬ�Ľ���
				std::vector<pt_coincide> tmp_vec;
				tmp_vec.push_back(pt_coincide(cpt1, 0));
				tmp_vec.push_back(pt_coincide(cpt2, 1));
				tmp_vec.push_back(pt_coincide(cpt3, 2));
				sort(tmp_vec.begin(), tmp_vec.end());
				result_segment.setEdge(tmp_vec[0].p, tmp_vec[2].p);
				//result_segment = dt::Edge<T >> {tmp_vec[0].p, tmp_vec[2].p};
				return_flag = true;
				return;
			}
			//û�н��㣬һ���������������
		}

	}

	template<typename T>
	void Delaunay<T>::line_in_tri(std::vector<dt::Edge<T>>line, dt::Triangle<T> tri, std::vector<dt::Edge<T>>& in_line) {
		// std::vector<dt::Edge<double>>&
		for (int i = 0; i<line.size(); i++) {
			bool return_flag;
			dt::Edge<T> result_segment;
			//dt::Edge<T> e = line[i];
			edge_in_tri(line[i], tri, return_flag, result_segment);
			if (return_flag)
				if(!almost_equal_pt(*(result_segment.v), *(result_segment.w)))
					in_line.push_back(result_segment);
		}

	}
	template<typename T>
	void Delaunay<T>::segments_sub_triangles(dt::Triangle<T> tri, std::vector<dt::Edge<T>>& segments, std::vector<dt::Triangle<T> >  &sub_triangles) {
		//�����߶����
		sub_triangles.push_back(tri);
		
		//tmp_triangles.push_back(tri);

		// sub_line_segments = convert_segments(sub_lines);
		for (int k = 0; k<segments.size(); k++) {
			//tmp_triangles.clear();
			std::vector<dt::Triangle<T>>tmp_triangles;
			process_tris_line(segments[k], sub_triangles, tmp_triangles);
			sub_triangles = tmp_triangles;
		}
	}
	template<typename T>
	bool Delaunay<T>::intersect_edge_tri(dt::Edge<T>  e, dt::Triangle<T>tri, int & ans1, int &ans2, int &ans3, dt::Vector2<T> &cpt1, dt::Vector2<T> & cpt2, dt::Vector2<T> &cpt3, dt::Edge<T> &ce1, dt::Edge<T> &ce2, dt::Edge<T> &ce3) {
		//�������һ��Ԫ�����ڲ���һ���н���
		/*if (pt_in_tri(*(e.v), tri) || pt_in_tri(*(e.w), tri))
		{
			return true;
		}			
		else {*/
			//https://blog.csdn.net/weixin_42736373/article/details/84587005
			// dt::Vector2<double> cpt1, cpt2, cpt3;
		ans1 = coincide_intersect_segments(*(e.v), *(e.w), *(tri.a), *(tri.b), cpt1, ce1);
		ans2 = coincide_intersect_segments(*(e.v), *(e.w), *(tri.b), *(tri.c), cpt2, ce2);
		ans3 = coincide_intersect_segments(*(e.v), *(e.w), *(tri.c), *(tri.a), cpt3, ce3);
		if (ans1 > 0 || ans2 > 0 || ans3>0)
			return true;
		else
			return false;

			//1.���e���������Ƿ��غ�
			//p:e	
			//(1). e��ab;
			//(2). e��bc;
			//(3). e��ca;

			//2.e����������
		//}


	}
	template<typename T>
	void Delaunay<T>::process_tris_line(dt::Edge<T>segment, std::vector<dt::Triangle<T>>triangles, std::vector<dt::Triangle<T>>&sub_triangles) {
		//�ҵ�segment��triangles�ڵĲ���;
		//sub_triangles.clear();
	/*	if (almost_equal_pt(*(segment.v), *(segment.w))) {
			sub_triangles = triangles;
			return;
		}*/
		for (int i = 0; i<triangles.size(); i++) {
			bool return_flag;
			dt::Edge<T> result_segment;
			std::vector<dt::Triangle<T>> result_triangles;
			edge_in_tri(segment, triangles[i], return_flag, result_segment);
			if (return_flag && almost_equal_pt(*(result_segment.v), *(result_segment.w))) {
				return_flag = false;
			}
			if (return_flag) {
				triangulation_restrict(result_segment, triangles[i], result_triangles);
				merge_triangles(sub_triangles, result_triangles);
			}
			else {
				//merge_triangles(sub_triangles, triangles[i]);
				sub_triangles.push_back(triangles[i]);
			}
		}
	}



	template<typename T>
	bool Delaunay<T>::outside_all(std::vector<dt::Edge<T>>line, dt::Triangle<T>tri) {

		for (int i = 0; i<line.size(); i++) {
			auto e = line[i];
			dt::Vector2<T> cpt1, cpt2, cpt3;
			int ans1, ans2, ans3;
			dt::Edge<T> ce1, ce2, ce3;
			bool flag = false;
			if (!pt_in_tri(*(tri.a), tri) && !pt_in_tri(*(tri.b), tri) && !pt_in_tri(*(tri.c), tri)) {
				if (!intersect_edge_tri(e, tri, ans1, ans2, ans3, cpt1, cpt2, cpt3, ce1, ce2, ce3) ) {
					flag = true;//������������
				}					
			}
			if (!flag)
			{

				return false;
			}
		}
		return true;
	}

	template<typename T>
	int Delaunay<T>::coincide_intersect_segments(dt::Vector2<T> p1, dt::Vector2<T> p2, dt::Vector2<T>q1, dt::Vector2<T>q2, dt::Vector2<T> &cross_pt, dt::Edge<T> &coincide_edge) {
		dt::Vector2<T> p{ p2.x - p1.x, p2.y - p1.y };
		dt::Vector2<T> q{ q2.x - q1.x, q2.y - q1.y };
		T D = cross_multi(p, q);
		dt::Vector2<T> difference_p1_q1{ q1.x - p1.x , q1.y - p1.y };
		T B1 = cross_multi(difference_p1_q1, q);
		if (almost_equal(fabs(D), 0) && almost_equal(fabs(B1), 0))//���� (p1��q����ֱ����)
		{
			//�ĸ����������
			//�����С��������ͬһ���߶Σ����غϡ�
			std::vector<pt_coincide> tmp_vec;
			tmp_vec.push_back(pt_coincide(p1, 0));
			tmp_vec.push_back(pt_coincide(p2, 0));
			tmp_vec.push_back(pt_coincide(q1, 1));
			tmp_vec.push_back(pt_coincide(q2, 1));
			sort(tmp_vec.begin(), tmp_vec.end());
			if (tmp_vec[0].category == tmp_vec[1].category) {
				return 0;//���غ�
			}
			else {
				// dt::Edge<double> e{ tmp_vec[1].p, tmp_vec[2].p};
			/*	coincide_edge.v = tmp_vec[1].p;
				coincide_edge.w = tmp_vec[2].p;*/
				coincide_edge.setEdge(tmp_vec[1].p, tmp_vec[2].p);



				return 1;//���غϲ���
			}

		}
		else if (almost_equal(fabs(D), 0) && !almost_equal(fabs(B1), 0)) {//ƽ�в�����
			return 0;

		}
		else {//pxq!=0, fabs(D)!=0
			T B2 = cross_multi(difference_p1_q1, p);
			T t = B1 / D;
			T u = B2 / D;
			if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {//�ཻ
				cross_pt.x = p1.x + t*p.x;

				cross_pt.y = p1.y + t*p.y;
				return 2;
			}
			else {//���뵫��ƽ��
				return 0;
			}
		}
	}


	template<typename T>
	std::vector<dt::Triangle<T>>& Delaunay<T>::polygon_lines(std::vector<std::vector<dt::Edge<T>>>lines){
		std::vector<dt::Triangle<T>> triangles = _triangles;
		//new_lines = [];
		for (int i = 0; i<lines.size(); i++) {
			std::vector<dt::Triangle<T>> result_triangles;
			for (int j = 0; j < triangles.size(); j++) {//���ÿ��triangles����������;//
														//triangles;
														//��Ҫ���ɶ����߽��зֱ��жϣ���ֳɶ����ߡ�
														//line->sub_lines;
				if (outside_all(lines[i], triangles[j])) {
					//ȫ�����ⲿ:
					//�������
					//merge_triangles(result_triangles, triangles[j])
					result_triangles.push_back(triangles[j]);
					continue;
				}
				//TODO?
				//����������������ڲ��Ĳ���;
				std::vector<dt::Edge<T>> segments;
				std::vector<dt::Triangle<T>> sub_triangles;
				line_in_tri(lines[i], triangles[j], segments);

				// convert2segments(in_lines, segments)
				segments_sub_triangles(triangles[j], segments, sub_triangles);

				 //else if(in_cnt == line.size()){
				 //	//ȫ�����ڲ�
				 //	//�Ե�ǰ
				 //	bool direct_solve_flag=sub_tri_line_func();				

				 //	if(!direct_solve_flag){//����ֱ�����
				 //		//�����߶����
				 //		//TODO?
				 //		segments_sub_triangles(triangles,  sub_triangles);                                                                                                                                                                                                                                                                                                                                                       )
				 //	}
				 //	else{
				 //		//ֱ�����

				 //		merge_triangles(result_triangles, sub_triangles);
				 //		continue;
				 //	}
				 //}
				 //else {
				 //	//�󽻵㣬��ɶ���߶�
				 //	//��黹ʣ�����������������ڣ�
				 //	int line_cnt=cross();
				 //	if(line_cnt>1){
				 //		//�����������ߣ������߶���⡣

				 //	}
				 //	else{
				 //		//��ʣ�����ڲ��Ĳ��ֽ�����⣺
				 //		//ȫ�����ڲ�
				 //		//�Ե�ǰ
				 //		bool direct_solve_flag=sub_tri_line_func();				

				 //		if(!direct_solve_flag){//����ֱ�����
				 //			//�����߶����
				 //			//TODO?
				 //			segments_sub_triangles(triangles,  sub_triangles);                                                                                                                                                                                                                                                                                                                                                       )
				 //		}
				 //		else{
				 //			//ֱ�����,	//ֱ�ӷֽ�������������⡣

				 //			merge_triangles(result_triangles, sub_triangles);
				 //			continue;
				 //		}

				 //	}

				 //}



				merge_triangles(result_triangles, sub_triangles);
			}
			triangles = result_triangles;

		}

		
		_tear_triangles = triangles;
		return _tear_triangles;
	}
	template<typename T>
	void  Delaunay<T>::setTriangles(std::vector<TriangleType> triangles, std::vector<VertexType> & new_pts) {
		_triangles = triangles;
		_vertices = new_pts;
		_orig_triangles = _triangles;
		_orig_vertices = _vertices;
		_total = _vertices.size();

		//Ϊsub_division��׼��
		_edges.clear();
		for (const auto t : _triangles)
		{
			_edges.push_back(Edge<T>{*t.a, *t.b});
			_edges.push_back(Edge<T>{*t.b, *t.c});
			_edges.push_back(Edge<T>{*t.c, *t.a});
		}

		//initialize _edge_linked_list;
		for (std::size_t i = 0; i < _vertices.size(); ++i) {
			_edge_linked_list.push_back(EdgeLinkedNode{ _vertices[i].x , _vertices[i].y });

		}
	}


	template<typename T>
	const std::vector<typename Delaunay<T>::TriangleType>&
		Delaunay<T>::getTriangles() const
	{
		return _triangles;
	}

	template<typename T>
	const std::vector<typename Delaunay<T>::EdgeType>&
		Delaunay<T>::getEdges() const
	{
		return _edges;
	}

	template<typename T>
	const std::vector<typename Delaunay<T>::VertexType>&
		Delaunay<T>::getVertices() const
	{
		return _vertices;
	}

	/*template class Delaunay<float>;*/
	template class Delaunay<double>;
	//template<typename T>
	//const int Delaunay<T>::remove_outlier_edges_and_triangle(std::vector<dt::Vector2<double>> points){
	

} // namespace dt
