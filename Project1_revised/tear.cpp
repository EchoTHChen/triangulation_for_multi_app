#include<map>
#include<set>
#include<iterator>//inserter��Ҫ
#include"tear.h"

namespace dt{
	
	bool TriSortCmp(int pt_1_index, int pt_2_index, const std::vector<dt::Triangle<double>>&triangles, dt::Vector2<double> tear_pt) {
		//
		/*if (a.x >= 0 && b.x < 0)
		return true;
		if (a.x == 0 && b.x == 0)
		return a.y > b.y;*/
		//����Vec1��Vec2�Ĳ��

		dt::Triangle<double> tri_1 = triangles[pt_1_index];
		double center1_x = (tri_1.a->x + tri_1.b->x + tri_1.c->x) / 3;
		double center1_y = (tri_1.a->y + tri_1.b->y + tri_1.c->y) / 3;
		dt::Vector2<double> Vec1{ center1_x - tear_pt.x, center1_y - tear_pt.y};
		
		dt::Triangle<double> tri_2 = triangles[pt_2_index];
		double center2_x = (tri_2.a->x + tri_2.b->x + tri_2.c->x) / 3;
		double center2_y = (tri_2.a->y + tri_2.b->y + tri_2.c->y) / 3;
		dt::Vector2<double> Vec2{ center2_x - tear_pt.x, center2_y - tear_pt.y };


		




		double det = Vec1.x*Vec2.y - Vec1.y*Vec2.x;
		if (det < 0)//
			return true;

		//det>0, A��B��˳ʱ�뷽��fanhui false,buyong xiugai weizhi
		if (det > 0)
			return false;
		//����OA������OB���ߣ��Ծ����жϴ�С��
		double d1 = Vec1.norm2();// +(a.y - center.y)*(a.y - center.y);
		double d2 = Vec2.norm2();// (b.x - center.x)*(b.x - center.y) + (b.y - center.y) *(b.y - center.y);
		return d1 > d2;
	}


	//
	//std::vector<dt::Vector2<T>>& dfs(std::vector<dt::Triangle<T>> triangles, int idx, int &total, std::vector<bool> have_vis_tri, std::set<int>have_appeared_pre_graph, std::vector<bool>& have_appeared_cur_graph, std::map<Edge_Node, int>edge_map, std::map<int, bool> if_tear_pt, std::map<int, int>map_tri2new_pt) {//�µĵ�����new_triangles��new_pts?

	//	//1.����ǰ��
	//	auto cur_tri = triangles[idx];
	//	//�ȴ����㣬�ٴ���ߣ�
	//	have_vis_tri[idx] = true;

	//	int pt_a = cur_tri.a->index;
	//	int pt_b = cur_tri.b->index;
	//	int pt_c = cur_tri.c->index;
	//	if (if_tear_pt.find(pt_a)!=if_tear_pt.end()) {//��˺�ѵ�
	//		if (have_appeared_pre_graph[pt_a] && have_appeared_cur_graph.find(pt_a])==have_appeared_cur_graph.end()){//��֮ǰ��ͨͼ����ֹ�
	//			//����¶��㣻��¼�¶����ԭ����Ķ�Ӧ��ϵ����¼���޸�������
	//			int new_pt_a = total;
	//			total++;
	//			//(3 * idx + 0)->new_pt_a; //���������ζ��㵽��index��ӳ��
	//			map_tri2new_pt[3*idx+0] = new_pt_a;
	//			have_appeared_cur_graph.insert(pt_a);
	//		}			
	//	}
	//	if (if_tear_pt.find(pt_b) != if_tear_pt.end()) {//��
	//		if (have_appeared_pre_graph[pt_b] && have_appeared_cur_graph.find(pt_b]) == have_appeared_cur_graph.end()) {//��֮ǰ��ͨͼ����ֹ�
	//																													//����¶��㣻��¼�¶����ԭ����Ķ�Ӧ��ϵ����¼���޸�������
	//			int new_pt_b = total;
	//			total++;
	//			//(3 * idx + 1)->new_pt_a; //���������ζ��㵽��index��ӳ��
	//			map_tri2new_pt[3 * idx + 1] = new_pt_b;
	//			have_appeared_cur_graph.insert(pt_b);
	//		}
	//	}
	//	if_tear_pt.find(pt_c) != if_tear_pt.end(){
	//	//if (pt_c��˺�ѵ�) {//��
	//		if (have_appeared_pre_graph[pt_c] && have_appeared_cur_graph.find(pt_c]) == have_appeared_cur_graph.end()) {//��֮ǰ��ͨͼ����ֹ�
	//																													//����¶��㣻��¼�¶����ԭ����Ķ�Ӧ��ϵ����¼���޸�������
	//			int new_pt_c = total;
	//			total++;
	//			//(3 * idx + 2)->new_pt_a; //���������ζ��㵽��index��ӳ��
	//			map_tri2new_pt[3 * idx + 2] = new_pt_c;
	//			have_appeared_cur_graph.insert(pt_c);
	//		}
	//	}

	//	//����˺�ѱ�
	//	//���ĳ������˺�ѱߣ����뵱ǰ��������������Ϊ�����ߵģ�δ���ʵģ����������μ���set<int> cur_chn_tri;
	//	//�����ǰ����������cur_chn_tri����ж����жϣ�˺�ѵ㣺���е����ӣ�

	//	//�����ǰ�����εı�;
	//	

	//	//2.��������������
	//	have_vis_tri[i] = true;
	//	//for (int i = 0; i < ; i++) {//������������
	//	// dfs()   
	//	//}
	//	std::vector<int> tri_vec1 = edge_map[Edge_Node{ pt_a, pt_b }];
	//	std::vector<int> tri_vec2 = edge_map[Edge_Node{ pt_b, pt_c }];
	//	std::vector<int> tri_vec3 = edge_map[Edge_Node{ pt_a, pt_c }];
	//	for (int i = 0; i < tri_vec1.size(); i++) {
	//		int tri_idx = tri_vec1[i];
	//		if (tri_idx != idx) {//��Ϊ��ǰ�����Σ�������������
	//			dfs( triangles, tri_idx, total, have_vis_tri, have_appeared_pre_graph, have_appeared_cur_graph, edge_map, if_tear_pt,  map_tri2new_pt);
	//		}
	//	}
	//	for (int i = 0; i < tri_vec2.size(); i++) {
	//		int tri_idx = tri_vec2[i];
	//		if (tri_idx != idx) {//��Ϊ��ǰ�����Σ�������������
	//			dfs(triangles, tri_idx, total, have_vis_tri, have_appeared_pre_graph, have_appeared_cur_graph, edge_map, if_tear_pt, map_tri2new_pt);
	//		}
	//	}
	//	for (int i = 0; i < tri_vec3.size(); i++) {
	//		int tri_idx = tri_vec3[i];
	//		if (tri_idx != idx) {//��Ϊ��ǰ�����Σ�������������
	//			dfs(triangles, tri_idx, total, have_vis_tri, have_appeared_pre_graph, have_appeared_cur_graph, edge_map, if_tear_pt, map_tri2new_pt);
	//		}
	//	}

	//}
	//template<typename T>
	//bool add_edge_map(int idx, dt::Vector2<T> p1, dt::Vector2<T> p2, std::map<Edge_Node<T>, std::vector<int>> &edge_map) {
	//	//���p1p2

	//	int idx1 = p1.index;
	//	int idx2 = p2.index;
	//	
	//	//idx1, idx2��˳����Ҫ
	//	dt::Edge_Node tmp_en = dt::Edge_Node{ idx1, idx2 };

	//	if (edge_map.find(tmp_en) == edge_map.end()) {//û�ҵ�
	//													//std::cout << "not found this edge in edge_map!" << std::endl;

	//		std::vector<TriangleType> map_tmp_triangles;
	//		map_tmp_triangles.clear();
	//		map_tmp_triangles.push_back(idx);
	//		edge_map[Edge_Node{ idx1, idx2 }] = map_tmp_triangles;
	//	}
	//	else {
	//		//std::cout << "found this edge in edge_map!" << std::endl;
	//		edge_map[Edge_Node{ idx1, idx2 }].push_back(idx);
	//	}
	//	return true;
	//	
	//}

	//template<typename T>
	//void tearing(std::vector<dt::Edge<T>> tear_edges, std::vector<dt::Vector2<T>> tear_pts, dt::Delaunay<T> &dela) {//���ѵ�Ҳ�ö��壻
	//	std::map<int, bool> if_tear_pt;//tearing_pts;��i�������Ƿ���tearing_point;
	//	//��ʼ��false
	//	
	//	for (int i = 0; i < tear_pts.size(); i++)
	//		if_tear_pt[tear_pts[i].index] = true;

	//	//�ȱ��ordinary_triangles;
	//	auto triangles = dela._triangles;
	//	//std::map<int, int> triangle_type;//0Ϊordinary, 1Ϊunchanged�� 2Ϊchanged;
	//	////0Ϊord;
	//	////��ʼ��Ϊ-1;


	//	//for (int i = 0; i < triangles.size(); i++)
	//	//{
	//	//	int pt_a = triangles[i].a->index;
	//	//	int pt_b = triangles[i].b->index;
	//	//	int pt_c = triangles[i].c->index;
	//	//	if (if_tear_pt[pt_a] == false && if_tear_pt[pt_b] == false && if_tear_pt[pt_c] == false) {
	//	//		triangle_type[i] = 0;
	//	//	}
	//	//}

	//	////
	//	std::vector<bool> have_appeared_pre_graph;//��¼������ǰ����ͨͼ���Ƿ���ֹ��������жϺ����������Ƿ�����¶���
	//	std::set<int> have_appeared_cur_graph;//��¼�ѷ��ʵĶ����ڵ�ǰ��ͨͼ���Ƿ���ֹ�;
	//	
	//	std::vector<bool> have_vis_tri;//��¼�����η��ʼ�¼������dfs
	//	for (int i = 0; i < dela._vertices.size(); i++) {
	//		have_appeared_pre_graph.push_back(false);
	//	}
	//	for (int i = 0; i < triangles.size(); i++) {
	//		have_vis_tri.push_back(false);
	//	}

	//	std::vector<TriangleType >> &edge_map;//��¼������Ӧ��������;
	//	for (int i = 0; i < triangles.size(); i++) {
	//		auto p1 = *(triangles.a);
	//		auto p2 = *(triangles.b);
	//		auto p3 = *(triangles.c);

	//		add_edge_map(i, p1, p2, edge_map);
	//		add_edge_map(i, p2, p3, edge_map);
	//		add_edge_map(i, p1, p3, edge_map);
	//	}


	//	int total = dela._vertices.size();
	//	for (int i = 0; i < triangles.size(); i++) {
	//		if(!have_vis_tri[i]){
	//			//triangles, idx, total, have_vis_tri, have_appeared_pre_graph, have_appeared_cur_graph, map_tri2new_pt
	//			vis_pts = dfs(triangles, i, total, have_vis_tri, have_appeared_pre_graph, have_appeared_cur_graph, edge_map, if_tear_pt, map_tri2new_pt);//have_appeared;//��dfs�з��ʹ��Ķ����¼����
	//			//TODO����ͨͼ�ж����¼��have_appeared_pre_graph��

	//		}
	//		have_appeared_cur_graph.clear();
	//	}
	//}
	template<typename T>
	bool judge_if_common_edge(int tri_index1, int tri_index2, std::vector<dt::Triangle<T>> triangles, std::vector<int> & common_edge, dt::Vector2<T> tear_pt) {
		dt::Triangle<T> tri1 = triangles[tri_index1];
		dt::Triangle<T> tri2 = triangles[tri_index2];
		
		std::set<int> s1;
		std::set<int> s2;
		std::set<int> s_result;
		
		s1.insert(tri1.a->index);
		s1.insert(tri1.b->index);
		s1.insert(tri1.c->index);
		
		s2.insert(tri2.a->index);
		s2.insert(tri2.b->index);
		s2.insert(tri2.c->index);

		std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(s_result, s_result.begin()));
		if (s_result.size() == 2) {
			for (auto it = s_result.begin(); it != s_result.end(); it++) {
				common_edge.push_back(*it);
			}

			return true;
		}
		else {
			return false;
		}
	}

	template<typename T>
     void tear_func(std::vector<dt::Edge<T>> tear_edges, std::vector<dt::Vector2<T>> tear_pts, dt::Delaunay<T> &dela, std::map<int, int>&map_tri2new_pt, int & r_total, std::vector<dt::Vector2<T>> & new_pts)
	 //void tear_func(std::vector<dt::Edge<double>> tear_edges, std::vector<dt::Vector2<double>> tear_pts, dt::Delaunay<double> &dela, std::map<int, int>&map_tri2new_pt, int & r_total)

	{//���ѵ�Ҳ�ö��壻
		std::map<int, bool> if_tear_pt;//tearing_pts;��i�������Ƿ���tearing_point;
									   //��ʼ��false
		for (int i = 0; i < tear_pts.size(); i++)
			if_tear_pt[tear_pts[i].index] = true;

		std::set<dt::Edge_Node> tear_edges_set;
		for (auto te: tear_edges) {
			dt::Edge_Node tmp_en{ te.v->index, te.w->index };
			tear_edges_set.insert(tmp_en);
		}
		//�ȱ��ordinary_triangles;
		auto triangles = dela._triangles;
		for (auto pt : dela._vertices) {//����һ��
			new_pts.push_back(pt);
		}

		std::map<int, std::vector<int>> map_pt2tri;
		for (int i = 0; i < triangles.size(); i++) {
			//map_pt2tri[i];
			int p1 = triangles[i].a->index;
			int p2 = triangles[i].b->index;
			int p3 = triangles[i].c->index;
			
			if (map_pt2tri.find(p1) == map_pt2tri.end()) {
				std::vector<int> tri_tmp;
				tri_tmp.clear();
				tri_tmp.push_back(i);
				map_pt2tri[p1] = tri_tmp;
			}
			else {
				map_pt2tri[p1].push_back(i);
			}

			if (map_pt2tri.find(p2) == map_pt2tri.end()) {
				std::vector<int> tri_tmp;
				tri_tmp.clear();
				tri_tmp.push_back(i);
				map_pt2tri[p2] = tri_tmp;
			}
			else {
				map_pt2tri[p2].push_back(i);
			}

			if (map_pt2tri.find(p3) == map_pt2tri.end()) {
				std::vector<int> tri_tmp;
				tri_tmp.clear();
				tri_tmp.push_back(i);
				map_pt2tri[p3] = tri_tmp;
			}
			else {
				map_pt2tri[p3].push_back(i);
			}

		}
		for (int m = 0; m < tear_pts.size(); m++) {
			std::vector<int> tri_indices=map_pt2tri[tear_pts[m].index];
			//ð������
			for (int i = 0; i<tri_indices.size() - 1; i++)
				for (int j = 0; j < tri_indices.size() - i - 1; ++j) {
					if (TriSortCmp(tri_indices[j], tri_indices[j + 1], triangles, tear_pts[m])) {
						int tmp = tri_indices[j];
						tri_indices[j] = tri_indices[j + 1];
						tri_indices[j + 1] = tmp;
					}
				}

			map_pt2tri[tear_pts[m].index] = tri_indices;
		}
		
		//total;��Ҫ���ؽ��
		int total = dela._total;
		for (int i = 0; i < tear_pts.size(); i++) {
			int cur_index = tear_pts[i].index;
			std::vector<int> tri_indices = map_pt2tri[cur_index];
			//�ҵ��õ�������ص�������
			for (int j = 1; j < tri_indices.size();j++) {//����360���жϣ�
				//�ж��Ƿ���Ҫ���Ѷ���tear_pt;
				//������
				std::vector<int> common_edge;
				bool if_common_edge = judge_if_common_edge(tri_indices[j-1], tri_indices[j], triangles, common_edge, tear_pts[i]);
				if (!if_common_edge) {
					cur_index = total;
					new_pts.push_back(dt::Vector2<double>{tear_pts[i].x, tear_pts[i].y, cur_index});
					total++;
				}
				else {
					//��ͬ���Ƿ�Ϊ˺�ѱ�
					if (tear_edges_set.find(Edge_Node{common_edge[0], common_edge[1]}) != tear_edges_set.end()) {
						cur_index = total;
						new_pts.push_back(dt::Vector2<double>{tear_pts[i].x, tear_pts[i].y, cur_index});
						total++;
						//˺�ѱ߹���
						//��˺�ѱ߹���
					}
				}
	
				if (cur_index != tear_pts[i].index) {
					int pt_1 = triangles[tri_indices[j]].a->index;
					int pt_2 = triangles[tri_indices[j]].b->index;
					int pt_3 = triangles[tri_indices[j]].c->index;
					

					//��¼�������ζ�Ӧ��tear_pt����Ӧ����pt;
					if (pt_1 == tear_pts[i].index) {
						map_tri2new_pt[3 * tri_indices[j] + 0] = cur_index;
					}
					else if (pt_2 == tear_pts[i].index) {
						map_tri2new_pt[3 * tri_indices[j] + 1] = cur_index;
					}
					else if (pt_3 == tear_pts[i].index) {
						map_tri2new_pt[3 * tri_indices[j] + 2] = cur_index;
					}
				}
			}
		}

		r_total = total;
	}
}