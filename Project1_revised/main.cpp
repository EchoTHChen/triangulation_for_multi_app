#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>

#include <SFML/Graphics.hpp>

#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"
#include "tear.h"
#include"intersection.h"


//namespace dt {
//	template<typename T>
//	void tear_func(std::vector<dt::Edge<T>> tear_edges, std::vector<dt::Vector2<T>> tear_pts, dt::Delaunay<T> &dela, std::map<int, int>&map_tri2new_pt, int & r_total);
//}
//#include "Polygon_cal.h"
//
//const std::vector<dt::Triangle<double>>  generate_thing(points) {
//
//}

//template<typename T>

void print_pts_index(std::vector<dt::Triangle<double>> triangles) {
	std::cout << "triangles:" << std::endl;

	for (int i = 0; i < triangles.size(); i++) {
		std::cout << triangles[i].a->index << std::endl;
		std::cout << triangles[i].b->index << std::endl;
		std::cout << triangles[i].c->index << std::endl;

	}
}


int main(int argc, char * argv[])
{
	int numberPoints = 6;
//	if (argc > 1)
//	{
//		numberPoints = atoi(argv[1]);
//	}
	// SFML window
	sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay triangulation");
	window.setFramerateLimit(1);

//	const auto start = std::chrono::high_resolution_clock::now();
//	//逆时针
	//生成背景的mesh triangulation:

	std::vector<double> arr_x_bg = { 720, 100, 100, 720 };
	std::vector<double> arr_y_bg = { 30,  30,  300, 300 };
	std::vector<dt::Vector2<double>> points_bg;
	//int numberPoints_bg = 4;

	for (int i = 0; i < arr_x_bg.size(); ++i) {
		points_bg.push_back(dt::Vector2<double>{arr_x_bg[i], arr_y_bg[i], i});
	}
	dt::Delaunay<double> triangulation_bg;

	//前景
	dt::Delaunay<double> triangulation;
	std::vector<dt::Vector2<double>> points;
	std::cout << "Generating " << numberPoints << " points" << std::endl;
	//逆时针
	std::vector<double> arr_x = { 600, 200, 350, 200, 600, 450 };//W
	std::vector<double> arr_y = { 60 , 60, 240, 400,  400, 240 };//H
	for (int i = 0; i < numberPoints; ++i) {
		points.push_back(dt::Vector2<double>{arr_x[i], arr_y[i], i});
	}
	int option = 3;
	if(option == 0){//ear-clipping三角剖分
		//三角剖分
		std::vector<dt::Triangle<double>> result_triangles;
		const std::vector<dt::Triangle<double>> triangles = triangulation.triangulate(points);
		//const std::vector<dt::Triangle<double>> triangles_bg= triangulation_bg.triangulate(points_bg);

		//(std::vector<dt::Vector2<T>> &points, std::vector<dt::Triangle<T>>&result_triangles) {
		std::cout << "triangles.size():" << triangles.size() << std::endl;
		std::vector<dt::Edge<double>> edges;
		for (const auto t : triangles)
		{
			edges.push_back(dt::Edge<double>{*t.a, *t.b});
			edges.push_back(dt::Edge<double>{*t.b, *t.c});
			edges.push_back(dt::Edge<double>{*t.c, *t.a});
		}

		//for (const auto t : triangles_bg)
		//{
		//	edges.push_back(dt::Edge<double>{*t.a, *t.b});
		//	edges.push_back(dt::Edge<double>{*t.b, *t.c});
		//	edges.push_back(dt::Edge<double>{*t.c, *t.a});
		//}

		for (auto edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) {
			auto e = *edge_iter;
	
			const std::array<sf::Vertex, 2> line{ {
					sf::Vertex(sf::Vector2f(
						static_cast<float>(e.v->x + 2.),
						static_cast<float>(e.v->y + 2.))),
				sf::Vertex(sf::Vector2f(
					static_cast<float>(e.w->x + 2.),
					static_cast<float>(e.w->y + 2.))),
				} };
			window.draw(line.data(), 2, sf::Lines);
		}

		window.display();

		while (window.isOpen())
		{
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();
			}
		}
	
	}
	else if (option == 1) {//多边形相交
		//三角剖分
		const std::vector<dt::Triangle<double>> triangles = triangulation.triangulate(points);
		const std::vector<dt::Triangle<double>> triangles_bg = triangulation_bg.triangulate(points_bg);

		std::vector<dt::Triangle<double>> triangles_inter;
		std::vector<dt::Vector2<double>> new_pts;
		/*dt::Intersection<double> inter;*/
		intersection(triangulation, triangulation_bg, triangles_inter, new_pts);
		std::cout << "triangles_inter.size():" << triangles_inter.size() << std::endl;

		std::vector<dt::Edge<double>> edges;
		for (const auto t : triangles_inter)
		{
			edges.push_back(dt::Edge<double>{*t.a, *t.b});
			edges.push_back(dt::Edge<double>{*t.b, *t.c});
			edges.push_back(dt::Edge<double>{*t.c, *t.a});
		}

		for (auto edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) {
			auto e = *edge_iter;

			const std::array<sf::Vertex, 2> line{ {
					sf::Vertex(sf::Vector2f(
						static_cast<float>(e.v->x + 2.),
						static_cast<float>(e.v->y + 2.))),
				sf::Vertex(sf::Vector2f(
					static_cast<float>(e.w->x + 2.),
					static_cast<float>(e.w->y + 2.))),
				} };
			window.draw(line.data(), 2, sf::Lines);
		}

		window.display();

		while (window.isOpen())
		{
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();
			}
		}



	}
	else if (option == 2) {//网格细分
		std::cout << "endl" << std::endl;
		triangulation.triangulate(points);
		const std::vector<dt::Triangle<double>> triangles = triangulation.subdivision(3);
		std::cout << "run" << std::endl;

		std::vector<dt::Edge<double>> edges;
		for (const auto t : triangles)
		{
			edges.push_back(dt::Edge<double>{*t.a, *t.b});
			edges.push_back(dt::Edge<double>{*t.b, *t.c});
			edges.push_back(dt::Edge<double>{*t.c, *t.a});
		}

		for (auto edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) {
			auto e = *edge_iter;

			const std::array<sf::Vertex, 2> line{ {
					sf::Vertex(sf::Vector2f(
						static_cast<float>(e.v->x + 2.),
						static_cast<float>(e.v->y + 2.))),
				sf::Vertex(sf::Vector2f(
					static_cast<float>(e.w->x + 2.),
					static_cast<float>(e.w->y + 2.))),
				} };
			window.draw(line.data(), 2, sf::Lines);
		}

		window.display();

		while (window.isOpen())
		{
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();
			}
		}

	}
	else if(option==3){//撕裂边受限情况下的三角剖分
		std::vector<dt::Edge<double>> line1;
		std::vector<dt::Edge<double>> line2;
		std::vector<std::vector<dt::Edge<double>>> lines;
		std::vector<double> test_case1_x = {400, 375, 425, 425};
		std::vector<double> test_case1_y = {160,  240, 240, 320};

		std::vector<double> test_case2_x = { 375, 400, 400 };
		std::vector<double> test_case2_y = { 160,  240, 400 };
		std::vector<dt::Vector2<double>> pts_vec1;
		std::vector<dt::Vector2<double>> pts_vec2;
		for (int i = 0; i < test_case1_x.size(); i++) {
			pts_vec1.push_back(dt::Vector2<double>{test_case1_x[i], test_case1_y[i]});

		}
		for (int i = 0; i < test_case2_x.size(); i++) {
			pts_vec2.push_back(dt::Vector2<double>{test_case2_x[i], test_case2_y[i]});
		}
		
		for (int i = 1; i < pts_vec1.size(); i++) {
			line1.push_back(dt::Edge<double>{pts_vec1[i-1], pts_vec1[i]});
		}

		for (int i = 1; i < pts_vec2.size(); i++) {
			line2.push_back(dt::Edge<double>{pts_vec2[i - 1], pts_vec2[i]});
		}


		//const std::vector<dt::Triangle<double>> triangles  = 
		triangulation.triangulate(points);
		//lines.push_back(line1);
		lines.push_back(line2);
		const std::vector<dt::Triangle<double>> triangles = triangulation.polygon_lines(lines);
		std::cout << "run" << std::endl;
		std::vector<dt::Edge<double>> edges;
		/*
		
		for (int  i = 0;i<line2.size(); i++)
			edges.push_back(line2[i]);
		int n = points.size();
		for (int i = 0; i < points.size(); i++)
			edges.push_back(dt::Edge<double>{points[i], points[(i + 1) % n]});*/



		std::cout << "triangle.size();" << triangles.size() << std::endl;
		for (int i = 0; i < triangles.size(); i++) {
			std::cout << "triangles[i]:" << i << std::endl;
			std::cout << "triangles[i].a:" << triangles[i].a->x << "," << triangles[i].a->y << std::endl;
			std::cout << "triangles[i].b:" << triangles[i].b->x << "," << triangles[i].b->y << std::endl;
			std::cout << "triangles[i].c:" << triangles[i].c->x << "," << triangles[i].c->y << std::endl;
		}
		for (const auto t : triangles)
		{
			edges.push_back(dt::Edge<double>{*t.a, *t.b});
			edges.push_back(dt::Edge<double>{*t.b, *t.c});
			edges.push_back(dt::Edge<double>{*t.c, *t.a});
		}

		for (auto edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) {
			auto e = *edge_iter;

			const std::array<sf::Vertex, 2> line{ {
					sf::Vertex(sf::Vector2f(
						static_cast<float>(e.v->x + 2.),
						static_cast<float>(e.v->y + 2.))),
				sf::Vertex(sf::Vector2f(
					static_cast<float>(e.w->x + 2.),
					static_cast<float>(e.w->y + 2.))),
				} };
			window.draw(line.data(), 2, sf::Lines);
		}

		window.display();

		while (window.isOpen())
		{
			sf::Event event;
			while (window.pollEvent(event))
			{
				if (event.type == sf::Event::Closed)
					window.close();
			}
		}

	}
		
		
	//
	//	// Transform each points of each vector as a rectangle
	//	for (const auto p : points) {
	//		sf::RectangleShape s{ sf::Vector2f(4, 4) };
	//		s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
	//		window.draw(s);
	//	}
	//
	//	std::vector<std::array<sf::Vertex, 2> > lines;

//	if (option == 1) {	
//
//		
//
//
//		//求交集
//		std::vector<dt::Triangle<double>> triangles_inter;
//		std::vector<dt::Vector2<double>> new_pts;
//		intersection(triangulation, triangulation_bg, triangles_inter, new_pts);
//
//		//print_pts_index(triangles_inter);
//	
//	
//		std::vector<dt::Edge<double>> edges_inter;// = triangulation_inter.getEdges();
//
//		dt::Delaunay<double> triangulation_inter;//进行loop_subdivision 和设置 _edges，original_vertices;
//		triangulation_inter.setTriangles(triangles_inter, new_pts);
//
//
//		const std::vector<dt::Triangle<double>> new_triangles_inter = triangulation_inter.subdivision(5);
//
//		for (const auto t : new_triangles_inter)
//		{
//			edges_inter.push_back(dt::Edge<double>{*t.a, *t.b});
//			edges_inter.push_back(dt::Edge<double>{*t.b, *t.c});
//			edges_inter.push_back(dt::Edge<double>{*t.c, *t.a});
//		}
//	
//		for (auto edge_iter = edges_inter.begin(); edge_iter != edges_inter.end(); edge_iter++) {
//			auto e = *edge_iter;
//
//			const std::array<sf::Vertex, 2> line{ {
//					sf::Vertex(sf::Vector2f(
//						static_cast<float>(e.v->x + 2.),
//						static_cast<float>(e.v->y + 2.))),
//				sf::Vertex(sf::Vector2f(
//					static_cast<float>(e.w->x + 2.),
//					static_cast<float>(e.w->y + 2.))),
//				} };
//			window.draw(line.data(), 2, sf::Lines);
//		}
//	}
//	else if(option == 2){//subdivision
//
//		////网格细分
//		//const std::vector<dt::Triangle<double>> new_triangles_bg = triangulation_bg.loop_subdivision(5);
//		//const std::vector<dt::Triangle<double>> new_triangles = triangulation.loop_subdivision(5);
//
//
//
//
//		
//
//		//std::default_random_engine eng(std::random_device{}());
//		/*std::uniform_real_distribution<double> dist_w(0, 800);
//		std::uniform_real_distribution<double> dist_h(0, 600);
//		*/
//
//		const auto end = std::chrono::high_resolution_clock::now();
//
//		const std::chrono::duration<double> diff = end - start;
//
//
//		//std::vector<std::array<sf::Vertex, 2> > lines_bg;
//		//triangulation_bg.loop_subdivision(5);
//
//		std::vector<dt::Edge<double>> edges_bg = triangulation_bg.getEdges();
//
//		for (auto edge_iter = edges_bg.begin(); edge_iter != edges_bg.end(); edge_iter++) {
//			auto e = *edge_iter;
//
//			const std::array<sf::Vertex, 2> line{ {
//					sf::Vertex(sf::Vector2f(
//						static_cast<float>(e.v->x + 2.),
//						static_cast<float>(e.v->y + 2.))),
//				sf::Vertex(sf::Vector2f(
//					static_cast<float>(e.w->x + 2.),
//					static_cast<float>(e.w->y + 2.))),
//				} };
//			window.draw(line.data(), 2, sf::Lines);
//		
//		}
//
//		
//		std::cout << triangles.size() << " triangles generated in " << diff.count()
//			<< "s\n";
//		triangulation.subdivision(5);
//
//		const std::vector<dt::Edge<double>> edges = triangulation.getEdges();
//		for (const auto &e : edges) {
//			const std::array<sf::Vertex, 2> line{ {
//					sf::Vertex(sf::Vector2f(
//						static_cast<float>(e.v->x + 2.),
//						static_cast<float>(e.v->y + 2.))),
//				sf::Vertex(sf::Vector2f(
//					static_cast<float>(e.w->x + 2.),
//					static_cast<float>(e.w->y + 2.))),
//				} };
//			window.draw(line.data(), 2, sf::Lines);
//		}
//
//
//	}
//	else  if(option==3){//撕裂
//		//给定撕裂点和撕裂边和待撕裂网格
//		//求每个三角形新的对应顶点
//		//将新的顶点试着调整坐标，看看效果；
//		
//		//std::vector<dt::Edge<double>> tear_edges, std::vector<dt::Vector2<double>> tear_pts, dt::Delaunay<double> &dela, std::map<int, int>&map_tri2new_pt, int & r_total
//		int total;
//		std::map<int, int> map_tri2new_pt;
//		std::vector<dt::Edge<double>> tear_edges;
//
//		std::vector<dt::Vector2<double>> tear_pts;
//		std::vector<double> tear_pts_pos_x{200, 450, 350};
//		std::vector<double> tear_pts_pos_y{80, 240, 240};
//
//		for (int i = 0; i < tear_pts_pos_x.size(); i++) {
//			bool flag = false;
//			for (auto v : triangulation._vertices) {
//				if (dt::almost_equal(v.x, tear_pts_pos_x[i]) && dt::almost_equal(v.y, tear_pts_pos_y[i])) {
//					tear_pts.push_back(dt::Vector2<double>{tear_pts_pos_x[i],tear_pts_pos_y[i], v.index});
//					flag = true;
//					break;
//				}
//				
//			}
//			if (!flag) {
//				std::cout << "x,y does not exist!" << tear_pts_pos_x[i] << "," << tear_pts_pos_y[i] << std::endl;
//				return 0;
//			}
//		}
//		//dt::Vector2<double> pt1{200, 80,};
//		std::vector<double> tear_pts_pos_s_x{ 200, 450 };
//		std::vector<double> tear_pts_pos_s_y{ 80, 240 };
//		std::vector<double> tear_pts_pos_d_x{ 450, 350 };
//		std::vector<double> tear_pts_pos_d_y{ 240, 240 };
//
//
//		for (int i = 0; i < tear_pts_pos_s_x.size(); i++) {
//			bool flag1 = false, flag2=false;
//			dt::Vector2<double>s, d;
//			for (auto v : triangulation._vertices) {
//
//				if (dt::almost_equal(v.x, tear_pts_pos_s_x[i]) && dt::almost_equal(v.y, tear_pts_pos_s_y[i])) {
//					s.x = tear_pts_pos_s_x[i];
//					s.y = tear_pts_pos_s_y[i];
//					s.index = v.index;
//					flag1 = true;
//				}
//
//				if (dt::almost_equal(v.x, tear_pts_pos_d_x[i]) && dt::almost_equal(v.y, tear_pts_pos_d_y[i])) {
//					d.x = tear_pts_pos_d_x[i];
//					d.y = tear_pts_pos_d_y[i];
//					d.index = v.index;
//					flag2 = true;
//				}
//			}
//			if((!flag1)|| (!flag2)){
//				std::cout << "tear_edge dose not exist:No." << i << std::endl;
//				if ((!flag1)) {
//					std::cout << "x,y does not exist!" << tear_pts_pos_s_x[i]<<","<<tear_pts_pos_s_y[i] << std::endl;
//					
//				}
//				if ((!flag2)) {
//					std::cout << "x,y does not exist!" << tear_pts_pos_d_x[i] << "," << tear_pts_pos_d_y[i] <<std::endl;
//					
//				}
//				return 0;
//			}
//			else {
//				tear_edges.push_back(dt::Edge<double>{s, d});
//			}
//		}
//		
//
//		//检查合法性check_legal TODO?
//		//template<typename T>
//		//void tear_func(std::vector<dt::Edge<T>> tear_edges, std::vector<dt::Vector2<T>> tear_pts, dt::Delaunay<T> &dela, std::map<int, int>&map_tri2new_pt, int & r_total)
//		std::vector < dt::Vector2<double>>  new_pts;
//		tear_func(tear_edges, tear_pts, triangulation, map_tri2new_pt, total, new_pts);
//		//std::vector<dt::Triangle<double>> new_triangles;
//
//		std::vector<dt::Triangle<double>> new_triangles = triangles;
//		for (auto it = map_tri2new_pt.begin(); it!=map_tri2new_pt.end();it++) {
//			int tmp_tri_idx = (it->first) / 3;
//			int internal_idx = (it->first) % 3;
//			int new_pt_idx = it->second;
//			std::cout << "triangle index:" << tmp_tri_idx << std::endl;
//			std::cout << "triangle:" << std::endl;
//			auto p1 = triangulation._triangles[tmp_tri_idx].a;
//			auto p2 = triangulation._triangles[tmp_tri_idx].b;
//			auto p3 = triangulation._triangles[tmp_tri_idx].c;
//			std::cout << "A" << p1->x << "," << p1->y << std::endl;
//			std::cout << "B" << p2->x << "," << p2->y << std::endl;
//			std::cout << "C" << p3->x << "," << p3->y << std::endl;
//			if (internal_idx == 0) {
//				std::cout << "A's index became " << new_pt_idx << std::endl;
//				new_triangles[tmp_tri_idx].a->index = new_pt_idx;
//			}
//			else if (internal_idx == 1) {
//				std::cout << "B's index became " << new_pt_idx << std::endl;
//				new_triangles[tmp_tri_idx].b->index = new_pt_idx;
//
//			}
//			else if(internal_idx == 2) {
//				std::cout << "C's index became " << new_pt_idx << std::endl;
//				new_triangles[tmp_tri_idx].c->index = new_pt_idx;
//			}
//			//std::cout << "internal index:" << tmp_tri_idx << std::endl;
//
//		}
//		//TODO修改坐标，产生撕裂效果
//		//pt1， x+20, y-20;
//		new_pts[1].x += 20;
//		new_pts[1].y -= 20;
//		new_pts[2].x -= 10;
//		new_pts[2].x -= 30;
//		//pt2, x-10, y-30;
//		new_pts[5].x += 50;
//		new_pts[5].y += 50;
//		//pt5, x+50, y+50;
//		//for (auto it = map_tri2new_pt.begin(); it != map_tri2new_pt.end(); it++) {
//		//	int tmp_tri_idx = (it->first) / 3;
//		//	int internal_idx = (it->first) % 3;
//		//	int new_pt_idx = it->second;
//		//	std::cout << "triangle index:" << tmp_tri_idx << std::endl;
//		//	std::cout << "triangle:" << std::endl;
//		//	auto p1 = triangulation._triangles[tmp_tri_idx].a;
//		//	auto p2 = triangulation._triangles[tmp_tri_idx].b;
//		//	auto p3 = triangulation._triangles[tmp_tri_idx].c;
//		//	//std::cout << "A" << p1->x << "," << p1->y << std::endl;
//		//	//std::cout << "B" << p2->x << "," << p2->y << std::endl;
//		//	//std::cout << "C" << p3->x << "," << p3->y << std::endl;
//		//	if (internal_idx == 0) {
//		//		std::cout << "A's index became " << new_pt_idx << std::endl;
//		//		new_triangles[tmp_tri_idx].a->x  = new_pts[new_pt_idx].x;
//		//		new_triangles[tmp_tri_idx].a->y = new_pts[new_pt_idx].y;
//
//		//	}
//		//	else if (internal_idx == 1) {
//		//		std::cout << "B's index became " << new_pt_idx << std::endl;
//		//		new_triangles[tmp_tri_idx].b->x = new_pts[new_pt_idx].x;
//		//		new_triangles[tmp_tri_idx].b->y = new_pts[new_pt_idx].y;
//
//		//	}
//		//	else if (internal_idx == 2) {
//		//		std::cout << "C's index became " << new_pt_idx << std::endl;
//		//		new_triangles[tmp_tri_idx].c->x = new_pts[new_pt_idx].x;
//		//		new_triangles[tmp_tri_idx].c->y = new_pts[new_pt_idx].y;
//		//	}
//		//	//std::cout << "internal index:" << tmp_tri_idx << std::endl;
//
//		//}
//	////	for (auto it : new_triangles) {
//	////		it->a->x = new_pts[it->a->index].x;
//	////		it->a->y = new_pts[it->a->index].x;
//
//
//	////	}
//	////	triangulation.setTriangles(new_triangles, new_pts);
//	////	std::vector<dt::Edge<double>> edges;// = triangulation.getEdges();
//
//	/////*	for (const auto &e : edges) {
//	////		const std::array<sf::Vertex, 2> line{ {
//	////				sf::Vertex(sf::Vector2f(
//	////					static_cast<float>(e.v->x + 2.),
//	////					static_cast<float>(e.v->y + 2.))),
//	////			sf::Vertex(sf::Vector2f(
//	////				static_cast<float>(e.w->x + 2.),
//	////				static_cast<float>(e.w->y + 2.))),
//	////			} };
//	////		window.draw(line.data(), 2, sf::Lines);
//	////	}*/
//	////	
//	////	for (const auto t : new_triangles)
//	////	{
//	////		edges.push_back(dt::Edge<double>{*t.a, *t.b});
//	////		edges.push_back(dt::Edge<double>{*t.b, *t.c});
//	////		edges.push_back(dt::Edge<double>{*t.c, *t.a});
//	////	}
//
//	////	for (auto edge_iter = edges.begin(); edge_iter != edges.end(); edge_iter++) {
//	////		auto e = *edge_iter;
//
//	////		const std::array<sf::Vertex, 2> line{ {
//	////				sf::Vertex(sf::Vector2f(
//	////					static_cast<float>(e.v->x + 2.),
//	////					static_cast<float>(e.v->y + 2.))),
//	////			sf::Vertex(sf::Vector2f(
//	////				static_cast<float>(e.w->x + 2.),
//	////				static_cast<float>(e.w->y + 2.))),
//	////			} };
//	////		window.draw(line.data(), 2, sf::Lines);
//	////	}
//	}
//
//
//
//	
//	
////
////	// Transform each points of each vector as a rectangle
////	for (const auto p : points) {
////		sf::RectangleShape s{ sf::Vector2f(4, 4) };
////		s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
////		window.draw(s);
////	}
////
////	std::vector<std::array<sf::Vertex, 2> > lines;
//	
////
//	window.display();
//
//	while (window.isOpen())
//	{
//		sf::Event event;
//		while (window.pollEvent(event))
//		{
//			if (event.type == sf::Event::Closed)
//				window.close();
//		}
//	}
//
//
//
//	return 0;
	system("pause");
	return 0;
}
