#pragma once
#ifndef H_TEAR
#define H_TEAR
#include"delaunay.h"
namespace dt{
	template<typename T>
	bool TriSortCmp(int pt_1_index, int pt_2_index, const std::vector<dt::Triangle<T>> &triangles, dt::Vector2<T> tear_pt);

	template<typename T>
	bool judge_if_common_edge(int tri_index1, int tri_index2, std::vector<dt::Triangle<T>> triangles, std::vector<int> & common_edge, dt::Vector2<T> tear_pt);

	
	template<typename T>
	void tear_func(std::vector<dt::Edge<T>> tear_edges, std::vector<dt::Vector2<T>> tear_pts, dt::Delaunay<T> &dela, std::map<int, int>&map_tri2new_pt, int & r_total, std::vector<dt::Vector2<T>> & new_pts);
	//ÉùÃ÷
	template void tear_func(std::vector<dt::Edge<double>> tear_edges, std::vector<dt::Vector2<double>> tear_pts, dt::Delaunay<double> &dela, std::map<int, int>&map_tri2new_pt, int & r_total, std::vector<dt::Vector2<double> >& new_pts);

	//void tear_func(std::vector<dt::Edge<double>> tear_edges, std::vector<dt::Vector2<double>> tear_pts, dt::Delaunay<double> &dela, std::map<int, int>&map_tri2new_pt, int & r_total);
	
}
#endif