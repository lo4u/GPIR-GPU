#include "occ.h"

// issue: inplement the data structure of hypergraph and the peeling algorithm of it
// We use a edge-vertex matrix and a pointer array which consists of 
// a node_id of vertex and a hash table of edges which the vertex is touched.
// Edge is consists of two node_id.
// Perfect data structure !!!

struct Edge
{
    size_t edge_id;
    size_t vertex[2]; // The vertex which it touches

    Edge(size_t e, size_t v1, size_t v2): edge_id(e)
    {
        vertex[0] = v1;
        vertex[1] = v2;
    }

    Edge() 
    {

    }
};

struct Vertex
{
    size_t vertex_id; // Preserve the vertex_id of this vertex
    std::unordered_map<size_t, Edge> edge_table; // Preserve the edge this vertex is touched
    
    Vertex(size_t v): vertex_id(v)
    {

    }

    Vertex()
    {

    }
};

struct VertexOriented
{
    size_t edge_id;
    size_t vertex_id;

    VertexOriented()
    {

    }

    VertexOriented(size_t e, size_t v): edge_id(e), vertex_id(v)
    {

    }
};

bool GCTHelper(
    std::unordered_map<size_t, size_t>& elements, // 输出结果：保存 position -> [current idx] 
    std::unordered_map<size_t, std::vector<size_t>>& positions, // 输入：保存 IDX -> [position 1, ..., position d] 的映射
    const size_t& idx, // 输入：待插入的键
    size_t attempt, // 尝试次数
    std::mt19937_64 rng // 随机性
);

size_t hash_and_mod_with_mt19937(
    size_t id, 
    size_t nonce, 
    const std::vector<unsigned char>& data, 
    size_t modulus
);

std::vector<bool> peel_the_hypergraph(
    std::stack<VertexOriented> &Stack,
    std::vector<Vertex> HyperGraph,
    std::vector<std::vector<size_t>> NeiborMatrix
);

void print_neibormatrix(
    std::vector<std::vector<size_t>> NeiborMatrix
);

void print_neibormatrix_without_index(
    std::vector<std::vector<size_t>> NeiborMatrix
);

void print_flag_list(
    std::vector<bool> flagList
);

void print_stack(
    std::stack<VertexOriented> Stack
);

void print_hypergraph(
    std::vector<Vertex> HyperGraph
);

void HyperGraphSetup(
    std::vector<Vertex> &HyperGraph
);

std::vector<std::vector<size_t>> choice_rows_in_neiborMatirx(
    std::vector<bool> rows_flag,
    std::vector<std::vector<size_t>> M
);

std::vector<MatPoly> choice_position_in_real_pt(
    std::vector<bool> rows_flag,
    std::vector<MatPoly> real_pt
);

void print_some_vars(
    std::vector<std::vector<size_t>> NeiborMatrix,
    std::vector<Vertex> HyperGraph,
    std::vector<bool> flagList,
    std::stack<VertexOriented> Stack,
    std::vector<std::vector<size_t>> core2
);

void print_solutions(
    std::vector<MatPoly> solutions
);

std::vector<MatPoly> get_whole_solution(
    std::vector<MatPoly> solution,
    std::stack<VertexOriented> Stack,
    std::vector<MatPoly> real_pt,
    std::unordered_map<size_t, std::vector<size_t>> positions,
    std::vector<size_t> realIndexes
);
