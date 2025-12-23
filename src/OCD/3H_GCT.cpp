#include "3H_GCT.h"
const size_t MAX_ATTEMPTS_NUM = 1000;
// 插入元素
bool GCTHelper(
    std::unordered_map<size_t, size_t>& elements, // 输出结果：保存 position -> [current idx] 
    std::unordered_map<size_t, std::vector<size_t>>& positions, // 输入：保存 IDX -> [position 1, ..., position d] 的映射
    const size_t& idx, // 输入：待插入的键
    size_t attempt = 0, // 尝试次数
    std::mt19937_64 rng = std::mt19937_64() // 随机性
) {
    if (attempt >= MAX_ATTEMPTS_NUM) {
        return false;
    }

    // Case 1: 检查是否有空桶可以插入
    for (auto position : positions[idx]) {
        if (elements.find(position) == elements.end()) {
            elements[position] = idx;
            return true;
        }
    }

    // Case 2: 所有可能的桶都满了，需要重新定位一个现有元素
    size_t index = std::uniform_int_distribution<size_t>(0, positions[idx].size() - 1)(rng);
    size_t chosen_position = positions[idx][index];

    // 插入新键，并获取之前插入的键
    auto old_idx = elements[chosen_position];
    elements[chosen_position] = idx;

    // 重新插入被移动的键
    return GCTHelper(elements, positions, old_idx, attempt + 1, rng);
}

size_t hash_and_mod_with_mt19937(
    size_t id, 
    size_t nonce, 
    const std::vector<unsigned char>& data, 
    size_t modulus
) {
    // Use SHA-256 for hashing
    SHA256_CTX sha256;
    unsigned char hash[SHA256_DIGEST_LENGTH];

    // Initialize SHA-256 context
    SHA256_Init(&sha256);

    // Hash the custom string (id and nonce)
    std::stringstream ss;
    ss << id << nonce;
    SHA256_Update(&sha256, ss.str().c_str(), ss.str().size());

    // Hash the data
    SHA256_Update(&sha256, data.data(), data.size());

    // Finalize the hash
    SHA256_Final(hash, &sha256);

    // Create a seed from the hash
    std::seed_seq seed(std::begin(hash), std::end(hash));
    std::mt19937 generator(seed);

    // Create a uniform distribution
    std::uniform_int_distribution<size_t> distribution(0, modulus - 1);

    // Generate a random number
    size_t random_number = distribution(generator);

    return random_number;
}

void print_neibormatrix(
    std::vector<std::vector<size_t>> NeiborMatrix
) {
    cout << "  ";
    for (size_t jj = 0; jj < NeiborMatrix[0].size(); jj++) {
        cout << jj << ' ';
    }
    cout << endl;
    for (size_t ii = 0; ii < NeiborMatrix.size(); ii++) {
        cout << ii << ' ';
        for (size_t jj = 0; jj < NeiborMatrix[0].size(); jj++) {
            cout << NeiborMatrix[ii][jj] << ' ';
        }
        cout << endl;
    }
}

void print_neibormatrix_without_index(
    std::vector<std::vector<size_t>> NeiborMatrix
) {
    for (size_t ii = 0; ii < NeiborMatrix.size(); ii++) {
        for (size_t jj = 0; jj < NeiborMatrix[0].size(); jj++) {
            cout << NeiborMatrix[ii][jj] << ' ';
        }
        cout << endl;
    }
}

void print_flag_list(
    std::vector<bool> flagList
) {
    for (size_t i = 0; i < flagList.size(); i++) {
        cout << flagList[i] << endl;
    }
}

void print_stack(
    std::stack<VertexOriented> Stack
) {
    while (!Stack.empty()) {
        cout << Stack.top().edge_id << " " << Stack.top().vertex_id << endl;
        Stack.pop();
    }
}

void print_hypergraph(
    std::vector<Vertex> HyperGraph
) {
    for (size_t i = 0 ; i < HyperGraph.size(); i++) {
        cout << HyperGraph[i].vertex_id << "-----> ";
        for (auto &item: HyperGraph[i].edge_table) {
            cout << item.first << ": ";
            cout << "\"" << item.second.edge_id << ": (" << item.second.vertex[0] << ", " << item.second.vertex[1] << ")\"  ";
        }
        cout << endl;
    }
}

void HyperGraphSetup(
    std::vector<Vertex> &HyperGraph
) {
    for (size_t i = 0 ; i < HyperGraph.size(); i++) {
        HyperGraph[i].vertex_id = i;
        HyperGraph[i].edge_table = std::unordered_map<size_t, Edge>();
    }
}

// 返回一个布尔值向量，保存每个 key （边）是否被完全剥离
std::vector<bool> peel_the_hypergraph(
    std::stack<VertexOriented> &Stack,
    std::vector<Vertex> HyperGraph,
    std::vector<std::vector<size_t>> NeiborMatrix
) {
    // New a queue
    std::queue<Vertex> VertexQueue;
    // Init the queue and fill the singleton to it
    for (size_t i = 0; i < HyperGraph.size(); i++) {
        if (HyperGraph[i].edge_table.size() == 1) {
            VertexQueue.push(HyperGraph[i]);
        }
    }
    // New the return value: a vector of bool
    std::vector<bool> res(NeiborMatrix.size(), false);
    
    // While there is a node j \in [m] such that the 
    // set ki \notin P | j \in {h1(ki), h2(ki), h3(ki)} is a singleton: 
    // Let ki be the element of that singleton, and push ki onto P.
    while (VertexQueue.size() > 0) {
        // 找到队首元素的相邻节点
        Vertex temp = VertexQueue.front();
        if (HyperGraph[temp.vertex_id].edge_table.size() != 1) {
            VertexQueue.pop();
            continue;
        }
        // assert(temp.edge_table.size() == 1);
        auto Neibor = temp.edge_table.begin();
        Edge NeiborEdge = Neibor->second;
        for (size_t i = 0; i < 2; i++) {
            // 删去对应的边
            HyperGraph[NeiborEdge.vertex[i]].edge_table.erase(NeiborEdge.edge_id);
            // 检查对应的边是否满足度为 1
            if (HyperGraph[NeiborEdge.vertex[i]].edge_table.size() == 1) {
                // 入队
                VertexQueue.push(HyperGraph[NeiborEdge.vertex[i]]);
            }
        }
        // 队首出队，相应边对应的 key 入栈，向量对应位置设置为 true
        Stack.push(VertexOriented(NeiborEdge.edge_id, temp.vertex_id));
        res[NeiborEdge.edge_id] = true;
        VertexQueue.pop();
    }

    return res;
}

std::vector<std::vector<size_t>> choice_rows_in_neiborMatirx(
    std::vector<bool> rows_flag,
    std::vector<std::vector<size_t>> M
) {
    assert(rows_flag.size() == M.size()); 
    std::vector<std::vector<size_t>> res;
    for (size_t ii = 0; ii < M.size(); ii++) {
        if (!rows_flag[ii]) {
            res.push_back(M[ii]);
        }
    }
    
    return res;
}

std::vector<MatPoly> choice_position_in_real_pt(
    std::vector<bool> rows_flag,
    std::vector<MatPoly> real_pt
) {
    assert(rows_flag.size() == real_pt.size()); 
    std::vector<MatPoly> res;
    for (size_t ii = 0; ii < real_pt.size(); ii++) {
        if (!rows_flag[ii]) {
            res.push_back(real_pt[ii]);
        }
    }
    
    return res;
}

std::vector<MatPoly> get_whole_solution(
    std::vector<MatPoly> solution,
    std::stack<VertexOriented> Stack,
    std::vector<MatPoly> real_pt,
    std::unordered_map<size_t, std::vector<size_t>> positions,
    std::vector<size_t> realIndexes
) {
    for (size_t ii = 0; ii < solution.size(); ii++) {
        solution[ii] = to_ntt(solution[ii]);
    }

    for (size_t jj = 0; jj < real_pt.size(); jj++) {
        real_pt[jj] = to_ntt(real_pt[jj]);
    }

    while (!Stack.empty()) {
        size_t idx = Stack.top().edge_id;
        size_t pos = Stack.top().vertex_id;
        size_t key = realIndexes[idx];
        size_t flag = -1;

        std::vector<size_t> position_choice = positions[key];
        for (size_t jj = 0; jj < position_choice.size(); ++jj) {
            if (position_choice[jj] == pos) {
                flag = jj;
                break;
            }
        }

        MatPoly temp = add(solution[position_choice[(flag+1)%3]], solution[position_choice[(flag+2)%3]]);
        solution[pos] = add(real_pt[idx], to_ntt(invert(from_ntt(temp))));

        Stack.pop();
    }

    for (size_t ii = 0; ii < solution.size(); ii++) {
        solution[ii] = from_ntt(solution[ii]);
    }

    return solution;
}

void print_some_vars(
    std::vector<std::vector<size_t>> NeiborMatrix,
    std::vector<Vertex> HyperGraph,
    std::vector<bool> flagList,
    std::stack<VertexOriented> Stack,
    std::vector<std::vector<size_t>> core2
) {
    print_neibormatrix(NeiborMatrix);
    print_hypergraph(HyperGraph);
    cout << "=======================================" << endl;

    print_flag_list(flagList);
    cout << "=======================================" << endl;

    print_stack(Stack);
    cout << "=======================================" << endl;

    print_neibormatrix_without_index(NeiborMatrix);
    cout << "=======================================" << endl;

    print_neibormatrix_without_index(core2);
}

void print_solutions(
    std::vector<MatPoly> solutions
) {
    for (size_t jj = 0; jj < solutions.size(); ++jj) {
        cout << jj << ": " << getValue(to_ntt(solutions[jj])) << endl;
    }
}

