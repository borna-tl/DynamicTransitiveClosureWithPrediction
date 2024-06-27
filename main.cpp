#include <algorithm>
#include <bits/stdc++.h>
#include <filesystem>
#include <fstream>
#include <queue>
#include <sstream>
#include <string>
#include <bitset>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


using namespace std;
using u32 = uint_least32_t;
using engine = std::mt19937;

typedef pair <uint32_t, uint32_t> ui_pair;

#define INF 0x3f3f3f3f
#define PROGRESS_STAMP 10 // define the progress bar count
#define PBSTR "++++++++++++++++++++++++++++++++++++++++++++++++++"
#define PBWIDTH 50





// add std::
// maybe use static inline?

void print_progress(double percentage) {
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

template <typename T>
T variance(const std::vector<T> &vec) {
	const size_t sz = vec.size();
	if (sz <= 1) {
		return 0.0;
	}

	// Calculate the mean
	const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

	// Now calculate the variance
	auto variance_func = [&mean, &sz](T accumulator, const T &val) {
		return accumulator + ((val - mean) * (val - mean) / (sz - 1));
	};

	return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template <typename T>
T sum(const std::vector<T> &vec) {
	return std::accumulate(vec.begin(), vec.end(), 0.0);
}

struct Logger {
	void reset() {
		num_reachable_queries = 0;
		curr_insertion_cnt = 0;
		curr_query_cnt = 0;
	}
	int test_id;
	string algorithm;
	uint32_t query_operations_cnt = 0;
	uint32_t insertion_operations_cnt = 0;
	string start_time;
	string end_time;
	vector<double_t> run_durations; // in milliseconds
	vector<double_t> query_durations;
	vector<double_t> insertion_durations;
	int64_t curr_query_cnt = 0;
	int64_t curr_insertion_cnt = 0;
	size_t hashed_output;
	uint32_t num_reachable_queries = 0;
};

struct Setting {

	string INPUT_FILE = "sample.txt";
	string META_FILE = "meta-sample.txt";
	string OUTPUT_FILE = "output.txt";
	string LOG_FILE = "log.txt";

	string ALGORITHM = "sv_1";
	uint32_t QUERY_PERCENTAGE = 50;
	uint32_t TEST_RUN_COUNT = 10;
	uint32_t TIMEOUT_SEC = 1800;
	uint32_t OPERATION_SEED = 1223;
	uint32_t QUERY_SEED = 2334;
	int64_t QUERY_TIMESTAMP = 0;

	uint32_t nodes = 0;
	uint32_t input_lines = 0;
};

struct Operation {
	Operation(){};
	~Operation(){};
	Operation(const bool is_query_, const pair<uint32_t, uint32_t> arguments_,
			  const int64_t timestamp_)
		: is_query(is_query_), arguments(arguments_), timestamp(timestamp_) {}
	
	Operation(const Operation &op) {
		is_query = op.is_query;
		arguments = op.arguments;
		timestamp = op.timestamp;

	} 

	bool const operator == (const Operation& p) const {
   		return is_query == p.is_query
		&& arguments.first == p.arguments.first
		&& arguments.second == p.arguments.second
		&& timestamp == p.timestamp;
	}

    bool const operator < (const Operation &p) const {
        return arguments.first < p.arguments.first || (arguments.first == p.arguments.first 
		&& arguments.second < p.arguments.second);
    }

	void const print() const {
		cout << ((is_query == true) ? "Query" : "Ins") << "| (" << arguments.first
				<< ", " << arguments.second << ") @ " << timestamp << endl;
	}

	void reset(){
		is_query = false;
		arguments.first = 0;
		arguments.second = 0;
		timestamp = 0;
	}
	
	bool is_query;
	pair<uint32_t, uint32_t> arguments;
	int64_t timestamp;

};

class reachabilityTree { // this is a simple incremental reachability tree for vertex s
public:
	reachabilityTree(const uint32_t id_, const uint32_t max_nodes_,
					 const vector<vector<uint32_t>> &out_edge,
					 const vector<vector<uint32_t>> &in_edge)
		: id(id_), max_nodes(max_nodes_) {
		if (max_nodes == 0) {
			// look for some better form of sending errors
			cerr << "Not expecting zero nodes" << endl;
			exit(0);
		}
		r_plus.resize(max_nodes);
		r_minus.resize(max_nodes);
		initialize(out_edge, r_plus);
		initialize(in_edge, r_minus);
	}

	~reachabilityTree() {}

	void initialize(const vector<vector<uint32_t>> &edge, vector<bool> &r) {
		vector<uint32_t> q;
		size_t pointer = 0;
		r[id] = true;
		q.push_back(id);
		uint32_t u;
		while (pointer < q.size()) {
			u = q[pointer];
			pointer++;
			for (const uint32_t i : edge[u]) {
				if (!r[i]) {
					r[i] = true;
					q.push_back(i);
				}
			}
		}
	}

	void update(const uint32_t u, const uint32_t v,
				const vector<vector<uint32_t>> &out_edge,
				const vector<vector<uint32_t>> &in_edge) {
		update_reachability(u, v, out_edge, r_plus); // source reachability
		update_reachability(v, u, in_edge, r_minus); // sink reachability
	}

	void print_reachability_list() {
		cout << "Reachable Nodes from " << id << ": ";
		for (size_t i = 0; i < max_nodes; i++)
			if (r_plus[i])
				cout << i << " ";
		cout << endl;
		cout << "Reachable Nodes to " << id << ": ";
		for (size_t i = 0; i < max_nodes; i++)
			if (r_minus[i])
				cout << i << " ";
		cout << endl;
	}

	bool reaches(const uint32_t u) { return r_plus[u]; }

	bool is_reachable_from(const uint32_t u) { return r_minus[u]; }

	const uint32_t id;

private:
	const uint32_t max_nodes;
	vector<bool> r_plus;
	vector<bool> r_minus;

	void update_reachability(const uint32_t u, const uint32_t v,
							 const vector<vector<uint32_t>> &edge,
							 vector<bool> &r) {
		if (r[v])
			return;
		if (!r[u])
			return;

		queue<uint32_t> q;
		uint32_t curr_node = v;

		r[curr_node] = true;
		q.push(curr_node);
		while (!q.empty()) {
			curr_node = q.front();
			q.pop();
			for (const auto& i : edge[curr_node]) {
				if (!r[i]) {
					r[i] = true;
					q.push(i);
				}
			}
		}
	}
};

class Algorithms {
public:
	Algorithms(const Setting &setting_, Logger &logg_)
		: setting(setting_), logg(logg_) {
		out_edge.assign(setting.nodes, vector<uint32_t>());
		in_edge.assign(setting.nodes, vector<uint32_t>());
	}

	virtual ~Algorithms(){};

	virtual bool answer_query(const Operation& op) = 0;

	void run(const vector<Operation> operations) {
		logg.reset();
		clock_t tStart = clock();
		size_t c_out = 0;
		int64_t query_time = 0;
		int64_t insertion_time = 0;

		for (const auto &x : operations) {
			if (x.is_query) {
				if (x.timestamp <= setting.QUERY_TIMESTAMP) {
					cerr << "Warning: No queries should be at time "
						 << setting.QUERY_TIMESTAMP << endl;
				}
				auto started = std::chrono::high_resolution_clock::now();
				bool result = answer_query(x);
				// if (result == true){
				// 	cout << "reachable query: " << x.arguments.first << " " <<
				// 	x.arguments.second << " @ " << x.timestamp << endl;
				// }
				query_time += chrono::duration_cast<std::chrono::nanoseconds>(
								  chrono::high_resolution_clock::now() - started)
								  .count();
				logg.num_reachable_queries += (result == true);
				logg.curr_query_cnt++;
				results.push_back(result);
			}
			else {
				auto started = std::chrono::high_resolution_clock::now();
				add_edge(x);
				insertion_time += chrono::duration_cast<std::chrono::nanoseconds>(
									  chrono::high_resolution_clock::now() - started)
									  .count();
				logg.curr_insertion_cnt++;
			}
			c_out++;
			if (c_out > PROGRESS_STAMP * operations.size() / 100) {
				c_out = 0;
				if ((double)(clock() - tStart) / CLOCKS_PER_SEC > setting.TIMEOUT_SEC) {
					cerr << "Timeout!" << endl;
					exit(0);
				}
				print_progress((double)(logg.curr_query_cnt + logg.curr_insertion_cnt) /
							   operations.size());
			}
		}
		cout << endl;
		
		size_t hash_output = hash<vector<bool>>{}(results);
		logg.query_durations.push_back(query_time / 1e9);
		logg.insertion_durations.push_back(insertion_time / 1e9);
		logg.hashed_output = hash_output;
	}

protected:
	const Setting &setting;
	Logger &logg;
	vector<vector<uint32_t>> out_edge;
	vector<vector<uint32_t>> in_edge;
	vector<bool> results;

	virtual void add_edge(const Operation& op) {
		out_edge[op.arguments.first].push_back(op.arguments.second);
		in_edge[op.arguments.second].push_back(op.arguments.first);
	}
};

class Bfs : public Algorithms {

public:
	Bfs(const Setting &setting_, Logger &logg_) : Algorithms(setting_, logg_) {
		visited_bfs.resize(setting.nodes, false);
	}

	bool calculate_bfs(const uint32_t u, const uint32_t v) {
		vector<uint32_t> q;
		size_t pointer = 0;
		uint32_t curr_node = u;
		visited_bfs[curr_node] = true;
		q.push_back(curr_node);
		while (pointer < q.size()) {
			curr_node = q[pointer];
			pointer++;
			for (const uint32_t i : out_edge[curr_node]) {
				if (!visited_bfs[i])
				{
					visited_bfs[i] = true;
					q.push_back(i);
				}
			}
		}
		bool ans = visited_bfs[v];
		for (const uint32_t i : q) {
			visited_bfs[i] = false;
		}
		return ans;
	}

	bool answer_query(const Operation& op) {
		return calculate_bfs(op.arguments.first, op.arguments.second);
	}

private:
	vector<bool> visited_bfs;
};

class Dfs : public Algorithms {

public:
	Dfs(const Setting &setting_, Logger &logg_) : Algorithms(setting_, logg_) {
		visited_dfs.resize(setting.nodes, false);
	}

	bool calculate_dfs(const uint32_t u, const uint32_t v) {
		visited_dfs[u] = true;
		if (visited_dfs[v]) {
			return true;
		}
		for (const auto& i : out_edge[u]) {
			if (!(visited_dfs[i])) {
				calculate_dfs(i, v);
			}
		}
		return visited_dfs[v];
	}

	// bool calculate_dfs_with_path(const uint32_t u, const uint32_t v) {
	// 	visited_dfs[u] = true;
	// 	current_path.push_back(u);  // Add the current node to the path

	// 	if (u == v) {
	// 		found_path = current_path;  // Store the path if the destination is found
	// 		return true;
	// 	}

	// 	for (const auto& i : out_edge[u]) {
	// 		if (!visited_dfs[i]) {
	// 			if (calculate_dfs(i, v)) {
	// 				return true;  // If the destination is found in the recursion, return true
	// 			}
	// 		}
	// 	}

	// 	current_path.pop_back();  // Remove the current node from the path if not part of the final path
	// 	return false;
	// }

	bool answer_query(const Operation& op) {
		visited_dfs.assign(visited_dfs.size(), false);
		return calculate_dfs(op.arguments.first, op.arguments.second);
	}

private:
	vector<bool> visited_dfs;
};

class mini_dfs {

public:
	mini_dfs(const uint32_t nodes) { //maybe mini bfs is faster
		visited_dfs.resize(nodes, false);
		out_edge.assign(nodes, vector<uint32_t>());

	}

	bool calculate_dfs(const uint32_t u, const uint32_t v) {
		visited_dfs[u] = true;
		if (visited_dfs[v]) {
			return true;
		}
		for (const auto& i : out_edge[u]) {
			if (!(visited_dfs[i])) {
				calculate_dfs(i, v);
			}
		}
		return visited_dfs[v];
	}
	bool answer_query(const uint32_t u, const uint32_t v) {
		visited_dfs.assign(visited_dfs.size(), false);
		return calculate_dfs(u, v);
	}
	void add_edge(const uint32_t u, const uint32_t v) {
		out_edge[u].push_back(v);
	}
	void print_edges(){
		for (uint32_t i = 0; i < out_edge.size(); i++){
			cout << "for node: " << i << endl;
			for (uint32_t j = 0; j < out_edge[i].size(); j++)
				cout << out_edge[i][j] << " ";
			cout << endl;
		}
	}
private:
	vector<bool> visited_dfs;
	vector<vector<uint32_t>> out_edge;

};


class Bibfs : public Algorithms {

public:
	Bibfs(const Setting &setting_, Logger &logg_) : Algorithms(setting_, logg_) {
		visited_bibfs_source.resize(setting.nodes, false);
		visited_bibfs_sink.resize(setting.nodes, false);
	}
	bool answer_query(const Operation& op) {
		return calculate_bibfs(op.arguments.first, op.arguments.second);
	}

private:
	vector<bool> visited_bibfs_source;
	vector<bool> visited_bibfs_sink;

	bool calculate_bibfs(const uint32_t u, const uint32_t v) {
		bool found_path = (u == v); //for special case where u and v are the same
		uint32_t curr_node;
		vector<uint32_t> source_queue, sink_queue;
		size_t source_pointer = 0;
		size_t sink_pointer = 0;
		visited_bibfs_source[u] = true;
		visited_bibfs_sink[v] = true;
		source_queue.push_back(u);
		sink_queue.push_back(v);
		while (!found_path && source_pointer < source_queue.size() &&
			   sink_pointer < sink_queue.size()) {
			// running bfs for the source queue one time
			curr_node = source_queue[source_pointer];
			source_pointer++;
			for (const auto& i : out_edge[curr_node]) {
				if (!visited_bibfs_source[i]) {
					visited_bibfs_source[i] = true;
					source_queue.push_back(i);
				}
				if (visited_bibfs_source[i] && visited_bibfs_sink[i]) {
					found_path = true;
				}
			}
			// running bfs for the back queue one time
			curr_node = sink_queue[sink_pointer];
			sink_pointer++;
			for (const auto& i : in_edge[curr_node]) {
				if (!visited_bibfs_sink[i]) {
					visited_bibfs_sink[i] = true;
					sink_queue.push_back(i);
				}
				if (visited_bibfs_source[i] && visited_bibfs_sink[i]) {
					found_path = true;
				}
			}
		}
		for (const uint32_t i : source_queue) {
			visited_bibfs_source[i] = false;
		}
		for (const uint32_t i : sink_queue) {
			visited_bibfs_sink[i] = false;
		}
		return found_path;
	}
};

class Sv : public Algorithms {

public:
	Sv(size_t count_, const Setting &setting_, Logger &logg_, uint32_t sv_seed_)
		: Algorithms(setting_, logg_), count(count_) {
		visited_bibfs_source.resize(setting.nodes, false);
		visited_bibfs_sink.resize(setting.nodes, false);
		sv_seed = sv_seed_;
	}
	virtual ~Sv() {}

	bool calculate_sv(const uint32_t u, const uint32_t v) {
		// instead of searching, we can use hash map as well.
		// since |sv| is little, i guess it's better to simply search
		auto it = find_sv(u);
		if (it != reachability_tree.end()) {
			return (*it)->reaches(v);
		}
		it = find_sv(v);
		if (it != reachability_tree.end()) {
			return (*it)->is_reachable_from(u);
		}
		for (const auto &rt : reachability_tree) {
			// obs. 1
			if (rt->is_reachable_from(u) && rt->reaches(v)) {
				return true;
			}
			// obs. 2
			if (rt->reaches(u) && !rt->reaches(v)) {
				return false;
			}
			// obs. 3
			if (rt->is_reachable_from(v) && !rt->is_reachable_from(u)) {
				return false;
			}
		}

		return calculate_bibfs(u, v);
	}
	bool answer_query(const Operation& op) {
		if (!generated_sv) {
			generate_sv_list();
			generated_sv = true;
		}
		return calculate_sv(op.arguments.first, op.arguments.second);

	}

private:
	vector<unique_ptr<reachabilityTree>> reachability_tree;
	// bringing bibfs fallback algorithm inside sv (because we need the same graph out/in edges)
	vector<bool> visited_bibfs_source;
	vector<bool> visited_bibfs_sink;
	size_t count;
	uint32_t sv_seed;
	bool generated_sv = false;

	const vector<unique_ptr<reachabilityTree>>::iterator find_sv(uint32_t sv) {
		return find_if(reachability_tree.begin(), reachability_tree.end(),
					   [&sv](const unique_ptr<reachabilityTree> &obj) {
						   return (*obj).id == sv;
					   });
	}

	void generate_candidate_svs(vector<uint32_t> &non_isolated_svs,
								vector<uint32_t> &half_isolated_svs) {
		for (size_t i = 0; i < setting.nodes; i++) {
			if (is_non_isolate(i))
				non_isolated_svs.push_back(i);
			else if (is_half_isolate(i))
				half_isolated_svs.push_back(i);
		}
	}

	void generate_sv_list() {
		// random_device os_seed;
		engine generator(sv_seed);
		vector<uint32_t> non_isolated_svs, half_isolated_svs;
		generate_candidate_svs(non_isolated_svs, half_isolated_svs);

		size_t i = 0;
		logg.algorithm += '{';
		while (i < count) {
			uint32_t sv = 0;

			if (non_isolated_svs.size() > 0) {
				uniform_int_distribution<u32> non_distribute(
					0, non_isolated_svs.size() - 1);
				size_t index = non_distribute(generator);
				sv = non_isolated_svs[index];
				non_isolated_svs.erase(non_isolated_svs.begin() + index);
			}
			else if (half_isolated_svs.size() > 0) {
				uniform_int_distribution<u32> half_distribute(
					0, half_isolated_svs.size() - 1);
				size_t index = half_distribute(generator);
				sv = half_isolated_svs[index];
				half_isolated_svs.erase(half_isolated_svs.begin() + index);
			}
			else {
				uniform_int_distribution<u32> distribute(0, setting.nodes - 1);
				sv = distribute(generator);
			}

			if (find_sv(sv) != reachability_tree.end()) {
				continue;
			}
			// reachability_tree[sv] = new reachabilityTree(sv);
			reachability_tree.push_back(unique_ptr<reachabilityTree>(
				new reachabilityTree(sv, setting.nodes, out_edge, in_edge)));
			// reachability_tree.push_back(make_unique<reachabilityTree>(sv));
			logg.algorithm += to_string(sv) + ", ";
			i++;
		}
		logg.algorithm += '}';
	}

	bool is_non_isolate(const uint32_t u) {
		return (out_edge[u].size() != 0 && in_edge[u].size() != 0);
	}

	bool is_half_isolate(const uint32_t u) {
		return (out_edge[u].size() != 0 || in_edge[u].size() != 0);
	}

	void update_sv(const uint32_t u, const uint32_t v) {
		for (const auto &rt : reachability_tree) {
			rt->update(u, v, out_edge, in_edge);
		}
	}

	void add_edge(const Operation& op) {
		uint32_t u, v;
		u = op.arguments.first;
		v = op.arguments.second;
		out_edge[u].push_back(v);
		in_edge[v].push_back(u);
		if (generated_sv)
			update_sv(u, v);
	}

	bool calculate_bibfs(const uint32_t u, const uint32_t v) {
		bool found_path = (u == v);
		uint32_t curr_node;
		vector<uint32_t> source_queue, sink_queue;
		size_t source_pointer = 0;
		size_t sink_pointer = 0;
		visited_bibfs_source[u] = true;
		visited_bibfs_sink[v] = true;
		source_queue.push_back(u);
		sink_queue.push_back(v);
		while (!found_path && source_pointer < source_queue.size() &&
			   sink_pointer < sink_queue.size()) {
			// running bfs for the source queue one time
			curr_node = source_queue[source_pointer];
			source_pointer++;
			for (const auto& i : out_edge[curr_node]) {
				if (!visited_bibfs_source[i]) {
					visited_bibfs_source[i] = true;
					source_queue.push_back(i);
				}
				if (visited_bibfs_source[i] && visited_bibfs_sink[i]) {
					found_path = true;
				}
			}
			// running bfs for the back queue one time
			curr_node = sink_queue[sink_pointer];
			sink_pointer++;
			for (const auto& i : in_edge[curr_node]) {
				if (!visited_bibfs_sink[i]) {
					visited_bibfs_sink[i] = true;
					sink_queue.push_back(i);
				}
				if (visited_bibfs_source[i] && visited_bibfs_sink[i]) {
					found_path = true;
				}
			}
		}
		for (const uint32_t i : source_queue) {
			visited_bibfs_source[i] = false;
		}
		for (const uint32_t i : sink_queue) {
			visited_bibfs_sink[i] = false;
		}
		return found_path;
	}
};


class Pred : public Algorithms {

public:
	Pred(const Setting &setting_, Logger &logg_, vector<Operation>& pred_insertions_)
			 : Algorithms(setting_, logg_), pred_insertions(pred_insertions_){
		

		last_seen_index = -1;

		indices_in_pred.clear();
		indices_in_pred.reserve(logg.insertion_operations_cnt);
		inserted.resize(logg.insertion_operations_cnt);

		for (size_t i = 0; i < logg.insertion_operations_cnt; i++){
			auto result = indices_in_pred.try_emplace(pred_insertions[i].arguments.first*(int64_t)1e9+pred_insertions[i].arguments.second,
				i);			
			
			if (result.second == false){
				inserted[i] = true;
			}
			else
				inserted[i] = false;
		}

	}
	void preheat() {
		calculate_d();
	}

	bool calculate_pred(const uint32_t u, const uint32_t v) {
		vector<ui_pair> nodes = {make_pair(u, u), make_pair(v, v)};

		for (auto x : edges_for_dfs){
			nodes.push_back(pred_insertions[x].arguments);
		}

		size_t nodes_s = nodes.size();
 		queue<int> q;
		vector<bool> visited(nodes_s, false);
		q.push(0);
		visited[0] = true;
		while (!q.empty()) {
			int u = q.front();
			q.pop();
			if (u == 1)
				return true;
			uint32_t start = nodes[u].second;
			for (size_t j = 0; j < nodes_s; j++) {
				if (!visited[j] && bottle_neck[start][nodes[j].first] <= last_seen_index) {
					q.push(j);
					visited[j] = true;
				}
				
			}
		}
		return false;
	}

	bool answer_query(const Operation& op) {
		uint32_t u, v;
		u = op.arguments.first;
		v = op.arguments.second;

		update_lcs();

		//makes it really fast!
		if (bottle_neck[u][v] <= last_seen_index){	
			return true;
		}

		return calculate_pred(u, v);
	}

private:
	
	vector<Operation>& pred_insertions; //(u, v) edit needed //change to pointer
	unordered_map<int64_t, int32_t> indices_in_pred;
	vector <bool> inserted;
	unordered_set <int> edges_for_dfs;
	vector<ui_pair>* edge_insertion_time;
	vector<vector<int32_t>> bottle_neck;
	int last_seen_index = -1; //is the index of the last seen element (starting from 0)

	bool has_op(const Operation& x, 
				const vector<Operation>::iterator start,
				const vector<Operation>::iterator end){
		
		for (auto i = start; i < end; i++){
			if (i->arguments.first == x.arguments.first 
				&& i->arguments.second == x.arguments.second 
				&& i->is_query == x.is_query)
				return true;
		}
		return false;
	}
	void update_lcs() { //returns the first position predicted wrong
		while (inserted[last_seen_index + 1] == true){ //(u, v)
			last_seen_index++;
			edges_for_dfs.erase(last_seen_index);
		}
	}

	void set_t(){ //t[u][v] when edge (u, v) will be inserted
		uint32_t u, v;
		for (size_t i = 0; i < pred_insertions.size(); ++i) {
			u = pred_insertions[i].arguments.first;
			v = pred_insertions[i].arguments.second;
			edge_insertion_time[u].push_back(make_pair(v, (uint32_t)i));
		}
	}

	void calculate_d(){ //calculate D[u][v] = bottle_neck (last edge addition) of earlieast u->v path
		edge_insertion_time = new vector<ui_pair> [setting.nodes];
		bottle_neck.assign(setting.nodes, vector<int32_t>());
		set_t();
		for (size_t i = 0; i < setting.nodes; i++){
			store_shortest_path(i);
			if (rand() % 100 == 1){
				cout << "calculated shortest paths from " << i << endl;
			}
		}
		
		delete[] edge_insertion_time; //return and dynamic pointer later
	}

	void store_shortest_path (uint32_t src) {
		priority_queue< ui_pair, vector <ui_pair> , greater<ui_pair> > pq;

		// Create a vector for distances and initialize all
		// distances as infinite (INF)
		// vector<uint32_t> dist(setting.nodes, INF);
		bottle_neck[src].assign(setting.nodes, INF);
		// Insert source itself in priority queue and initialize its distance as 0.
		pq.push(make_pair(0, src));
		bottle_neck[src][src] = 0;
		
		vector<bool> f(setting.nodes, false);

	
		while (!pq.empty()) {
			
			uint32_t u = pq.top().second;
			pq.pop();
			if (f[u])
				continue;
			f[u] = true;

			vector< ui_pair >::iterator i;
			for (i = edge_insertion_time[u].begin(); i != edge_insertion_time[u].end(); ++i) {
				
				uint32_t v = (*i).first;
				int32_t timestamp = (*i).second;

				// If there is shorted path to v through u.
				if (bottle_neck[src][v] > max(bottle_neck[src][u], timestamp)) {
					// Updating distance of v
					bottle_neck[src][v] = max(bottle_neck[src][u], timestamp);
					pq.push(make_pair(bottle_neck[src][v], v));
				
				}
			}
		}

	}


	void add_edge(const Operation& op) {
		uint32_t index = indices_in_pred[op.arguments.first*(int64_t)1e9+op.arguments.second];

		if (inserted[index] != true) {
			inserted[index] = true;
			edges_for_dfs.insert(index);
		}
	}
};

void set_time(string &t) {
	auto timepoint = chrono::system_clock::now();
	auto coarse = chrono::system_clock::to_time_t(timepoint);
	auto fine = chrono::time_point_cast<std::chrono::milliseconds>(timepoint);

	char buffer[sizeof "9999-12-31 23:59:59.999"];
	std::snprintf(buffer + std::strftime(buffer, sizeof buffer - 3, "%F %T.",
										 std::localtime(&coarse)),
				  4, "%03lu", fine.time_since_epoch().count() % 1000);
	t = buffer;
}

class Program {
public:
	Program() {}
	~Program(){}
	void read_parse_input(const int argc, const char *argv[]) {
		for (int i = 1; i < argc; i += 2) {
			if (!strcmp(argv[i], "-alg")) {
				if (strcmp(argv[i + 1], "dfs") && strcmp(argv[i + 1], "bfs") &&
					strcmp(argv[i + 1], "bibfs") &&
					strcmp(argv[i + 1], "pred") &&
					string(argv[i + 1]).substr(0, 3) != "sv_") {
					cerr << "Wrong Input for Algorithm.\n";
					exit(0);
				}
				setting.ALGORITHM = argv[i + 1];
			}
			if (!strcmp(argv[i], "-qp")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1]) ||
					stoi(argv[i + 1]) < 0 || stoi(argv[i + 1]) > 100) {
					cerr << "Wrong Input for Query Percentage.\n";
					exit(0);
				}
				setting.QUERY_PERCENTAGE = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-trc")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1])) {
					cerr << "Wrong Input for Test Run Count.\n";
					exit(0);
				}
				setting.TEST_RUN_COUNT = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-ts")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1])) {
					cerr << "Wrong Input for Timeout Seconds.\n";
					exit(0);
				}
				setting.TIMEOUT_SEC = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-os")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1])) {
					cerr << "Wrong Input for Operation Seed.\n";
					exit(0);
				}
				setting.OPERATION_SEED = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-qs")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1])) {
					cerr << "Wrong Input for Query Seed.\n";
					exit(0);
				}
				setting.QUERY_SEED = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-qt")) {
				if (strspn(argv[i + 1], "-.0123456789") != strlen(argv[i + 1])) {
					cerr << "Wrong Input for Query Timestamp.\n";
					exit(0);
				}
				setting.QUERY_TIMESTAMP = stoi(argv[i + 1]);
			}
			if (!strcmp(argv[i], "-inp")) {
				std::ifstream infile(argv[i + 1]);
				if (!infile.good()) {
					cerr << "Input File Does Not Exist.\n";
					exit(0);
				}
				setting.INPUT_FILE = argv[i + 1];
			}
			if (!strcmp(argv[i], "-meta")) {
				std::ifstream infile(argv[i + 1]);
				if (!infile.good()) {
					cerr << "Meta File Does Not Exist.\n";
					exit(0);
				}
				setting.META_FILE = argv[i + 1];
			}
			if (!strcmp(argv[i], "-out")) {
				std::ifstream infile(argv[i + 1]);
				if (!infile.good()) {
					cerr << "Warning: Output File Did Not Exist. One Will Be Created.\n";
				}
				setting.OUTPUT_FILE = argv[i + 1];
			}
			if (!strcmp(argv[i], "-log")) {
				std::ifstream infile(argv[i + 1]);
				if (!infile.good()) {
					cerr << "Warning: Log File Did Not Exist. One Will Be Created.\n";
				}
				setting.LOG_FILE = argv[i + 1];
			}
		}

		read_input_file(); // read dataset and build input
		operations.reserve(setting.input_lines *
						   ((100 + setting.QUERY_PERCENTAGE) / 100));
		generate_operations(insertion_operations, operations); // add query operatios to input and store result in
							   // "operations"
		cout << "operations generated" << endl;
		logg.test_id = get_test_id();
		logg.algorithm = setting.ALGORITHM;
	}

	void execute_test(){
		if (setting.ALGORITHM == "pred")
			execute_pred();
		else
			execute_reg();
	}

	void execute_pred() {
		set_time(logg.start_time);
		for (size_t i = 0; i < setting.TEST_RUN_COUNT; i++) {
			unique_ptr<Pred> pred;
			permute_insertions(i);

			// cout << "OPERATIONS\n\n" << endl;
			// for (auto x : insertion_operations)
			// 	x.print();
			// cout << "PERMUTED OPERATIONS\n\n" << endl;
			// for (auto x : insertion_operations_permuted)
			// 	x.print();
			// cout << "\n\n";
			pred = unique_ptr<Pred>(new Pred(setting, logg, insertion_operations_permuted));
			pred->preheat();			

			pred->run(operations);
			logg.run_durations.push_back(logg.query_durations.back() +
										 logg.insertion_durations.back());
		}
		set_time(logg.end_time);
		// cout << "Search time: " << search_time << endl <<
		// "Clear and Reserve time: " << clear_reserve_time << endl <<
		// "Emplace time:" << emplace_time << endl <<
		// "LCS time:" << update_lcs_time << endl <<
		// "Adding Minidfs nodes:" << adding_minidfs_nodes << endl <<
		// "Declaring Minidfs:" << declaring_minidfs << endl <<
		// "Mini DFS 1st edge:" << minidfs_first_edge_set << endl <<
		// "Mini DFS 2nd edge:" << minidfs_second_edge_set << endl <<
		// "Building Mini DFS time:" << building_minidfs_time << endl;

	}

	void execute_reg() {
		set_time(logg.start_time);
		for (size_t i = 0; i < setting.TEST_RUN_COUNT; i++) {
			unique_ptr<Algorithms> alg;
			if (setting.ALGORITHM == "dfs")
				alg = unique_ptr<Algorithms>(new Dfs(setting, logg));
			else if (setting.ALGORITHM == "bfs")
				alg = unique_ptr<Algorithms>(new Bfs(setting, logg));
			else if (setting.ALGORITHM == "bibfs")
				alg = unique_ptr<Algorithms>(new Bibfs(setting, logg));
			else if (setting.ALGORITHM.substr(0, 3) == "sv_") {
				alg = unique_ptr<Algorithms>(
					new Sv(stoi(setting.ALGORITHM.substr(3)), setting, logg, i));
			}
			else
				assert(false);
			alg->run(operations);
			logg.run_durations.push_back(logg.query_durations.back() +
										 logg.insertion_durations.back());
		}
		set_time(logg.end_time);
	}

	void write_to_log() {
		double_t duration = sum(logg.run_durations);
		double_t queries = sum(logg.query_durations);
		double_t insertions = sum(logg.insertion_durations);

		ofstream log_file;
		log_file.open(setting.LOG_FILE, std::ios_base::app);
		log_file << "test id: " << logg.test_id << '\n'
				 << "input file name: " << setting.INPUT_FILE << '\n'
				 << "input file lines: " << setting.input_lines << '\n'
				 << "seed: " << setting.OPERATION_SEED << ", " << setting.QUERY_SEED
				 << '\n'
				 << "#run: " << setting.TEST_RUN_COUNT << '\n'
				 << "algorithm: " << logg.algorithm << '\n'
				 << "#query opertions: " << logg.query_operations_cnt << '\n'
				 << "#insertion operations: " << logg.insertion_operations_cnt
				 << '\n'
				 << "start time: " << logg.start_time << '\n'
				 << "end time: " << logg.end_time << '\n'
				 << "duration: " << duration << "{";
		for (const auto& x : logg.run_durations) {
			log_file << x << " ";
		}
		log_file << '}' << ' ' << '[' << variance(logg.run_durations) << ']' << '\n'
				 << "queries: " << queries << "{";
		for (const auto& x : logg.query_durations) {
			log_file << x << " ";
		}
		log_file << '}' << ' ' << '[' << variance(logg.query_durations) << ']'
				 << '\n'
				 << "insertions: " << insertions << "{";
		for (const auto& x : logg.insertion_durations) {
			log_file << x << " ";
		}
		log_file << '}' << ' ' << '[' << variance(logg.insertion_durations) << ']'
				 << '\n'
				 << "hashed output: " << logg.hashed_output << '\n'
				 << "#reachable queries: " << logg.num_reachable_queries << '\n'
				 << string(50, '*') << "\n";
	}
	
	
private:
	Setting setting;
	Logger logg;
	vector<Operation> insertion_operations;
	vector<Operation> insertion_operations_permuted;
	vector<Operation> operations; //insertions or queries

	void convert_input() {
		
		ifstream file(setting.INPUT_FILE);
		std::filesystem::path p(setting.INPUT_FILE);

		if (filesystem::exists("build/" + string(p.filename()))) {
			cout << "Directory already exists." << endl;
		}
		else {
			filesystem::create_directories("build");
		}

		ofstream new_file("build/" + string(p.filename()), std::ios::trunc);

		uint32_t max_nodes = 0;
		uint32_t lines = 0;

		uint32_t u, v;
		string type, timestamp;
		while (file >> u >> v >> type >> timestamp) {

			max_nodes = max(max_nodes, max(u, v));
			if (type == "+1") {
				new_file << u << ' ' << v << ' ' << timestamp << endl;
				lines++;
			}
		}
		setting.nodes = max_nodes + 1;
		setting.input_lines = lines;

		new_file.close();
	}

	void read_input_file() {
		convert_input();
		cout << "Input generated!" << endl;

		std::filesystem::path p(setting.INPUT_FILE);
		ifstream infile("build/" + string(p.filename()));
		uint32_t u, v;
		int64_t timestamp;

		while (infile >> u >> v >> timestamp) {
			insertion_operations.push_back(
				Operation(false, make_pair(u, v), timestamp));
		}
		infile.close();
		cout << "input done" << endl;
	}

	size_t mod (size_t x, size_t y){
		size_t m = x % y;
		if (m < 0)
			m +=  y;
		return m; 
	}

	void permute_insertions(uint32_t sv_seed) {
		engine generator(sv_seed);
		normal_distribution<double> distribute(0, 3.0); //might have to calibrate

		insertion_operations_permuted.resize((int)logg.insertion_operations_cnt);
		vector <pair<double, int>> new_order; 
		for (int i = 0; i < (int)logg.insertion_operations_cnt ; i++){
			insertion_operations_permuted[i].reset();
		}
		insertion_operations_permuted = insertion_operations;

		int i;
		for (i = 0; i < (int)logg.insertion_operations_cnt; i ++){
			double noise = distribute(generator);
			double new_index = i + noise;
			new_order.push_back({new_index, i});
		}
		sort(new_order.begin(), new_order.end());
		i = 0;
		for (auto x : new_order){ 
			insertion_operations_permuted[i].arguments.first = insertion_operations[x.second].arguments.first;
			insertion_operations_permuted[i].arguments.second = insertion_operations[x.second].arguments.second;
			i++;
		}
	}
	

	void generate_operations(vector<Operation>& insertions, vector<Operation>& operations) {
		engine generator(setting.OPERATION_SEED);
		engine query_generator(setting.QUERY_SEED);
		uniform_int_distribution<u32> distribute(0, 99);
		uniform_int_distribution<u32> query_chance_distribute(0, setting.nodes - 1);
		uint32_t u, v;
		logg.query_operations_cnt = 0;
		logg.insertion_operations_cnt = 0;
		for (const auto &x : insertions) {
			u = x.arguments.first;
			v = x.arguments.second;
			// currently each insertion and query have the same timestamp
			// we can modify for different implementations later


			operations.push_back(Operation(false, make_pair(u, v), x.timestamp));
			//revert back to while for lower %qp

			while (x.timestamp > setting.QUERY_TIMESTAMP &&
				distribute(generator) < setting.QUERY_PERCENTAGE) {
				uint32_t u_q = query_chance_distribute(query_generator);
				uint32_t v_q = query_chance_distribute(query_generator);
				operations.push_back(Operation(true, make_pair(u_q, v_q), x.timestamp));
				logg.query_operations_cnt++;
			}
		}
		logg.insertion_operations_cnt = insertions.size();
	}

	int get_test_id() { // can probably read backwards later to enhance speed
		string line;
		ifstream log_file(setting.LOG_FILE);
		int id = 0;
		while (getline(log_file, line)) {
			if (line.find('*') != string::npos) {
				id++;
			}
		}
		return id;
	}
};

int main(int argc, const char *argv[]) {

	Program program;
	program.read_parse_input(argc, argv);

	program.execute_test();

	program.write_to_log();
	
}