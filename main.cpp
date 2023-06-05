#include <bits/stdc++.h>
#include <fstream>
#include <queue>
#include <string>
#include <sstream> 
#include <filesystem>
#include <algorithm>


using namespace std;
using u32    = uint_least32_t;
using engine = std::mt19937;


#define PROGRESS_STAMP 10 //define the progress bar count
#define PBSTR "++++++++++++++++++++++++++++++++++++++++++++++++++"
#define PBWIDTH 50


//add const where applicable
//add std::
//maybe use static inline?


void print_progress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

template<typename T>
T variance(const std::vector<T> &vec) {
    const size_t sz = vec.size();
    if (sz <= 1) {
        return 0.0;
    }

    // Calculate the mean
    const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    // Now calculate the variance
    auto variance_func = [&mean, &sz](T accumulator, const T& val) {
        return accumulator + ((val - mean)*(val - mean) / (sz - 1));
    };

    return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
T sum (const std::vector<T> &vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0);
}


struct Logger {
    void reset(){
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
    vector<double_t> run_durations; //in milliseconds
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
    uint32_t QUERY_PERCENTAGE = 100;
    uint32_t TEST_RUN_COUNT = 10;
    uint32_t TIMEOUT_SEC = 1800;
    uint32_t OPERATION_SEED = 1223;
    uint32_t QUERY_SEED = 2334;
    int64_t QUERY_TIMESTAMP = 0;


    uint32_t nodes = 0;
    uint32_t input_lines = 0;
};


struct Operation {
    Operation(const bool is_query_, const pair <uint32_t, uint32_t> arguments_,
                const int64_t timestamp_): is_query(is_query_), arguments(arguments_),
                timestamp(timestamp_){}
    bool is_query;
    pair <uint32_t, uint32_t> arguments;
    int64_t timestamp;
};


class reachabilityTree{ //this is a simple incremental reachability tree for vertex s
public:
    reachabilityTree(const uint32_t id_, const uint32_t max_nodes_, 
                        const vector<vector<uint32_t>> &out_edge,
                        const vector<vector<uint32_t>> &in_edge):
                        id(id_), max_nodes(max_nodes_){
        if (max_nodes == 0){
            //look for some better form of sending errors
            cerr << "Not expecting zero nodes" << endl;
            exit(0);
        }
        r_plus.resize(max_nodes);
        r_minus.resize(max_nodes);
        initialize(out_edge, r_plus);
        initialize(in_edge, r_minus);

    }
    
    ~reachabilityTree(){}

    void initialize(const vector<vector<uint32_t>> &edge, vector<bool>& r){
        vector<uint32_t> q;
        size_t pointer = 0;
        r[id] = true;
        q.push_back(id); 
        uint32_t u;
        while (pointer < q.size()) {
            u = q[pointer];
            pointer++;
            for (const uint32_t i : edge[u]){
                if (!r[i]){
                    r[i] = true;
                    q.push_back(i);
                }
            }
        }
        
    }
    
    void update(const uint32_t u, const uint32_t v,
                const vector<vector<uint32_t>> &out_edge,
                const vector<vector<uint32_t>> &in_edge){
        update_reachability(u, v, out_edge, r_plus); //source reachability
        update_reachability(v, u, in_edge, r_minus); //sink reachability
    }
    
    void print_reachability_list(){
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
    
    bool reaches(const uint32_t u){
        return r_plus[u];
    }
    
    bool is_reachable_from(const uint32_t u){
        return r_minus[u];
    }
    
    const uint32_t id; 
private:
    const uint32_t max_nodes; 
    vector <bool> r_plus;
    vector <bool> r_minus;

    void update_reachability(const uint32_t u, const uint32_t v,
                                const vector<vector<uint32_t>> &edge, vector<bool>& r){
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
            for (const auto i : edge[curr_node]){
                if (!r[i]){
                    r[i] = true;
                    q.push(i);
                }
            }
        }
    }
};

class Algorithms{
public:
    Algorithms(const Setting& setting_, Logger& logg_) : setting(setting_), logg(logg_){
        out_edge.assign(setting.nodes, vector<uint32_t>());
        in_edge.assign(setting.nodes, vector<uint32_t>());
    }
    
    virtual ~Algorithms(){};
    
    virtual bool answer_query(const uint32_t u, const uint32_t v) = 0; 
    
    void run(const vector<Operation> operations){
        logg.reset();
        uint32_t u, v;
        clock_t tStart = clock();
        size_t c_out = 0;
        int64_t query_time = 0;
        int64_t insertion_time = 0;

        for (const auto &x : operations){
            u = x.arguments.first;
            v = x.arguments.second;
            if (x.is_query){
                if (x.timestamp <= setting.QUERY_TIMESTAMP){
                    cerr << "Warning: No queries should be at time " << setting.QUERY_TIMESTAMP << endl;
                    cerr << u << " " << v << endl;
                }
                auto started = std::chrono::high_resolution_clock::now();
                bool result = answer_query(u, v);  

                query_time += chrono::duration_cast<std::chrono::nanoseconds>
                            (chrono::high_resolution_clock::now()-started).count();
                logg.num_reachable_queries += (result == true); 
                logg.curr_query_cnt ++;
                results.push_back(result);
            }
            else{
                auto started = std::chrono::high_resolution_clock::now();
                add_edge(u, v);
                insertion_time += chrono::duration_cast<std::chrono::nanoseconds>
                            (chrono::high_resolution_clock::now()-started).count();  
                logg.curr_insertion_cnt ++;              
            }
            c_out++;
            if (c_out > PROGRESS_STAMP * operations.size() / 100){
                c_out = 0;
                if ((double)(clock() - tStart)/CLOCKS_PER_SEC > setting.TIMEOUT_SEC){
                    cerr << "Timeout!" << endl;
                    exit(0);
                }
                print_progress((double)(logg.curr_query_cnt+logg.curr_insertion_cnt)/operations.size());
            }
        }
        cout << endl;
        // printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        // cout << "Queries answered: " << queries_answered << endl;
        // cout << "Reachable queries: " << true_q << endl;
        // ofstream outfile(OUTPUT_FILE, ios_base::app);
        // ostream_iterator<string> output_iterator(outfile, "");
        size_t hash_output = hash<vector<bool>>{}(results);
        // outfile << hash<vector<bool>>{}(results) << "\n";
        // outfile.close();
        logg.query_durations.push_back(query_time / 1e9);
        logg.insertion_durations.push_back(insertion_time / 1e9);
        logg.hashed_output = hash_output;
    }

protected:
    const Setting& setting;
    Logger& logg;
    vector<vector<uint32_t>> out_edge;
    vector<vector<uint32_t>> in_edge;
    vector <bool> results;
    
    virtual void add_edge(const uint32_t u, const uint32_t v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
    }
    
};


class Bfs : public Algorithms{

public:
    Bfs(const Setting& setting_, Logger& logg_) : Algorithms(setting_, logg_){
        visited_bfs.resize(setting.nodes, false);
    }
    
    bool calculate_bfs(const uint32_t u, const uint32_t v){  
        vector<uint32_t> q;
        size_t pointer = 0;
        uint32_t curr_node = u;
        visited_bfs[curr_node] = true;
        q.push_back(curr_node); 
        while (pointer < q.size()) {
            curr_node = q[pointer];
            pointer++;
            for (const uint32_t i : out_edge[curr_node]){
                if (!visited_bfs[i]){
                    visited_bfs[i] = true;
                    q.push_back(i);
                }
            }
        }
        bool ans = visited_bfs[v];
        for (const uint32_t i : q){
            visited_bfs[i] = false;
        }
        return ans;
    }
    
    bool answer_query(const uint32_t u, const uint32_t v){
        return calculate_bfs(u, v);                
    }
private:
    vector<bool> visited_bfs; 
};

class Dfs : public Algorithms{

public:
    Dfs(const Setting& setting_, Logger& logg_) : Algorithms(setting_, logg_){
        visited_dfs.resize(setting.nodes, false);
    }

    bool calculate_dfs(const uint32_t u, const uint32_t v){ 
        visited_dfs[u] = true;
        if (visited_dfs[v]){
            return true;
        }
        for (const auto i : out_edge[u]){
            if (!(visited_dfs[i])){
                calculate_dfs(i, v);
            }
        }
        return visited_dfs[v];
    }
    bool answer_query(const uint32_t u, const uint32_t v){       
        visited_dfs.assign(visited_dfs.size(), false);
        return calculate_dfs(u, v);
    }
private:
    vector<bool> visited_dfs;
};

class Bibfs : public Algorithms{

public:
    Bibfs(const Setting& setting_, Logger& logg_) : Algorithms(setting_, logg_){
        visited_bibfs_source.resize(setting.nodes, false);
        visited_bibfs_sink.resize(setting.nodes, false);
    }
    bool answer_query(const uint32_t u, const uint32_t v){
        return calculate_bibfs(u, v);                
    }
private:
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;
    
    bool calculate_bibfs(const uint32_t u, const uint32_t v){  
        bool found_path = false;
        uint32_t curr_node;
        vector<uint32_t> source_queue, sink_queue;
        size_t source_pointer = 0;
        size_t sink_pointer = 0;
        visited_bibfs_source[u] = true;
        visited_bibfs_sink[v] = true;
        source_queue.push_back(u);
        sink_queue.push_back(v);
        while (!found_path && source_pointer < source_queue.size()
                            && sink_pointer < sink_queue.size()) {
            //running bfs for the source queue one time
            curr_node = source_queue[source_pointer];
            source_pointer++;
            for (const auto i : out_edge[curr_node]){
                if (!visited_bibfs_source[i]){
                    visited_bibfs_source[i] = true;
                    source_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            curr_node = sink_queue[sink_pointer];
            sink_pointer++;
            for (const auto i : in_edge[curr_node]){
                if (!visited_bibfs_sink[i]){
                    visited_bibfs_sink[i] = true;
                    sink_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
        }
        for (const uint32_t i : source_queue){
            visited_bibfs_source[i] = false;
        }
        for (const uint32_t i : sink_queue){
            visited_bibfs_sink[i] = false;
        }
        return found_path;
    }
};

class Sv : public Algorithms{

public:
    Sv(size_t count_, const Setting& setting_, Logger& logg_, uint32_t sv_seed_) : Algorithms(setting_, logg_), count(count_){
        visited_bibfs_source.resize(setting.nodes, false);
        visited_bibfs_sink.resize(setting.nodes, false);
        sv_seed = sv_seed_;
    }
    virtual ~Sv(){}
    
    bool calculate_sv(const uint32_t u, const uint32_t v){  
        // instead of searching, we can use hash map as well.
        // since |sv| is little, i guess it's better to simply search
        auto it = find_sv(u);
        if (it != reachability_tree.end()){
            return (*it)->reaches(v);
        }
        it = find_sv(v);
        if (it != reachability_tree.end()){
            return (*it)->is_reachable_from(u);
        }
        for (const auto &rt: reachability_tree){
            //obs. 1
            if (rt->is_reachable_from(u) && rt->reaches(v)){
                return true;
            }
            //obs. 2
            if (rt->reaches(u) && !rt->reaches(v)){
                return false;
            }
            //obs. 3
            if (rt->is_reachable_from(v) && !rt->is_reachable_from(u)){
                return false;
            }
        }
        
        return calculate_bibfs(u, v);
    }
    bool answer_query (const uint32_t u, const uint32_t v){
        if (!generated_sv){
            generate_sv_list();
            generated_sv = true;
        }
        return calculate_sv(u, v);                
    }
private:
    vector<unique_ptr<reachabilityTree>> reachability_tree;
    //bringing bibfs fallback algorithm inside sv (because we need the same graph out/in edges)
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;
    size_t count;
    uint32_t sv_seed;
    bool generated_sv = false;
    
    const vector<unique_ptr<reachabilityTree>>::iterator find_sv(uint32_t sv){
        return find_if(reachability_tree.begin(), reachability_tree.end(),
                        [&sv](const unique_ptr<reachabilityTree>& obj) {return (*obj).id == sv;});
    }

    void generate_candidate_svs(vector<uint32_t>& non_isolated_svs, vector<uint32_t>& half_isolated_svs){
        for (size_t i = 0; i < setting.nodes; i++){
            if (is_non_isolate(i))
                non_isolated_svs.push_back(i);
            else if (is_half_isolate(i))
                half_isolated_svs.push_back(i);
        }
    }

    void generate_sv_list(){
        // random_device os_seed;
        engine generator(sv_seed);
        vector <uint32_t> non_isolated_svs, half_isolated_svs;
        generate_candidate_svs(non_isolated_svs, half_isolated_svs); 

        size_t i = 0;
        logg.algorithm += '{';
        while (i < count){
            uint32_t sv = 0;

            if (non_isolated_svs.size() > 0){
                uniform_int_distribution< u32 > non_distribute(0, non_isolated_svs.size() - 1);
                size_t index = non_distribute(generator);
                sv = non_isolated_svs[index];
                non_isolated_svs.erase(non_isolated_svs.begin() + index);
            }
            else if (half_isolated_svs.size() > 0){
                uniform_int_distribution< u32 > half_distribute(0, half_isolated_svs.size() - 1);
                size_t index = half_distribute(generator);
                sv = half_isolated_svs[index];
                half_isolated_svs.erase(half_isolated_svs.begin() + index);
            }
            else{
                uniform_int_distribution< u32 > distribute(0, setting.nodes - 1);
                sv = distribute(generator);
            }

            if (find_sv(sv) != reachability_tree.end()) {
                continue;
            }
            // reachability_tree[sv] = new reachabilityTree(sv);
            reachability_tree.push_back(unique_ptr<reachabilityTree>(new reachabilityTree(sv, setting.nodes, out_edge, in_edge)));
            // reachability_tree.push_back(make_unique<reachabilityTree>(sv));
            logg.algorithm += to_string(sv) + ", ";
            i++;
        }
        logg.algorithm += '}';
    }

    bool is_non_isolate(const uint32_t u){
        return (out_edge[u].size() != 0 && in_edge[u].size() != 0);
    }

    bool is_half_isolate(const uint32_t u){
        return (out_edge[u].size() != 0 || in_edge[u].size() != 0);
    }
    
    void update_sv(const uint32_t u, const uint32_t v){
        for (const auto &rt : reachability_tree){
            rt->update(u, v, out_edge, in_edge);
        }
    }

    void add_edge(const uint32_t u, const uint32_t v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
        if (generated_sv)
            update_sv(u, v);
    }

    bool calculate_bibfs(const uint32_t u, const uint32_t v){  
        bool found_path = false;
        uint32_t curr_node;
        vector<uint32_t> source_queue, sink_queue;
        size_t source_pointer = 0;
        size_t sink_pointer = 0;
        visited_bibfs_source[u] = true;
        visited_bibfs_sink[v] = true;
        source_queue.push_back(u);
        sink_queue.push_back(v);
        while (!found_path && source_pointer < source_queue.size()
                            && sink_pointer < sink_queue.size()) {
            //running bfs for the source queue one time
            curr_node = source_queue[source_pointer];
            source_pointer++;
            for (const auto i : out_edge[curr_node]){
                if (!visited_bibfs_source[i]){
                    visited_bibfs_source[i] = true;
                    source_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            curr_node = sink_queue[sink_pointer];
            sink_pointer++;
            for (const auto i : in_edge[curr_node]){
                if (!visited_bibfs_sink[i]){
                    visited_bibfs_sink[i] = true;
                    sink_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
        }
        for (const uint32_t i : source_queue){
            visited_bibfs_source[i] = false;
        }
        for (const uint32_t i : sink_queue){
            visited_bibfs_sink[i] = false;
        }
        return found_path;
    }
};


void set_time(string& t){
    auto timepoint = chrono::system_clock::now();
    auto coarse = chrono::system_clock::to_time_t(timepoint);
    auto fine = chrono::time_point_cast<std::chrono::milliseconds>(timepoint);

    char buffer[sizeof "9999-12-31 23:59:59.999"];
    std::snprintf(buffer + std::strftime(buffer, sizeof buffer - 3, 
                    "%F %T.", std::localtime(&coarse)), 4, "%03lu",
                    fine.time_since_epoch().count() % 1000);
    t = buffer;
}



class Program{
public:
    Program(){}
    
    void read_parse_input(const int argc, const char* argv[]){
        for (int i = 1; i < argc; i += 2){
            if (!strcmp(argv[i], "-alg")){
                if (strcmp(argv[i+1], "dfs") && strcmp(argv[i+1], "bfs") && 
                    strcmp(argv[i+1], "bibfs") && string(argv[i+1]).substr(0, 3) != "sv_") {
                    cerr << "Wrong Input for Algorithm.\n";
                    exit(0);
                }
                setting.ALGORITHM = argv[i+1];
            }
            if (!strcmp(argv[i], "-qp")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1]) || 
                    stoi(argv[i+1]) < 0 || stoi(argv[i+1]) > 100){
                    cerr << "Wrong Input for Query Percentage.\n";
                    exit(0);
                }
                setting.QUERY_PERCENTAGE = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-trc")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1])){
                    cerr << "Wrong Input for Test Run Count.\n";
                    exit(0);
                }
                setting.TEST_RUN_COUNT = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-ts")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1])){
                    cerr << "Wrong Input for Timeout Seconds.\n";
                    exit(0);
                }
                setting.TIMEOUT_SEC = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-os")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1])){
                    cerr << "Wrong Input for Operation Seed.\n";
                    exit(0);
                }
                setting.OPERATION_SEED = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-qs")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1])){
                    cerr << "Wrong Input for Query Seed.\n";
                    exit(0);
                }
                setting.QUERY_SEED = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-qt")){
                if (strspn(argv[i+1], "-.0123456789" ) != strlen(argv[i+1])){
                    cerr << "Wrong Input for Query Timestamp.\n";
                    exit(0);
                }
                setting.QUERY_TIMESTAMP = stoi(argv[i+1]);
            }
            if (!strcmp(argv[i], "-inp")){
                std::ifstream infile(argv[i+1]);
                if (!infile.good()){
                    cerr << "Input File Does Not Exist.\n";
                    exit(0);
                }
                setting.INPUT_FILE = argv[i+1];
            }
            if (!strcmp(argv[i], "-meta")){
                std::ifstream infile(argv[i+1]);
                if (!infile.good()){
                    cerr << "Meta File Does Not Exist.\n";
                    exit(0);
                }
                setting.META_FILE = argv[i+1];
            }
            if (!strcmp(argv[i], "-out")){
                std::ifstream infile(argv[i+1]);
                if (!infile.good()){
                    cerr << "Warning: Output File Did Not Exist. One Will Be Created.\n";
                }
                setting.OUTPUT_FILE = argv[i+1];
            }
            if (!strcmp(argv[i], "-log")){
                std::ifstream infile(argv[i+1]);
                if (!infile.good()){
                    cerr << "Warning: Log File Did Not Exist. One Will Be Created.\n";
                }
                setting.LOG_FILE = argv[i+1];
            }
        }

        read_input_file(); //read dataset and build input
        operations.reserve(setting.input_lines * ((100 + setting.QUERY_PERCENTAGE) / 100));
        generate_operations(); //add query operatios to input and store result in "operations"

        logg.test_id = get_test_id();
        logg.algorithm = setting.ALGORITHM; 
    }

    void execute_test(){
        set_time(logg.start_time);
        for (size_t i = 0; i < setting.TEST_RUN_COUNT; i++){
            unique_ptr<Algorithms> alg; 
            if (setting.ALGORITHM == "dfs")
                alg = unique_ptr<Algorithms>(new Dfs(setting, logg));
            else if (setting.ALGORITHM == "bfs")
                alg = unique_ptr<Algorithms>(new Bfs(setting, logg));
            else if (setting.ALGORITHM == "bibfs")
                alg = unique_ptr<Algorithms>(new Bibfs(setting, logg));
            else if (setting.ALGORITHM.substr(0, 3) == "sv_"){
                alg = unique_ptr<Algorithms>(new Sv(stoi(setting.ALGORITHM.substr(3)), setting, logg, i));
            }
            else
                assert(false);
            alg->run(operations);
            logg.run_durations.push_back(logg.query_durations.back() + logg.insertion_durations.back());
        }
        
        set_time(logg.end_time);
    }


    void write_to_log(){
        double_t duration = sum(logg.run_durations);
        double_t queries = sum(logg.query_durations);
        double_t insertions = sum(logg.insertion_durations);
        
        ofstream log_file;
        log_file.open(setting.LOG_FILE, std::ios_base::app);
        log_file << "test id: " << logg.test_id << '\n' <<
                    "input file name: " << setting.INPUT_FILE << '\n' <<
                    "input file lines: " << setting.input_lines << '\n' <<
                    "seed: " << setting.OPERATION_SEED << ", " << setting.QUERY_SEED << '\n' <<
                    "#run: " << setting.TEST_RUN_COUNT << '\n' << 
                    "algorithm: " << logg.algorithm << '\n' <<
                    "#query opertions: " << logg.query_operations_cnt << '\n' <<
                    "#insertion operations: " << logg.insertion_operations_cnt << '\n' <<
                    "start time: " << logg.start_time << '\n' <<
                    "end time: " << logg.end_time << '\n' <<
                    "duration: " << duration << "{";
                    for (const auto x : logg.run_durations){
                        log_file << x << " ";
                    }
                    log_file << '}' << ' ' << '[' << variance(logg.run_durations) << ']' << '\n' <<
                    "queries: " << queries << "{";
                    for (const auto x : logg.query_durations){
                        log_file << x << " ";
                    }
                    log_file << '}' << ' ' << '[' << variance(logg.query_durations) << ']' << '\n' <<
                    "insertions: " << insertions << "{";
                    for (const auto x : logg.insertion_durations){
                        log_file << x << " ";
                    }
                    log_file << '}' << ' ' << '[' << variance(logg.insertion_durations) << ']' << '\n' <<
                    "hashed output: " << logg.hashed_output << '\n' <<
                    "#reachable queries: " << logg.num_reachable_queries << '\n' <<
                    string(50, '*') << "\n";
    }


private:
    Setting setting;
    Logger logg;
    vector<Operation> insertion_operations;
    vector<Operation> operations;

    void convert_input(){
        ifstream file(setting.INPUT_FILE);
        std::filesystem::path p(setting.INPUT_FILE);
        if (filesystem::exists("build/" + string(p.filename()))){
            cout << "Directory already exists." << endl;
        }
        else{
            filesystem::create_directories("build");
        }

        ofstream new_file("build/" + string(p.filename()) , std::ios::trunc);

        uint32_t max_nodes = 0;
        uint32_t lines = 0;
       
        uint32_t u, v;
        string type, timestamp;
        while (file >> u >> v >> type >> timestamp){
            
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
    void read_input_file(){
        convert_input(); 
        cout << "Input generated!" << endl;

        std::filesystem::path p(setting.INPUT_FILE);
        ifstream infile("build/" + string(p.filename()));
        uint32_t u, v;
        int timestamp;

        while (infile >> u >> v >> timestamp){
            insertion_operations.push_back(Operation(false, make_pair(u, v), timestamp));
        }
        infile.close();
    }
    void generate_operations (){
        engine generator(setting.OPERATION_SEED);
        engine query_generator(setting.QUERY_SEED);
        uniform_int_distribution< u32 > distribute(0, 99);
        uniform_int_distribution< u32 > query_chance_distribute(0, setting.nodes-1);
        uint32_t u, v;
        for (const auto &x : insertion_operations){
            u = x.arguments.first;
            v = x.arguments.second;
            // currently each insertion and query have the same timestamp
            // we can modify for different implementations later 
            if (x.timestamp > setting.QUERY_TIMESTAMP && distribute(generator) < setting.QUERY_PERCENTAGE){
                uint32_t u_q = query_chance_distribute(query_generator);
                uint32_t v_q = query_chance_distribute(query_generator);
                operations.push_back(Operation(true, make_pair(u_q, v_q), x.timestamp));
                logg.query_operations_cnt ++;
            }
            operations.push_back(Operation(false, make_pair(u, v), x.timestamp));
        }
        logg.insertion_operations_cnt = insertion_operations.size();
    }
    int get_test_id(){ //can probably read backwards later to enhance speed
        string line;
        ifstream log_file (setting.LOG_FILE);
        int id = 0;
        while(getline(log_file, line)) { 
            if (line.find('*') != string::npos) {
                id ++;
            }
        }
        return id;
    }

};

int main(int argc, const char* argv[]){

    Program program; 
    program.read_parse_input(argc, argv);

    program.execute_test();

    program.write_to_log();

}