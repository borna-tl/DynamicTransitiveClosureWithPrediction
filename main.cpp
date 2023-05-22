#include <bits/stdc++.h>
#include <fstream>
#include <queue>
#include <string>
#include <sstream> 

using namespace std;
using u32    = uint_least32_t; 
using engine = std::mt19937;


#define PROGRESS_STAMP 10 //define the progress bar count
#define PBSTR "++++++++++++++++++++++++++++++++++++++++++++++++++"
#define PBWIDTH 50


struct Settings{
    Settings(){}

    vector<pair<uint32_t, uint32_t>> input_file_operations;
    uint32_t nodes = 0;
    int input_num_lines = 0;
    string INPUT_FILE = "sample.txt";
    string META_FILE = "meta-sample.txt";
    string OUTPUT_FILE = "output.txt";
    string LOG_FILE = "log.txt";
    uint32_t QUERY_PERCENTAGE = 50;
    uint32_t TEST_RUN_COUNT = 10;
    uint32_t TIMEOUT_SEC = 1800;
    uint32_t OPERATION_SEED = 1223;
    uint32_t QUERY_SEED = 2334;
} setting;


void print_progress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


struct Logger {
    int test_id;
    string input_file_name = setting.INPUT_FILE;
    pair <int, int> seed = make_pair(setting.OPERATION_SEED, setting.QUERY_SEED);
    int num_queries = -1;
    int num_insertions = -1;
    string algorithm; 
    string start_time;
    string end_time;
    vector<int64_t> run_duration; //in mili-seconds
    size_t hashed_output;
    int num_reachable_queries = -1;
} logg;

struct Operation {
    Operation(bool is_query_, pair <uint32_t, uint32_t> arguments_) :
                is_query(is_query_), arguments(arguments_){}
    bool is_query;
    pair <uint32_t, uint32_t> arguments;
};

class reachabilityTree{ //this is a simple incremental reachability tree for vertex s
public:
    reachabilityTree(){}
    reachabilityTree(uint32_t id_):id(id_){
        if (setting.nodes == 0){
            //look for some better form of sending errors
            cout << "Not expecting zero nodes" << endl;
            exit(0);
        }
        r_plus.resize(setting.nodes);
        r_minus.resize(setting.nodes);
        r_plus[id] = true;
        r_minus[id] = true;
    }
    ~reachabilityTree(){}

    void update(uint32_t u, uint32_t v, const vector<vector<uint32_t>> &out_edge, const vector<vector<uint32_t>> &in_edge){
        update_reachability(u, v, out_edge, r_plus); //source reachability
        update_reachability(v, u, in_edge, r_minus); //sink reachability
    }
    
    void print_reachability_list(){
        cout << "Reachable Nodes from " << id << ": ";
        for (size_t i = 0; i < setting.nodes; i++)
            if (r_plus[i])
                cout << i << " ";
        cout << endl;
        cout << "Reachable Nodes to " << id << ": ";
        for (size_t i = 0; i < setting.nodes; i++)
            if (r_minus[i])
                cout << i << " ";
        cout << endl;
    }
    bool reaches(uint32_t u){
        return r_plus[u];
    }
    bool is_reachable_from(uint32_t u){
        return r_minus[u];
    }
    uint32_t id; 
private:
    vector <bool> r_plus;
    vector <bool> r_minus;

    void update_reachability(uint32_t u, uint32_t v, const vector<vector<uint32_t>> &edge, vector<bool>& r){ //read on refrence vs. pointers
        if (r[v])
            return;
        if (!r[u])
            return;
        queue<uint32_t> q;
        r[v] = true;
        q.push(v);
        while (!q.empty()) {
            v = q.front();
            q.pop();
            for (auto i : edge[v]){
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
    Algorithms(){
        out_edge.assign(setting.nodes, vector<uint32_t>());
        in_edge.assign(setting.nodes, vector<uint32_t>());
    }
    virtual ~Algorithms(){};
    virtual bool answer_query(uint32_t u, uint32_t v) = 0; 

    void run(const vector<Operation> operations){
        uint32_t u, v;
        clock_t tStart = clock();
        int queries_answered = 0, true_q = 0, num_insertions = 0;
        size_t c_out = 0;
        for (auto x : operations){
            u = x.arguments.first;
            v = x.arguments.second;
            if (x.is_query){
                bool result = answer_query(u, v);   
                queries_answered ++;          
                true_q += (result == true); 
                results.push_back(result);
            }
            else{
                add_edge(u, v);
                num_insertions++;
            }
            c_out++;
            if (c_out > PROGRESS_STAMP * operations.size() / 100){
                c_out = 0;
                if ((double)(clock() - tStart)/CLOCKS_PER_SEC > setting.TIMEOUT_SEC){
                    cerr << "Timeout!" << endl;
                    exit(0);
                }
                // print_progress((double)(queries_answered+num_insertions)/operations.size());
            }
        }
        // cout << endl;
        // printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        // cout << "Queries answered: " << queries_answered << endl;
        // cout << "Reachable queries: " << true_q << endl;
        // ofstream outfile(OUTPUT_FILE, ios_base::app);
        // ostream_iterator<string> output_iterator(outfile, "");
        size_t hash_output = hash<vector<bool>>{}(results);
        // outfile << hash<vector<bool>>{}(results) << "\n";
        // outfile.close();
        logg.num_queries = queries_answered;
        logg.num_insertions = num_insertions;
        logg.hashed_output = hash_output;
        logg.num_reachable_queries = true_q;
    }
protected:
    vector<vector<uint32_t>> out_edge;
    vector<vector<uint32_t>> in_edge;
    vector <bool> results;
    virtual void add_edge(uint32_t u, uint32_t v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
    }
    
};


class Bfs : public Algorithms{

public:
    Bfs() : Algorithms(){
        visited_bfs.resize(setting.nodes, false);
    }
    bool calculate_bfs(uint32_t u, uint32_t v){  
        vector<uint32_t> q;
        size_t pointer = 0;
        visited_bfs[u] = true;
        q.push_back(u); 
        while (pointer < q.size()) {
            u = q[pointer];
            pointer++;
            for (uint32_t i : out_edge[u]){
                if (!visited_bfs[i]){
                    visited_bfs[i] = true;
                    q.push_back(i);
                }
            }
        }
        bool ans = visited_bfs[v];
        for (uint32_t i : q){
            visited_bfs[i] = false;
        }
        return ans;
    }
    bool answer_query(uint32_t u, uint32_t v){
        return calculate_bfs(u, v);                
    }
private:
    vector<bool> visited_bfs; 
};

class Dfs : public Algorithms{

public:
    Dfs() : Algorithms(){
        visited_dfs.resize(setting.nodes, false);
    }
    bool calculate_dfs(uint32_t u, uint32_t v){ 
        visited_dfs[u] = true;
        if (visited_dfs[v]){
            return true;
        }
        for (auto i : out_edge[u]){
            if (!(visited_dfs[i])){
                calculate_dfs(i, v);
            }
        }
        return visited_dfs[v];
    }
    bool answer_query(uint32_t u, uint32_t v){       
        visited_dfs.assign(visited_dfs.size(), false);
        return calculate_dfs(u, v);
    }
private:
    vector<bool> visited_dfs;
};

class Bibfs : public Algorithms{

public:
    Bibfs() : Algorithms(){
        visited_bibfs_source.resize(setting.nodes, false);
        visited_bibfs_sink.resize(setting.nodes, false);
    }
    bool answer_query(uint32_t u, uint32_t v){
        return calculate_bibfs(u, v);                
    }
private:
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;
    bool calculate_bibfs(uint32_t u, uint32_t v){  
        bool found_path = false;
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
            u = source_queue[source_pointer];
            source_pointer++;
            for (auto i : out_edge[u]){
                if (!visited_bibfs_source[i]){
                    visited_bibfs_source[i] = true;
                    source_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            v = sink_queue[sink_pointer];
            sink_pointer++;
            for (auto i : in_edge[v]){
                if (!visited_bibfs_sink[i]){
                    visited_bibfs_sink[i] = true;
                    sink_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
        }
        for (uint32_t i : source_queue){
            visited_bibfs_source[i] = false;
        }
        for (uint32_t i : sink_queue){
            visited_bibfs_sink[i] = false;
        }
        return found_path;
    }
};

class Sv : public Algorithms{

public:
    Sv(int count_) : Algorithms(), count(count_){
        visited_bibfs_source.resize(setting.nodes, false);
        visited_bibfs_sink.resize(setting.nodes, false);
        generate_sv_list();
    }
    virtual ~Sv(){
        // for (auto sv : sv_list){
        //     delete reachability_tree[sv];
        // }
    }
    bool calculate_sv(uint32_t u, uint32_t v){  
        // instead of searching, we can use hash map as well.
        // since |sv| is little, i guess it's better to simply search
        // if (find(sv_list.begin(), sv_list.end(), u) != sv_list.end()){
        auto it = find_sv(u);
        if (it != reachability_tree.end()){
            return (*it)->reaches(v);
            // return reachability_tree[u]->reaches(v);
        }
        it = find_sv(v);
        if (it != reachability_tree.end()){
            return (*it)->is_reachable_from(u);
        }
        // if (find(sv_list.begin(), sv_list.end(), v) != sv_list.end()){
            // return reachability_tree[v]->is_reachable_from(u);
        // }
        for (auto &rt: reachability_tree){
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
    bool answer_query (uint32_t u, uint32_t v){
        return calculate_sv(u, v);                
    }
private:
    //decided to use normal pointers, can use smart pointers later 
    vector<unique_ptr<reachabilityTree>> reachability_tree;
    // reachabilityTree* reachability_tree[MAX_NODES];
    // vector <uint32_t> sv_list;
    //bringing bibfs fallback algorithm inside sv (because we need the same graph out/in edges)
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;
    int count;
    const vector<unique_ptr<reachabilityTree>>::iterator find_sv(uint32_t sv){ //check return type later
        return find_if(reachability_tree.begin(), reachability_tree.end(),
                        [&sv](const unique_ptr<reachabilityTree>& obj) {return (*obj).id == sv;});
    }
    void generate_sv_list(){
        vector<uint32_t> svs = {0,2,8,9,10,12,16,17,19,29}; //from katherin's work for answers
        // vector<uint32_t> svs = {0,1,2,4,5,6,7,8,12,16}; //for bio-protein
        // vector<uint32_t> svs = {0,1,2,5,10,12,16,17,19,20}; //for blog-nat05-6m
        // vector<uint32_t> svs = {0,2,4,5,8,9,10,16,17,24}; //for ca-dblp
        // vector<uint32_t> svs = {0,2,3,8,11,16,17,19,20,28}; //for gnutella-25
        // program crashes
        // vector<uint32_t> svs = {0,1,2,4,5,8,16,18,20,24}; //for email-inside

        // sv_list.clear();

        random_device os_seed;
        //generate const seed instead of totally rand 
        engine generator(os_seed());
        uniform_int_distribution< u32 > distribute(0, svs.size() - 1);

        logg.algorithm += '{';
        int i = 0;
        while (i < count){
            uint32_t sv = svs[distribute(generator)];
            if (find_sv(sv) != reachability_tree.end()) 
                continue;
            // reachability_tree[sv] = new reachabilityTree(sv);
            reachability_tree.push_back(unique_ptr<reachabilityTree>(new reachabilityTree(sv)));
            // reachability_tree.push_back(make_unique<reachabilityTree>(sv));
            
            // sv_list.push_back(sv);
            logg.algorithm += to_string(sv) + ", ";
            i++;
        }
        logg.algorithm += '}';

    }
    void update_sv(uint32_t u, uint32_t v){
        for (auto &rt : reachability_tree){
            rt->update(u, v, out_edge, in_edge);
        }
        // for (int i = 0; i < reachability_tree.size(); ++i){
            // rt->update(u, v, out_edge, in_edge);
            // (*reachability_tree[i]).update(u, v, out_edge, in_edge);
        // }
    }
    void add_edge(uint32_t u, uint32_t v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
        update_sv(u, v);
    }
    bool calculate_bibfs(uint32_t u, uint32_t v){  
        bool found_path = false;
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
            u = source_queue[source_pointer];
            source_pointer++;
            for (auto i : out_edge[u]){
                if (!visited_bibfs_source[i]){
                    visited_bibfs_source[i] = true;
                    source_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            v = sink_queue[sink_pointer];
            sink_pointer++;
            for (auto i : in_edge[v]){
                if (!visited_bibfs_sink[i]){
                    visited_bibfs_sink[i] = true;
                    sink_queue.push_back(i);
                }
                if (visited_bibfs_source[i] && visited_bibfs_sink[i]){
                    found_path = true;
                }
            }
        }
        for (uint32_t i : source_queue){
            visited_bibfs_source[i] = false;
        }
        for (uint32_t i : sink_queue){
            visited_bibfs_sink[i] = false;
        }
        return found_path;
    }
};

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

void write_to_log(){
    ofstream log_file;
    log_file.open(setting.LOG_FILE, std::ios_base::app);
    log_file << "test id: " << logg.test_id << '\n' <<
                "input file name: " << logg.input_file_name << '\n' <<
                "seed: " << logg.seed.first << ", " << logg.seed.second << '\n' <<
                "#run: " << setting.TEST_RUN_COUNT << '\n' << 
                "#queries: " << logg.num_queries << '\n' <<
                "#insertions: " << logg.num_insertions << '\n' <<
                "algorithm: " << logg.algorithm << '\n' <<
                "start time: " << logg.start_time << '\n' <<
                "end time: " << logg.end_time << '\n' <<
                "duration: " << std::accumulate(logg.run_duration.begin(),
                                    logg.run_duration.end(), decltype(logg.run_duration)::
                                    value_type(0)) << "{";
                for (auto x : logg.run_duration)
                    log_file << x << " ";
                log_file << '}' << '\n' <<
                "hashed output: " << logg.hashed_output << '\n' <<
                "#reachable queries: " << logg.num_reachable_queries << '\n' <<
                string(50, '*') << "\n";
}

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

void read_meta_file(){
    ifstream infile(setting.META_FILE);
    string command;
    infile >> command >> setting.nodes;
    infile >> command >> setting.input_num_lines;
    infile >> command >> logg.input_file_name;
    setting.nodes++; //because input file is zero-based
    infile.close();
}

void read_input_file(){
    ifstream infile(setting.INPUT_FILE);
    uint32_t u, v;
    while (infile >> u >> v){
        setting.input_file_operations.push_back(make_pair(u, v));
    }
    infile.close();
}

void generate_operations (vector<Operation>& operations){
    // random_device os_seed; //can use later for seeding engine.
    engine generator(setting.OPERATION_SEED);
    engine query_generator(setting.QUERY_SEED);
    uniform_int_distribution< u32 > distribute(0, 99);
    uniform_int_distribution< u32 > query_chance_distribute(0, setting.nodes-1);
    uint32_t u, v;

    for (auto x : setting.input_file_operations){
        u = x.first;
        v = x.second;
        while (distribute(generator) < setting.QUERY_PERCENTAGE){
            uint32_t u_q = query_chance_distribute(query_generator);
            uint32_t v_q = query_chance_distribute(query_generator);
            operations.push_back(Operation(true, make_pair(u_q, v_q)));
        }
        operations.push_back(Operation(false, make_pair(u, v)));
    }
}

void execute_test(const vector<Operation>& operations, string mode){
    
    for (size_t i = 0; i < setting.TEST_RUN_COUNT; i++){
        unique_ptr<Algorithms> alg;  // use smart pointer
        if (mode == "dfs")
            alg = unique_ptr<Algorithms>(new Dfs());
        else if (mode == "bfs")
            alg = unique_ptr<Algorithms>(new Bfs());
        else if (mode == "bibfs")
            alg = unique_ptr<Algorithms>(new Bibfs());
        else if (mode == "sv_1")
            alg = unique_ptr<Algorithms>(new Sv(1));
        else if (mode == "sv_2")
            alg = unique_ptr<Algorithms>(new Sv(2));
        else
            assert(false);
        auto started = std::chrono::high_resolution_clock::now();
        alg->run(operations);
        logg.run_duration.push_back(chrono::duration_cast<std::chrono::milliseconds>
                         (chrono::high_resolution_clock::now()-started).count());
    }
   
}


int main(int argc, char* argv[]){

    
    if (strcmp(argv[1], "dfs") && strcmp(argv[1], "bfs") && 
        strcmp(argv[1], "bibfs") && strcmp(argv[1], "sv_1") && strcmp(argv[1], "sv_2")){
            cerr << "Wrong Input\n";
            exit(0);
    }
    
    read_meta_file();
    vector<Operation> operations;
    operations.reserve(setting.input_num_lines * ((100 + setting.QUERY_PERCENTAGE) / 100));
    read_input_file(); //read and store input file in "input_file_operations"
    generate_operations(operations); //add query operatios to "input_file_operations" and store result in "operations"

    logg.test_id = get_test_id();
    logg.algorithm = argv[1]; 
    set_time(logg.start_time);

    execute_test(operations, argv[1]);
    
    set_time(logg.end_time);
    write_to_log();

}