#include <bits/stdc++.h>
#include <fstream>
#include <queue>

using namespace std;
using u32    = uint_least32_t; 
using engine = std::mt19937;

#define MAX_NODES 900000+7

#define INPUT_FILE "sample.txt"
#define META_FILE "meta-sample.txt"


class reachabilityTree{ //this is a simple incremental reachability tree for vertex s
public:
    reachabilityTree(){}
    reachabilityTree(int id_):id(id_){
        //should run bfs to initialize: do we consider the initialization in the timing as well?
        r_plus[id] = true;
        r_minus[id] = true;
    }

    void update(int u, int v, const vector<vector<int>> &out_edge, const vector<vector<int>> &in_edge){
        update_reachability(u, v, out_edge, r_plus); //source reachability
        update_reachability(v, u, in_edge, r_minus); //sink reachability
    }
    void update_reachability(int u, int v, const vector<vector<int>> &edge, bool* r){ //read on refrence vs. pointers
        if (r[v])
            return;
        if (!r[u])
            return;
        queue<int> q;
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
    
    void print_reachability_list(){
        cout << "Reachable Nodes from " << id << ": ";
        for (int i = 0; i < MAX_NODES; i++)
            if (r_plus[i])
                cout << i << " ";
        cout << endl;
        cout << "Reachable Nodes to " << id << ": ";
        for (int i = 0; i < MAX_NODES; i++)
            if (r_minus[i])
                cout << i << " ";
        cout << endl;
    }
    bool reaches(int u){
        return r_plus[u];
    }
    bool is_reachable_from(int u){
        return r_minus[u];
    }
private:
    int id;

    //for u in [MAX_NODES]:
    bool r_plus[MAX_NODES]; // is node u reachable from s? s--->?--->u 
    bool r_minus[MAX_NODES]; // is node s reachable from u? u--->?--->s 
};

class Algorithms{
public:
    Algorithms(){
        ifstream infile(META_FILE);
        infile >> nodes;
        nodes++; //because input file is zero-based
        out_edge.assign(nodes, vector<int>());
        in_edge.assign(nodes, vector<int>());
        infile.close();
    }
    virtual bool answer_query(int u, int v) = 0; 
    void run(){
        ifstream infile(INPUT_FILE);
        
        // random_device os_seed; //can use later for seeding engine.
        engine generator(2334);
        engine query_generator(1223);
        uniform_int_distribution< u32 > distribute(0, 99);
        uniform_int_distribution< u32 > query_chance_distribute(0, nodes-1);

        int u, v;
        clock_t tStart = clock();
        int queries_answered = 0, true_q = 0;
        while (infile >> u >> v){
            while (distribute(generator) < 33){
                int u_q = query_chance_distribute(query_generator);
                int v_q = query_chance_distribute(query_generator);
                bool result = answer_query(u_q, v_q);   
                queries_answered ++;          
                true_q += (result == true); 
                //TBD: add correctness checking
            }
            add_edge(u, v);
        }

        printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        cout << "Queries answered: " << queries_answered << endl;
        cout << "Reachable queries: " << true_q << endl;
        infile.close();
    }
protected:
    int nodes;
    vector<vector<int>> out_edge;
    vector<vector<int>> in_edge;
    virtual void add_edge(int u, int v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
    }
};


class Bfs : public Algorithms{

public:
    Bfs() : Algorithms(){
        visited_bfs.resize(nodes, false);
    }
    bool calculate_bfs(int u, int v){  
        vector<int> q;
        size_t pointer = 0;
        visited_bfs[u] = true;
        q.push_back(u); 
        while (pointer < q.size()) {
            u = q[pointer];
            pointer++;
            for (int i : out_edge[u]){
                if (!visited_bfs[i]){
                    visited_bfs[i] = true;
                    q.push_back(i);
                }
            }
        }
        bool ans = visited_bfs[v];
        for (int i : q){
            visited_bfs[i] = false;
        }
        return ans;
    }
    bool answer_query(int u, int v){
        return calculate_bfs(u, v);                
    }
private:
    vector<bool> visited_bfs; 
};

class Dfs : public Algorithms{

public:
    Dfs() : Algorithms(){
        visited_dfs.resize(nodes, false);
    }
    bool calculate_dfs(int u, int v){ 
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
    bool answer_query(int u, int v){       
        visited_dfs.assign(visited_dfs.size(), false);
        return calculate_dfs(u, v);
    }
private:
    vector<bool> visited_dfs;
};

class Bibfs : public Algorithms{

public:
    Bibfs() : Algorithms(){
        visited_bibfs_source.resize(nodes, false);
        visited_bibfs_sink.resize(nodes, false);
    }
    virtual bool calculate_bibfs(int u, int v){  
        bool found_path = false;
        vector<int> source_queue, sink_queue;
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
        for (int i : source_queue){
            visited_bibfs_source[i] = false;
        }
        for (int i : sink_queue){
            visited_bibfs_sink[i] = false;
        }
        return found_path;
    }
    bool answer_query(int u, int v){
        return calculate_bibfs(u, v);                
    }
private:
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;

};

class Sv : public Algorithms{

public:
    Sv() : Algorithms(){
        generate_sv_list();
        fallback = new Bibfs();

    }
    bool calculate_sv(int u, int v){  
        //instead of searching, we can use hash map as well.
        // since |sv| is little, i guess it's better to simply search
        if (find(sv_list.begin(), sv_list.end(), u) != sv_list.end()){
            return reachability_tree[u]->reaches(v);
        }
        if (find(sv_list.begin(), sv_list.end(), v) != sv_list.end()){
            return reachability_tree[v]->is_reachable_from(u);
        }
        for (auto sv: sv_list){
            //obs. 1
            if (reachability_tree[sv]->is_reachable_from(u) && reachability_tree[sv]->reaches(v)){
                return true;
            }
            //obs. 2
            if (reachability_tree[sv]->reaches(u) && !reachability_tree[sv]->reaches(v)){
                return false;
            }
            //obs. 3
            if (reachability_tree[sv]->is_reachable_from(v) && !reachability_tree[sv]->is_reachable_from(u)){
                return false;
            }
        }
        //fallback to bfs
        // return Bibfs(input_file, output_file).calculate_bibfs(u, v);
        return fallback->calculate_bibfs(u, v);
    }
    bool answer_query (int u, int v){
        return calculate_sv(u, v);                
    }
private:
    reachabilityTree* reachability_tree[MAX_NODES];
    vector <int> sv_list;
    Bibfs* fallback;
    void generate_sv_list(){
        // vector<int> svs = {0,2,8,16,64,256,2048,8192,32768,65536};
        vector<int> svs = {0};
        for (auto sv : svs){
            reachability_tree[sv] = new reachabilityTree(sv); //using node sv as new SV
            sv_list.push_back(sv);
        }
    }
    void update_sv(int u, int v){
        for (auto sv : sv_list){
            reachability_tree[sv]->update(u, v, out_edge, in_edge);
        }
    }
    void add_edge(int u, int v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
        update_sv(u, v);
    }
};


int main(int argc, char* argv[]){
    if (!strcmp(argv[1], "dfs")){
        Dfs alg;
        alg.run();
    }
    else if (!strcmp(argv[1], "bfs")){
        Bfs alg;
        alg.run();
    }
    else if (!strcmp(argv[1], "bibfs")){
        Bibfs alg;
        alg.run();
    }
    else if (!strcmp(argv[1], "sv")){
        Sv alg;
        alg.run();
    }
    else{
        cerr << "Wrong input" << endl;
    }
}