#include <bits/stdc++.h>
#include <fstream>
#include <queue>

using namespace std;
using u32    = uint_least32_t; 
using engine = std::mt19937;

#define MAX_NODES 900000+7
#define BFS_MODE 1
#define DFS_MODE 2
#define BIBFS_MODE 3
#define SV_MODE 4


class reachabilityTree{//this is a simple incremental reachability tree for vertex s
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
        cout << "Reachable Nodes from s: ";
        for (int i = 0; i < MAX_NODES; i++)
            if (r_plus[i])
                cout << i << " ";
        cout << endl;
        cout << "Reachable Nodes to s: ";
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
    Algorithms(){cout << "Algorithm Default constructor called! Check something...\n";}
    Algorithms(string input_file_, string output_file_):
        input_file(input_file_), output_file(output_file_){
        out_edge.assign(MAX_NODES, vector<int>());
        in_edge.assign(MAX_NODES, vector<int>());
        
    }
    virtual bool answer_query(int u, int v) = 0; 
    void run(){
        ifstream infile(input_file);
        ofstream *file = new ofstream();
        file->open(output_file, ios::ate); 
        
        random_device os_seed;
        engine generator(2334);
        uniform_int_distribution< u32 > distribute(1, 10000);
        
        bool operation;
        int u, v;
        // while (infile >> operation>> u >> v) {
        clock_t tStart = clock();
        int queries_answered = 0, true_q = 0;
        while (infile >> u >> v){
            if (distribute(generator) < 2){
                bool result = answer_query(u, v);   
                queries_answered++;          
                true_q += (result == true); 
                // (*file) << "(" << u << ", " << v << ") is: " << result << endl;
            }
            else{
                add_edge(u, v);
            }
        }
        printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
        cout << "Queries answered: " << queries_answered << endl;
        cout << "Reachable queries: " << true_q << endl;
        file->close();
        infile.close();
    }
protected:
    string input_file, output_file;
    vector<vector<int>> out_edge;
    vector<vector<int>> in_edge;
    virtual void add_edge(int u, int v){
        out_edge[u].push_back(v);
        in_edge[v].push_back(u);
    }
};


class Bfs : public Algorithms{

public:
    Bfs(){"BFS Default constructor called! Check something...\n";}
    Bfs(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        visited_bfs.resize(MAX_NODES, false);
    }
    bool calculate_bfs(int u, int v){  
        vector<int> q;
        int pointer = 0;
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
    vector<bool> visited_bfs; //in order to reduce resizing cost
};

class Dfs : public Algorithms{

public:
    Dfs(){"DFS Default constructor called! Check something...\n\n";}
    Dfs(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        visited_dfs.resize(MAX_NODES, false);
    }
    bool calculate_dfs(int u, int v){  //maybe we should try non-iterative dfs as well
        if (u >= visited_dfs.size() || u < 0 || v < 0 || v >= visited_dfs.size()){
            exit(0);
        }
        visited_dfs[u] = true;
        if (visited_dfs[v]){
            return true;
        }
        for (auto i : out_edge[u]){
            if (i < 0 || i >= visited_dfs.size()){
                exit(0);
            }
            if (!(visited_dfs[i])){
                calculate_dfs(i, v);
            }
        }
        return visited_dfs[v];
    }
    bool answer_query(int u, int v){
        // visited_dfs->clear();
        // fill(visited_dfs->begin(), visited_dfs->end(), 0);
        cout << "vector size is: " << visited_dfs.size() << endl;
        for (auto v: visited_dfs)
            v = false;
        // visited_dfs->assign(visited_dfs->size(), false);
        // visited_dfs.resize(MAX_NODES, false); //maybe not so efficient
        return calculate_dfs(u, v);
    }
private:
    vector<bool> visited_dfs;
};

class Bibfs : public Algorithms{

public:
    Bibfs(){"BiBFS Default constructor called! Check something...\n";}
    Bibfs(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        visited_bibfs_source.resize(MAX_NODES, false);
        visited_bibfs_sink.resize(MAX_NODES, false);
    }
    virtual bool calculate_bibfs(int u, int v){  
        bool found_path = false;
        vector<int> source_queue, sink_queue;
        int source_pointer = 0;
        int sink_pointer = 0;
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
    Sv(){"Sv Default constructor called! Check something...\n";}
    Sv(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        generate_sv_list();
        fallback = new Bibfs("sample.txt", "output.txt");

    }
    bool calculate_sv(int u, int v){  
        // cout << "For insertion on (" << u << ", " << v << "): " << endl;
        // reachability_tree[6]->print_reachability_list(); //for testing

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
        // cout << "sv not successful, falling back..." << endl;
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
        Dfs alg("sample.txt", "output.txt");
        alg.run();
    }
    else if (!strcmp(argv[1], "bfs")){
        Bfs alg("sample.txt", "output.txt");
        alg.run();
    }
    else if (!strcmp(argv[1], "bibfs")){
        Bibfs alg("sample.txt", "output.txt");
        alg.run();
    }
    else if (!strcmp(argv[1], "sv")){
        Sv alg("sample.txt", "output.txt");
        alg.run();
    }
    else{
        cerr << "Wrong input" << endl;
    }
}