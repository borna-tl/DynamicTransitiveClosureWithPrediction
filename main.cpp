#include <bits/stdc++.h>
#include <fstream>
#include <queue>

using namespace std;


#define MAX_NODES 1000+7
#define BFS_MODE 1
#define DFS_MODE 2
#define BIBFS_MODE 3
#define SV_MODE 4

//class algorithm with child that handle insertion and queries (like bfs, dfs, etc)

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
            for (auto i = edge[v].begin(); i != edge[v].end(); ++i){
                if (!r[*i]){
                    r[*i] = true;
                    q.push(*i);
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
            cout << "cinstructing algoihms" << endl;
        out_edge.assign(MAX_NODES, vector<int>());
        in_edge.assign(MAX_NODES, vector<int>());
        
    }
    virtual void run(int mode){
        //it might be faster to read the file once, instead of keeping it open
        ifstream infile(input_file);
        bool operation;
        int u, v;
        if (mode == BFS_MODE){ //fastest way to compute variables?
            
        }
        else if (mode == DFS_MODE){
           
        } 
        else if (mode == BIBFS_MODE){ //note, maybe use multi-directional BFS? can we figure out what that is?
            
        }
        if (mode == SV_MODE){ //fastest way to compute variables?
            
        } 
        infile.close();
    }

protected:
    string input_file, output_file;
    vector<vector<int>> out_edge;
    vector<vector<int>> in_edge;
    
    
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
    void run(){
        ifstream infile(input_file);
        ofstream *file = new ofstream();
        file->open(output_file, ios::ate); 
        bool operation;
        int u, v;
        while (infile >> operation >> u >> v) {
            if (operation){
                bool result = calculate_bfs(u, v);                
                (*file) << "(" << u << ", " << v << ") is: " << result << endl;
            }
            else{
                out_edge[u].push_back(v);
                in_edge[v].push_back(u);
            }
        }
        file->close();
    }
private:
    vector<bool> visited_bfs; //in order to reduce resizing cost
};

class Dfs : public Algorithms{

public:
    Dfs(){"DFS Default constructor called! Check something...\n\n";}
    Dfs(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        visited_dfs.resize(MAX_NODES, false);
        round_number = false;
    }
    bool calculate_dfs(int u, int v){  //maybe we should try non-iterative dfs as well
        visited_dfs[u] = true;
        if (visited_dfs[v]){
            return true;
        }
        for (auto i : out_edge[u]){
            if (!visited_dfs[i]){
                calculate_dfs(i, v);
            }
        }
        return visited_dfs[v];
    }
    void run(){
        ifstream infile(input_file);
        ofstream *file = new ofstream();
        file->open(output_file, ios::ate); 
        bool operation;
        int u, v;
        while (infile >> operation >> u >> v) {
            if (operation){
                visited_dfs.clear();
                visited_dfs.resize(MAX_NODES, false); //maybe not so efficient
                bool result = calculate_dfs(u, v);
                (*file) << "(" << u << ", " << v << ") is: " << result << endl;
            }
            else{
                out_edge[u].push_back(v);
                in_edge[v].push_back(u); 
            }
        }
        file->close();
    }
private:
    vector<bool> visited_dfs;
    bool round_number;
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
    void run(){
        ifstream infile(input_file);
        ofstream *file = new ofstream();
        file->open(output_file, ios::ate); 
        bool operation;
        int u, v;
        while (infile >> operation >> u >> v) {
            if (operation){
                bool result = calculate_bibfs(u, v);                
                (*file) << "(" << u << ", " << v << ") is: " << result << endl;
            }
            else{
                out_edge[u].push_back(v);
                in_edge[v].push_back(u);
            }
        }
        file->close();
    }
private:
    vector<bool> visited_bibfs_source;
    vector<bool> visited_bibfs_sink;

};

class Sv : public Algorithms{

public:
    Sv(){"Sv Default constructor called! Check something...\n";}
    Sv(string input_file_, string output_file_) : Algorithms(input_file_, output_file_){
        
    }
    bool calculate_sv(int u, int v, const vector<int>& sv_list){  
        cout << "For insertion on (" << u << ", " << v << "): " << endl;
        reachability_tree[6]->print_reachability_list(); //for testing

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
        cout << "sv not successful, falling back..." << endl;
        return Bfs(input_file, output_file).calculate_bfs(u, v);
        //get back later
    
    }
    void run(){
        ifstream infile(input_file);
        ofstream *file = new ofstream();
        file->open(output_file, ios::ate); 
        bool operation;
        int u, v;
        vector <int> sv_list = generate_sv_list();
        while (infile >> operation >> u >> v) {
            if (operation){
                bool result = calculate_sv(u, v, sv_list);                
                (*file) << "(" << u << ", " << v << ") is: " << result << endl;
            }
            else{
                out_edge[u].push_back(v);
                in_edge[v].push_back(u);
                update_sv(u, v, sv_list);
            }
        }
        file->close();

    }
private:
    reachabilityTree* reachability_tree[MAX_NODES];
    vector<int> generate_sv_list(){
        vector<int> sv_list;
        reachability_tree[6] = new reachabilityTree(6); //using node 6 as SV (will add rand later)
        sv_list.push_back(6);
        return sv_list;
    }
    void update_sv(int u, int v, const vector<int>& sv_list){
        for (auto sv : sv_list){
            reachability_tree[sv]->update(u, v, out_edge, in_edge);
        }
    }
};



int main(int argc, char* argv[]){
    Algorithms G("sample.txt", "output.txt");
    Sv b("sample.txt", "output.txt");
    b.run();
    
    // G.run(!strcmp(argv[1], "dfs") ? DFS_MODE :
    //         !strcmp(argv[1], "bfs") ? BFS_MODE :
    //             !strcmp(argv[1], "bibfs") ? BIBFS_MODE :
    //                 !strcmp(argv[1], "sv") ? SV_MODE : BFS_MODE);

    
}