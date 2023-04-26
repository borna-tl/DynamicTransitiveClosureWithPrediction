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

class graph{
public:
    graph(string input_file_, string output_file_):
        input_file(input_file_), output_file(output_file_){
        out_edge.assign(MAX_NODES, vector<int>());
        in_edge.assign(MAX_NODES, vector<int>());
        visited_bfs.resize(MAX_NODES, false);
        visited_bibfs_source.resize(MAX_NODES, false);
        visited_bibfs_sink.resize(MAX_NODES, false);
    }
    void run(int mode){
        //it might be faster to read the file once, instead of keeping it open
        ifstream infile(input_file);
        bool operation;
        int u, v;
        if (mode == BFS_MODE){ //fastest way to compute variables?
            while (infile >> operation >> u >> v) {
                if (operation){
                    bool result = run_bfs(u, v);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    out_edge[u].push_back(v);
                    in_edge[v].push_back(u);
                }
            }
        }
        else if (mode == DFS_MODE){
           while (infile >> operation >> u >> v) {
                if (operation){
                    visited_dfs.clear();
                    bool result = run_dfs(u, v);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    out_edge[u].push_back(v);
                    in_edge[v].push_back(u); 
                }
            }
        } 
        else if (mode == BIBFS_MODE){ //note, maybe use multi-directional BFS? can we figure out what that is?
            while (infile >> operation >> u >> v) {
                if (operation){
                    bool result = run_bibfs(u, v);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    out_edge[u].push_back(v);
                    in_edge[v].push_back(u);                     
                }
            }
        }
        if (mode == SV_MODE){ //fastest way to compute variables?
            vector <int> sv_list = generate_sv_list();
            while (infile >> operation >> u >> v) {
                if (operation){
                    bool result = run_sv(u, v, sv_list);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    out_edge[u].push_back(v);
                    in_edge[v].push_back(u);
                    update_sv(u, v, sv_list);
                }
            }
        } 
        infile.close();
    }

private:
    string input_file, output_file;
    vector<vector<int>> out_edge;
    vector<vector<int>> in_edge;
    map<int, bool> visited_dfs;
    vector<bool> visited_bfs; //in order to reduce resizing cost
    vector<bool> visited_bibfs_source, visited_bibfs_sink;
    reachabilityTree* reachability_tree[MAX_NODES];
    bool run_bfs(int u, int v){  
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
    bool run_dfs(int u, int v){
        visited_dfs[u] = true;
        if (visited_dfs[v])
            return true;
        vector<int>::iterator i;
        for (i = out_edge[u].begin(); i != out_edge[u].end(); ++i){
            if (!visited_dfs[*i]){
                run_dfs(*i, v);
            }
        }
        return visited_dfs[v];
    }
    bool run_bibfs(int u, int v){   
        
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
            vector<int>::iterator i;
            for (i = out_edge[u].begin(); i != out_edge[u].end(); ++i){
                if (!visited_bibfs_source[*i]){
                    visited_bibfs_source[*i] = true;
                    source_queue.push_back(*i);
                }
                if (visited_bibfs_source[*i] && visited_bibfs_sink[*i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            v = sink_queue[sink_pointer];
            sink_pointer++;
            for (i = in_edge[v].begin(); i != in_edge[v].end(); ++i){
                if (!visited_bibfs_sink[*i]){
                    visited_bibfs_sink[*i] = true;
                    sink_queue.push_back(*i);
                }
                if (visited_bibfs_source[*i] && visited_bibfs_sink[*i]){
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
    bool run_sv(int u, int v, const vector<int>& sv_list){  
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
        return run_bfs(u, v);
    }
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
    graph G("sample.txt", "output.txt");
    G.run(!strcmp(argv[1], "dfs") ? DFS_MODE :
            !strcmp(argv[1], "bfs") ? BFS_MODE :
                !strcmp(argv[1], "bibfs") ? BIBFS_MODE :
                    !strcmp(argv[1], "sv") ? SV_MODE : BFS_MODE);

    
}