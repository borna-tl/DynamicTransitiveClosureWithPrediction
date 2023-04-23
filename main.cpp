#include <bits/stdc++.h>
#include <fstream>
#include <queue>

using namespace std;


#define MAX_NODES 1000+7
#define BFS_MODE 1
#define DFS_MODE 2
#define BIBFS_MODE 3

class reachabilityTree{//this is a simple incremental reachability tree for vertex s
public:
    reachabilityTree(){}
    reachabilityTree(int vertice_number_):vertice_number(vertice_number_){
        //should run bfs to initialize: do we consider the initialization in the timing as well?
        reachable[vertice_number] = true;
    }
    void add_edge(int u, int v, const vector<vector<int>> &out_edge){ //read on refrence vs. pointers
        if (reachable[v])
            return;
        if (!reachable[u])
            return;
        queue<int> q;
        reachable[v] = true;
        q.push(v);
        while (!q.empty()) {
            v = q.front();
            q.pop();
            for (auto i = out_edge[v].begin(); i != out_edge[v].end(); ++i){
                if (!reachable[*i]){
                    reachable[*i] = true;
                    q.push(*i);
                }
            }
        }
    }
    void print_reachability_list(){
        cout << "Reachable Nodes: ";
        for (int i = 0; i < MAX_NODES; i++)
            if (reachable[i])
                cout << i << " ";
        cout << endl;
    }
private:
    int vertice_number;
    bool reachable[MAX_NODES]; // is node u reachable from s? s--->?--->u 
    // int parent[MAX_NODES]; // (parent[v], v) is the edge 
                            // in the reachability tree whose head is v and reachable
};

class graph{
public:
    graph(string input_file_, string output_file_):
        input_file(input_file_), output_file(output_file_){
        // source_reachability[6] = new reachabilityTree(6);
        out_edge.assign(MAX_NODES, vector<int>());
        in_edge.assign(MAX_NODES, vector<int>());

    }
    void run(int mode){
        //it might be faster to read the file once, instead of keeping it open
        ifstream infile(input_file);
        bool operation;
        int u, v;
        if (mode == BFS_MODE){ //fastest way to compute variables?
            while (infile >> operation >> u >> v) {
                if (operation){
                    // cout << "operation is query on (" << u << ", " << v << ")\n";
                    bool result = run_bfs(u, v);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    // cout << "operation is edge insertion on (" << u << ", " << v << ")\n";
                    out_edge[u].push_back(v);
                    in_edge[v].push_back(u);
                    // source_reachability[6]->add_edge(u, v, out_edge); //sample
                    // source_reachability[6]->add_edge(v, u, adj); //because it's both ways
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
                    // cout << "operation is query on (" << u << ", " << v << ")\n";
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
        // source_reachability[6]->print_reachability_list();
        infile.close();
    }

private:
    string input_file, output_file;
    vector<vector<int>> out_edge; //maybe i should use vector?
    vector<vector<int>> in_edge; //maybe i should use vector?
    map<int, bool> visited_dfs;
    reachabilityTree* source_reachability[MAX_NODES];
    reachabilityTree sink_reachability[MAX_NODES]; //add functions tomorrow
    bool run_bfs(int u, int v){  
        vector<bool> visited; //maybe this should be map<int, bool>?
        visited.resize(MAX_NODES, false);
        queue<int> q;
        visited[u] = true;
        q.push(u); 
        while (!q.empty()) {
            u = q.front();
            q.pop();
            // vector<int>::iterator i;

            for (int i : out_edge[u]){
                if (!visited[i]){
                    visited[i] = true;
                    q.push(i);
                }
            }
        }
        return visited[v];
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
        vector<bool> visited_source, visited_sink; //maybe we should use map<int, bool>
        visited_source.resize(MAX_NODES, false);
        visited_sink.resize(MAX_NODES, false);
        bool found_path = false;
        queue<int> source_queue, sink_queue;
        visited_source[u] = true;
        visited_sink[v] = true;
        source_queue.push(u);
        sink_queue.push(v);
        while (!found_path && !source_queue.empty() && !sink_queue.empty()) { //i think if one queue is empty we can terminate
            //running bfs for the source queue one time
            u = source_queue.front();
            source_queue.pop();
            vector<int>::iterator i;
            for (i = out_edge[u].begin(); i != out_edge[u].end(); ++i){
                if (!visited_source[*i]){
                    visited_source[*i] = true;
                    source_queue.push(*i);
                }
                if (visited_source[*i] && visited_sink[*i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            v = sink_queue.front();
            sink_queue.pop();
            for (i = in_edge[v].begin(); i != in_edge[v].end(); ++i){
                if (!visited_sink[*i]){
                    visited_sink[*i] = true;
                    sink_queue.push(*i);
                }
                if (visited_source[*i] && visited_sink[*i]){
                    found_path = true;
                }
            }
        }
        return found_path;
    }   
};



int main(int argc, char* argv[]){
    graph G("sample.txt", "output.txt");
    G.run(!strcmp(argv[1], "dfs") ? DFS_MODE :
            !strcmp(argv[1], "bfs") ? BFS_MODE :
                !strcmp(argv[1], "bibfs") ? BIBFS_MODE : BFS_MODE);
    
}