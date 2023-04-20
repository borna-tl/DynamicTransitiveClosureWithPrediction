#include <bits/stdc++.h>
#include <fstream>

using namespace std;


#define MAX_NODES 1000+7
#define BFS_MODE 1
#define DFS_MODE 2
#define BIBFS_MODE 3

class reachabilityTree{//this is a simple incremental reachability tree for vertex s
public:
    reachabilityTree(){}
    reachabilityTree(int vertice_number_):vertice_number(vertice_number_){
        //should run bfs to initialize
        reachable[vertice_number] = true;
    }
    void add_edge(int u, int v, const map<int, list<int> > &adj){ //read on refrence vs. pointers
        if (reachable[v])
            return;
        if (!reachable[u])
            return;
  
        list<int> queue;
        reachable[v] = true;
        queue.push_back(v);
        while (!queue.empty()) {
            v = queue.front();
            queue.pop_front();
            for (auto i = adj.at(v).begin(); i != adj.at(v).end(); ++i){
                if (!reachable[*i]){
                    // cout << "u is " << u << " and adj is: " << i << endl; 
                    reachable[*i] = true;
                    // parent[*i] = v;
                    queue.push_back(*i);
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
        source_reachability[6] = new reachabilityTree(6);
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
                    adj[u].push_back(v);
                    adj[v].push_back(u);
                    source_reachability[6]->add_edge(u, v, adj); //sample
                    source_reachability[6]->add_edge(v, u, adj); //because it's both ways
                }
            }
        }
        else if (mode == DFS_MODE){
           while (infile >> operation >> u >> v) {
                if (operation){
                    // cout << "operation is query on (" << u << ", " << v << ")\n";
                    visited_dfs.clear();
                    bool result = run_dfs(u, v);
                    //maybe optimize writing to file using fwrite
                    ofstream *file = new ofstream();
                    file->open(output_file, ios::app); 
                    (*file) << "(" << u << ", " << v << ") is: " << result << endl;
                    file->close();
                }
                else{
                    // cout << "operation is edge insertion on (" << u << ", " << v << ")\n";
                    adj[u].push_back(v);
                    adj[v].push_back(u);
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
                    // cout << "operation is edge insertion on (" << u << ", " << v << ")\n";
                    adj[u].push_back(v);
                    adj[v].push_back(u);
                }
            }
        } 
        source_reachability[6]->print_reachability_list();
        infile.close();
    }

private:
    string input_file, output_file;
    // bool adj[MAX_NODES][MAX_NODES]; //this will be redundant
    // bool visited_dfs[MAX_NODES]; 
    map<int, list<int> > adj; //maybe i should use vector?
    map<int, bool> visited_dfs;
    reachabilityTree* source_reachability[MAX_NODES];
    reachabilityTree sink_reachability[MAX_NODES]; //add functions tomorrow
    bool run_bfs(int u, int v){   
        // Mark all the vertices as not visited
        vector<bool> visited; //maybe this should be map<int, bool>?
        visited.resize(MAX_NODES, false);
    
        // Create a queue for BFS
        list<int> queue;
    
        // Mark the current node as visited and enqueue it
        visited[u] = true;
        queue.push_back(u);
 
        while (!queue.empty()) {
            // Dequeue a vertex from queue
            u = queue.front();
            queue.pop_front();
    
            // Get all adjacent vertices of the dequeued
            // vertex s. If a adjacent has not been visited,
            // then mark it visited and enqueue it
            list<int>::iterator i;
            for (i = adj[u].begin(); i != adj[u].end(); ++i){
                if (!visited[*i]){
                    // cout << "u is " << u << " and adj is: " << i << endl; 
                    visited[*i] = true;
                    queue.push_back(*i);
                }
            }
        }

        return visited[v];
    }
    bool run_dfs(int u, int v){
        // Mark the current node as visited and
        visited_dfs[u] = true;
        if (visited_dfs[v])
            return true;
        // Recur for all the vertices adjacent
        // to this vertex
        list<int>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i){
            if (!visited_dfs[*i]){
                cout << "dfs going from " << u << " to " << *i << endl;
                run_dfs(*i, v);
            }
        }
        return visited_dfs[v];
    }
    bool run_bibfs(int u, int v){   
        // Mark all the vertices as not visited
        vector<bool> visited_front, visited_back; //maybe we should use map<int, bool>
        visited_front.resize(MAX_NODES, false);
        visited_back.resize(MAX_NODES, false);
        bool found_path = false;
        list<int> frontqueue, backqueue;
    
        // Mark the current node as visited and enqueue it
        visited_front[u] = true;
        visited_back[v] = true;
        frontqueue.push_back(u);
        backqueue.push_back(v);
 
        while (!found_path && !frontqueue.empty() && !backqueue.empty()) { //i think if one queue is empty we can terminate
            //running bfs for the front queue one time
            u = frontqueue.front();
            frontqueue.pop_front();
            list<int>::iterator i;
            for (i = adj[u].begin(); i != adj[u].end(); ++i){
                if (!visited_front[*i]){
                    visited_front[*i] = true;
                    frontqueue.push_back(*i);
                }
                if (visited_front[*i] && visited_back[*i]){
                    found_path = true;
                }
            }
            //running bfs for the back queue one time
            v = backqueue.front();
            backqueue.pop_front();
            for (i = adj[v].begin(); i != adj[v].end(); ++i){
                if (!visited_back[*i]){
                    visited_back[*i] = true;
                    backqueue.push_back(*i);
                }
                if (visited_front[*i] && visited_back[*i]){
                    found_path = true;
                }
            }
        }

        return found_path;
    }
        
};



int main(){
    graph G("sample.txt", "output.txt");
    G.run(BFS_MODE);
    
}