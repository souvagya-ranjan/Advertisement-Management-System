#include <fstream>
#include <iostream>
#include <bits/stdc++.h>

#include <mpi.h>
#include <omp.h>
using namespace std;


void DFS(unordered_map<int, set<int>>& graph, int node, unordered_set<int>& visited, unordered_set<int>& component) {
    visited.insert(node);
    component.insert(node);
    for (auto it=graph[node].begin(); it!=graph[node].end(); it++) {
        int neighbor = *it;
        if (visited.find(neighbor) == visited.end()) {
            DFS(graph, neighbor, visited, component);
        }
    }
}

vector<unordered_set<int>> findConnectedComponents(unordered_map<int, set<int>>& graph) {
    vector<unordered_set<int>> connectedComponents;
    unordered_set<int> visited;
    for (auto it = graph.begin(); it!=graph.end();it++) {
        int node = it->first;
        if (visited.find(node) == visited.end()) {
            unordered_set<int> component;
            DFS(graph, node, visited, component);
            connectedComponents.push_back(component);
        }
    }
    return connectedComponents;
}


int main(int argc, char* argv[]){
    string a1 = argv[1];
    string a2 = argv[2];
    string a3 = argv[3];
    string a4 = argv[4];
    string a5 = argv[5];
    string a6 = argv[6];
    string a7 = argv[7];
    string a8 = argv[8];
    int taskid = stoi(a1.substr(a1.find("=")+1));
    string inputpath = a2.substr(a2.find("=")+1);
    string headerpath = a3.substr(a3.find("=")+1);
    string outputpath = a4.substr(a4.find("=")+1);
    int verbose = stoi(a5.substr(a5.find("=")+1));
    int startk = stoi(a6.substr(a6.find("=")+1));
    int endk = stoi(a7.substr(a7.find("=")+1));
    int p  = stoi(a8.substr(a8.find("=")+1));
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //open the inputfile
    ifstream inputfile(inputpath);
    ifstream headerfile(headerpath);

    int n = 0, m = 0;
    inputfile.read((char*)&n, 4);
    inputfile.read((char*)&m, 4);

    //divide the vertices among the processes
    int vertices_per_process = n/comm_size;
    int start_n = rank*vertices_per_process;
    int end_n = (rank+1)*vertices_per_process;
    if(rank == comm_size-1){
        end_n = n;
    }

    //read the subgraph for each file 
    unordered_map <int, set<int>> subgraph;
    for(int i = 0; i < n; i++){
        int vertex = 0;
        inputfile.read((char*)&vertex, 4);
        if(vertex >= start_n && vertex < end_n){
            int degree = 0;
            inputfile.read((char*)&degree, 4);
            set <int> neighbours;
            for(int j = 0; j < degree; j++){
                int neighbour = 0;
                inputfile.read((char*)&neighbour, 4);
                neighbours.insert(neighbour);
            }
            subgraph[vertex] = neighbours;
        }
        else{
            int degree = 0;
            inputfile.read((char*)&degree, 4);
            for(int j = 0; j < degree; j++){
                int neighbour = 0;
                inputfile.read((char*)&neighbour, 4);
            }
        }
    }

    //triangle enumeration using openmp
    map< pair <int, int>, set<int>>triangle_map;
    for(auto it = subgraph.begin(); it != subgraph.end(); it++){
        int vertex = it->first;
        set <int> neighbours = it->second;
        for(auto it1 = neighbours.begin(); it1 != neighbours.end(); it1++){
            int vertex_v = *it1;
            if(vertex > vertex_v){
                continue;
            }
            triangle_map[make_pair(vertex, vertex_v)] = set<int>();
            for(auto it2 = neighbours.begin(); it2 != neighbours.end(); it2++){
                int vertex_w = *it2;
                if(vertex_v != vertex_w){
                    if(vertex_v>=start_n && vertex_v < end_n){
                        if(subgraph[vertex_v].find(vertex_w) != subgraph[vertex_v].end()){
                        //triangle found
                            triangle_map[make_pair(vertex, vertex_v)].insert(vertex_w);
                        }
                    }
                    else if(vertex_w >= start_n && vertex_w < end_n){
                        if(subgraph[vertex_w].find(vertex_v) != subgraph[vertex_w].end()){
                            //triangle found
                            triangle_map[make_pair(vertex, vertex_v)].insert(vertex_w);
                        }
                    }
                    else{
                        headerfile.seekg(4*vertex_v, ios::beg);
                        int offset = 0;
                        headerfile.read((char*)&offset, 4);
                        inputfile.seekg(offset, ios::beg);
                        int node_id = 0;
                        inputfile.read((char*)&node_id, 4);
                        int degree = 0;
                        inputfile.read((char*)&degree, 4);
                        set<int> neighbours_v;
                        for(int i = 0; i < degree; i++){
                            int neighbour = 0;
                            inputfile.read((char*)&neighbour, 4);
                            neighbours_v.insert(neighbour);
                        }
                        if (neighbours_v.find(vertex_w) != neighbours_v.end()){
                            //triangle found
                            triangle_map[make_pair(vertex, vertex_v)].insert(vertex_w);
                        }
                    }
                }
            }
        }
    }
    inputfile.close();
    headerfile.close();
    if(rank != 0){
        int size_of_triangle_map = triangle_map.size();
        int total_map_size = 3*size_of_triangle_map; //initalised the total size 

        for(auto it = triangle_map.begin(); it != triangle_map.end(); it++){
            total_map_size += it->second.size();
        }
        auto triangle_map_array = new int[total_map_size];
        int index = 0;
        for(auto it = triangle_map.begin(); it != triangle_map.end(); it++){
            triangle_map_array[index++] = it->first.first;
            triangle_map_array[index++] = it->first.second;
            triangle_map_array[index++] = it->second.size();
            for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
                triangle_map_array[index++] = *it1;
            }
        }
        MPI_Send(&total_map_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(triangle_map_array, total_map_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    if(rank == 0){
        for(int rank_no = 1; rank_no < comm_size; rank_no++){
            int size_of_recv_triangle_map = 0;
            MPI_Recv(&size_of_recv_triangle_map, 1, MPI_INT, rank_no, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            auto triangle_map_array = new int[size_of_recv_triangle_map];
            MPI_Recv(triangle_map_array, size_of_recv_triangle_map, MPI_INT, rank_no, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int index = 0;
            while(index < size_of_recv_triangle_map){
                int vertex = triangle_map_array[index++];
                int neighbour = triangle_map_array[index++];
                int size = triangle_map_array[index++];
                set <int> neighbours;
                for(int j = 0; j < size; j++){
                    neighbours.insert(triangle_map_array[index++]);
                }
                triangle_map[make_pair(vertex, neighbour)] = neighbours;
            }
        }

        //create support map 
        map< pair <int, int>, int>support_map;
        for(auto it = triangle_map.begin(); it != triangle_map.end(); it++){
            support_map[it->first] = it->second.size();
        }
        cout<<support_map.size()<<endl;
        //read graph for k-trruss computation 
        unordered_map < int, set<int>>graph;
        ifstream inputfile(inputpath, ios::binary);
        int nn, mm;
        inputfile.read((char*)&nn, 4);
        inputfile.read((char*)&mm, 4);
        for(int i = 0; i < n; i++){
            int vertex = 0;
            inputfile.read((char*)&vertex, 4);
            int degree = 0;
            inputfile.read((char*)&degree, 4);
            set <int> neighbours;
            for(int j = 0; j < degree; j++){
                int neighbour = 0;
                inputfile.read((char*)&neighbour, 4);
                neighbours.insert(neighbour);
            }
            graph[vertex] = neighbours;
        }
        inputfile.close();

        //k-truss computation
        ofstream outputfile(outputpath);
        if(taskid == 1){
            for(int k = startk; k <= endk; k++){
                queue<pair<int, int>> deletable;
                for(auto it = support_map.begin(); it != support_map.end(); it++){
                    if(it->second < k){
                        deletable.push(it->first);
                    }
                }
                while(!deletable.empty()){
                    pair<int, int> edge = deletable.front();
                    deletable.pop();
                    int u = edge.first;
                    int v = edge.second;
                    if(graph[u].find(v)==graph[u].end()){
                        continue;
                    }
                    graph[u].erase(v);
                    graph[v].erase(u);

                    for(auto it1 = triangle_map[edge].begin(); it1 != triangle_map[edge].end(); it1++){
                        int w = *it1;
                        if(graph[u].find(w) != graph[u].end()){
                            pair<int, int> edge_uw = make_pair(min(u, w), max(u, w));
                            support_map[edge_uw]--;
                            triangle_map[edge_uw].erase(v);
                            if(support_map[edge_uw] < k){
                                deletable.push(edge_uw);
                            }
                        }
                        if(graph[v].find(w) != graph[v].end()){
                            pair<int, int> edge_vw = make_pair(min(v, w), max(v, w));
                            support_map[edge_vw]--;
                            triangle_map[edge_vw].erase(u);
                            if(support_map[edge_vw] < k){
                                deletable.push(edge_vw);
                            }   
                        }
                    }

                }

                vector< pair< int, int>> edge_remove;
                for(auto it = support_map.begin(); it != support_map.end(); it++){
                    if(it->second < k){
                        edge_remove.push_back(it->first);
                    }
                }
                for(auto it = edge_remove.begin(); it != edge_remove.end(); it++){
                    pair<int, int> edge = *it;
                    int u = edge.first;
                    int v = edge.second;
                    if(graph[u].find(v) != graph[u].end()){
                        graph[u].erase(v);
                        graph[v].erase(u);
                    }
                }
                vector<int> verticesToRemove;
                for (auto it = graph.begin(); it != graph.end(); it++) {
                    if (it->second.empty()) {
                        verticesToRemove.push_back(it->first);
                    }
                }
                for (auto v : verticesToRemove) {
                    graph.erase(v);
                }

                //writing part 
                string s = "0\n";
                if(verbose == 0){
                    if(graph.size()!=0){
                        s = "1\n";
                    }
                }
                else{
                    if(graph.size()!=0){
                        vector< unordered_set<int>> kTruss = findConnectedComponents(graph);
                        s = "1\n";
                        s+= to_string(kTruss.size()) + "\n";
                        for(auto it = kTruss.begin(); it != kTruss.end(); it++){
                            for(auto it1 = it->begin(); it1 != it->end(); it1++){
                                s += to_string(*it1) + " ";
                            }
                            s += "\n";
                        }
                    }
                }
                outputfile <<s;
            }

        }
        else{
            int k = endk;
            queue<pair<int, int>> deletable;
            for(auto it = support_map.begin(); it != support_map.end(); it++){
                if(it->second < k){
                    deletable.push(it->first);
                }
            }
            while(!deletable.empty()){
                pair<int, int> edge = deletable.front();
                deletable.pop();
                int u = edge.first;
                int v = edge.second;
                if(graph[u].find(v)==graph[u].end()){
                    continue;
                }
                graph[u].erase(v);
                graph[v].erase(u);

                for(auto it1 = triangle_map[edge].begin(); it1 != triangle_map[edge].end(); it1++){
                    int w = *it1;
                    if(graph[u].find(w) != graph[u].end()){
                        pair<int, int> edge_uw = make_pair(min(u, w), max(u, w));
                        support_map[edge_uw]--;
                        triangle_map[edge_uw].erase(v);
                        if(support_map[edge_uw] < k){
                            deletable.push(edge_uw);
                        }
                    }
                    if(graph[v].find(w) != graph[v].end()){
                        pair<int, int> edge_vw = make_pair(min(v, w), max(v, w));
                        support_map[edge_vw]--;
                        triangle_map[edge_vw].erase(u);
                        if(support_map[edge_vw] < k){
                            deletable.push(edge_vw);
                        }   
                    }
                }

            }

            vector< pair< int, int>> edge_remove;
            for(auto it = support_map.begin(); it != support_map.end(); it++){
                if(it->second < k){
                    edge_remove.push_back(it->first);
                }
            }
            for(auto it = edge_remove.begin(); it != edge_remove.end(); it++){
                pair<int, int> edge = *it;
                int u = edge.first;
                int v = edge.second;
                if(graph[u].find(v) != graph[u].end()){
                    graph[u].erase(v);
                    graph[v].erase(u);
                }
            }
            vector<int> verticesToRemove;
            for (auto it = graph.begin(); it != graph.end(); it++) {
                if (it->second.empty()) {
                    verticesToRemove.push_back(it->first);
                }
            }
            for (auto v : verticesToRemove) {
                graph.erase(v);
            }

            vector<unordered_set<int>> KTruss = findConnectedComponents(graph);
            unordered_map<int,int> trussMap;
            for(int idx = 0; idx< KTruss.size(); idx++){
                for(auto it1 = KTruss[idx].begin(); it1 != KTruss[idx].end(); it1++){
                    trussMap[*it1] = idx;
                }
            }

            //read new influencer graph
            unordered_map<int, set<int>> influencerGraph;
            ifstream infile;
            infile.open(inputpath, ios::in | ios::binary);
            int nnn = 0,mmm = 0;
            infile.read((char*)&nnn, 4);
            infile.read((char*)&mmm, 4);
            for(int i = 0; i < n; i++){
                int vertex = 0;
                infile.read((char*)&vertex, 4);
                int degree = 0;
                infile.read((char*)&degree, 4);
                set<int> influencedGroups;
                for(int j = 0; j < degree; j++){
                    int neighbour = 0;
                    infile.read((char*)&neighbour, 4);
                    if(trussMap.find(neighbour)!= trussMap.end()){
                        influencedGroups.insert(trussMap[neighbour]);
                    }
                }
                if(influencedGroups.size()>=p){
                    influencerGraph[vertex] = influencedGroups;
                }
            }
            infile.close();

            string s;
            if(verbose == 0)
            {
                s = "0 \n";
                if(influencerGraph.size() > 0){
                    s = to_string(influencerGraph.size()) + "\n";
                    for(auto it = influencerGraph.begin(); it != influencerGraph.end(); it++){
                        s += to_string(it->first) + " ";
                    }
                    s += "\n";
                }
            }
            else
            {
                s = "0 \n";
                if(influencerGraph.size() > 0){
                    s = to_string(influencerGraph.size()) + "\n";
                    for(auto it = influencerGraph.begin(); it != influencerGraph.end(); it++){
                        s += to_string(it->first) + "\n";
                        for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                            int component_no = *it2;
                            for(auto it3 = KTruss[component_no].begin(); it3!= KTruss[component_no].end(); it3++){
                                s += to_string(*it3) + " ";
                            }
                        }
                        s += "\n";
                    }
                }
            }
            outputfile << s;
        }
        outputfile.close();
    }
    
    MPI_Finalize();
    return 0;

}