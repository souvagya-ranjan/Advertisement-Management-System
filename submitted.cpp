#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <mpi.h>
#include "omp.h"

using namespace std;

bool cmp_for_vertexSort(pair <int, int> a, pair <int, int> b){ 
    if (a.second == b.second){
        return a.first < b.first;
    }
    return a.second < b.second;
}

void printGraph(unordered_map<int, set<int>>& graph){
    stringstream ss;
    cout<<graph.size()<<endl;
    // iterate over each vertex in the graph
    for (auto& v : graph) {
        // add the vertex id to the stringstream
        ss << v.first << " -> ";

        // add the adjacent vertices to the stringstream
        for (auto& adj_v : v.second) {
            ss << adj_v << " ";
        }

        // add a newline character to separate each vertex
        ss << endl;
    }

    // print the graph string
    cout << ss.str() << endl;
}

void readGraph(int n, int m, unordered_map<int,set<int>>& graph, string iInputPath, string iHeaderPath){
    ifstream input_file(iInputPath, ios::binary);
    if (!input_file.is_open()) {
        cout << "Error opening input file"<<endl;
        return;
    }
    ifstream header_file(iHeaderPath, ios::binary);
    if (!header_file.is_open()) {
        cout << "Error opening header file"<<endl;
        return;
    }

    int a = 0;
    input_file.read((char*) &a, 4);
    n = a;
    input_file.read((char*) &a, 4);
    m = a;

    // cout<<"n: "<< n << " "<< "m:"<<m<<endl;
    
    vector <int> offset(n);
    for (int i = 0; i < n; i++){
        int offset_value = 0;
        header_file.read((char*) &offset_value, 4);
        offset[i] = offset_value;
    }
    
    for(int i = 0; i < n; i++){
        input_file.seekg(offset[i], ios::beg);
        
        int degree = 0, node_id = 0;
        
        input_file.read((char*) &node_id, 4);
        input_file.read((char*) &degree, 4);
        
        set <int> neighbour_list;
        for (int i = 0; i < degree; i++){
            int neighbour = 0;
            input_file.read((char*) &neighbour, 4);
            neighbour_list.insert(neighbour);
        }
        graph[node_id] = neighbour_list;
    }
    input_file.close();
    header_file.close();
}

void writeGraph(vector<string> toPrint, string iOutputPath, int iVerbose){
    ofstream outfile(iOutputPath);
    if(iVerbose == 0){
        for (int i=0; i<toPrint.size(); i++){
            if (i!=toPrint.size()-1){
                outfile << toPrint[i] << " ";
            }
            else{
                outfile << toPrint[i];
            }
        }
    }
    else{
        for(int i = 0; i < toPrint.size(); i++){
            outfile<<toPrint[i];
        }
    }
    outfile.close();  
}

void convertSupportMap(map <pair<int, int>, int>& support, int* arrayOfSupport){
    int i = 0;
    for (auto it = support.begin(); it != support.end(); it++)
    {
        int u = it->first.first;
        int v = it->first.second;
        int w = it->second;
        
        arrayOfSupport[i] = u; arrayOfSupport[i+1] = v; arrayOfSupport[i+2] = w;
        i+=3;
    }
}

void convertTriangleMap(map<pair<int,int>,set<int>>& triangle, int* triangleMapArray, int totalSizeOfTriangleMap){
    int i = 0;
    for (auto it = triangle.begin(); it != triangle.end(); it++)
    {
        int u = it->first.first;
        int v = it->first.second;
        set<int> w = it->second;
        
        triangleMapArray[i] = u; triangleMapArray[i+1] = v;
        i+=2;
        for (auto it = w.begin(); it != w.end(); it++)
        {
            triangleMapArray[i] = *it;
            i++;
        }
        triangleMapArray[i] = -1;
        i++;
    }
}

void Initialize(int k, vector<pair<int,int>>& deletable, map<pair<int,int>,int>& support_map) {
    for (auto it = support_map.begin(); it != support_map.end(); it++)
    {
        pair<int, int> edge = it->first;
        int w = it->second;
        if (w < k-2)
        {
            deletable.push_back(edge);
        }
    }    
}

void FilterEdges(unordered_map<int,set<int>>& graph, int k, vector<pair<int,int>>& deletable, map<pair<int,int>,int>& support_map, map<pair<int,int>,set<int>>& triangle_map, int num_t) {

    while (1)
    {
        if (deletable.empty())
        {
            break;
        }

        int currSize = deletable.size();
        vector<pair<int,int>> newDeletable;

        #pragma omp parallel for num_threads(num_t)
        for (int i = 0; i < currSize; i++)
        {
            pair<int,int> uv = deletable[i];
            int u = uv.first, v = uv.second;
            // deletable.pop();
            if (graph[u].find(v)==graph[u].end() || graph[v].find(u)==graph[v].end())
            {
                continue;
            }
            #pragma omp critical
            {
                graph[u].erase(v);
                graph[v].erase(u);
            }
            // graph[u].erase(v);
            // graph[v].erase(u);

            set<int> vertices_UV = triangle_map[uv];

            for (auto it=vertices_UV.begin(); it!=vertices_UV.end(); it++)
            {
                int w = *it;
                pair<int,int> uw = make_pair(u,w);
                pair<int,int> vw = make_pair(v,w);

                pair<int,int> wu = make_pair(w,u);
                pair<int,int> wv = make_pair(w,v);

                #pragma omp critical
                {
                    if ((graph[u].find(w)!=graph[u].end() && graph[w].find(u)!=graph[w].end()))
                    {
                        if (support_map.find(uw) != support_map.end())
                        {   
                            support_map[uw]--;
                            triangle_map[uw].erase(v);
                            if (support_map[uw] < k-2)
                            {
                                newDeletable.push_back(uw);
                            }
                        }
                        else if (support_map.find(wu) != support_map.end())
                        {
                            support_map[wu]--;
                            triangle_map[wu].erase(v);
                            if (support_map[wu] < k-2)
                            {
                                newDeletable.push_back(wu);
                            }
                        
                        }
                    }

                    if ((graph[v].find(w)!=graph[v].end() && graph[w].find(v)!=graph[w].end()))
                    {
                        if (support_map.find(vw) != support_map.end())
                        {
                            support_map[vw]--;
                            triangle_map[vw].erase(u);
                            if (support_map[vw] < k-2)
                            {
                                newDeletable.push_back(vw);
                            }
                        }
                        else if (support_map.find(wv) != support_map.end())
                        {
                            support_map[wv]--;
                            triangle_map[wv].erase(u);
                            if (support_map[wv] < k-2)
                            {
                                newDeletable.push_back(wv);
                            }
                        }
                    }
                }                  
            }  
        }

        #pragma omp crictical
        {
            deletable = newDeletable; 
        }
    }
    
}

// , int k, map<pair<int,int>,int>& support_map, map<pair<int,int>,set<int>>& triangle_map
void modifyGraph(unordered_map<int,set<int>>& graph) 
{

    vector<int> verticesToRemove;
    for (auto it = graph.begin(); it != graph.end(); it++) {
        if (it->second.empty()) {
            verticesToRemove.push_back(it->first);
        }
    }
    for (auto v : verticesToRemove) {
        graph.erase(v);
    }
}

void findKTruss(unordered_map < int, set<int>>& graph, int k){   
    queue<pair<int,int>> edgeToDelete;
    
    for (auto it = graph.begin(); it!=graph.end(); it++) {
        int u = it->first;
        for (auto it = graph[u].begin(); it!=graph[u].end(); it++) {
            int v = *it;
            if (u < v) { 
                set<int> vertices_UV;
                set_intersection(graph[u].begin(), graph[u].end(), graph[v].begin(), graph[v].end(), inserter(vertices_UV, vertices_UV.begin()));         
                int nsupport = vertices_UV.size();

                if (nsupport < k-2) {
                    edgeToDelete.push(make_pair(u,v));
                }
            }
        }
    }

    while (!edgeToDelete.empty()) {
        pair<int,int> edge = edgeToDelete.front();
        edgeToDelete.pop();

        int u = edge.first, v = edge.second;
        set<int> vertices_UV;
        set_intersection(graph[u].begin(), graph[u].end(), graph[v].begin(), graph[v].end(), inserter(vertices_UV, vertices_UV.begin()));         

        if (vertices_UV.size() < k-2) {
            graph[u].erase(v);
            graph[v].erase(u);
        }    
    }
}


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

void convertInfluencerArrayToMap(int* influencerArray, int size, unordered_map<int, set<int>>& influencerMap) {
    int i = 0;
    while (i < size) {
        int node = influencerArray[i];
        influencerMap[node] = set<int>();
        int count = influencerArray[i+1];
        for(int j = 0; j < count; j++) {
            influencerMap[node].insert(influencerArray[i+2+j]);
        }
        i = i + 2 + count;
    }
}

int main(int argc, char* argv[])
{
    
    string taskid = argv[1];
    string inputpath = argv[2];
    string headerpath = argv[3];
    string outputpath = argv[4];
    string verbose = argv[5];
    string startK = argv[6];
    string endK = argv[7];
    string p = argv[8];

    int iTaskId        = stoi(taskid.substr(taskid.find('=') + 1));
    string iInputPath  = inputpath.substr(inputpath.find('=') + 1);
    string iHeaderPath = headerpath.substr(headerpath.find('=') + 1);
    string iOutputPath = outputpath.substr(outputpath.find('=') + 1);
    int iVerbose       = stoi(verbose.substr(verbose.find('=') + 1));
    int istartK        = stoi(startK.substr(startK.find('=')+1));
    int iendK          = stoi(endK.substr(endK.find('=')+1));
    int iP             = stoi(p.substr(p.find('=')+1));

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype tuple_type_4;
    MPI_Type_contiguous(4, MPI_INT, &tuple_type_4);
    MPI_Type_commit(&tuple_type_4);

    ifstream infile(iInputPath, ios::binary);
    ifstream headerfile(iHeaderPath, ios::binary);

    vector <pair<int, int>> vertex_degree_list;  // first value of each element in the vector is node and second is degree
    unordered_map <int, int > process_map; // key = node_id and value  = process_id that will process it 
    
    //this part is processes by all the processes
    int n = 0, m = 0;
    infile.read((char*)&n, 4);
    infile.read((char*)&m, 4);
    
    for(int i = 0; i < n; i++){
        int offset = 0;
        //vertex = i;
        headerfile.seekg(4*i, ios::beg);
        headerfile.read((char*)&offset, 4);
        
        int node_id = 0;
        int degree = 0;
        
        infile.seekg(offset, ios::beg);
        infile.read((char*)&node_id, 4);
        infile.read((char*)&degree, 4);
        
        vertex_degree_list.push_back(make_pair(node_id, degree));
    }
    // now it will sort the vertices based on degree
    sort(vertex_degree_list.begin(), vertex_degree_list.end(), cmp_for_vertexSort);
        
    // now it will distribute the vertices to the processes
    int process = 0;
    for(int i = 0; i < vertex_degree_list.size(); i++){
        process_map[vertex_degree_list[i].first] = process;
        process++;
        if(process == size){
            process = 0;
        }
    }    

    //this part is for each individual process 
    unordered_map <int , set<int>> process_graph;
    map <pair<int, int>, int> support_map;
    map < pair<int, int> , set<int> > triangle_map;
    for (int i = 0; i < n; i++){
        if (process_map[i] == rank){
            int offset = 0;
            headerfile.seekg(4*i, ios::beg);
            headerfile.read((char*)&offset, 4);
            
            int node_id = 0;
            int degree = 0;
            
            infile.seekg(offset, ios::beg);
            infile.read((char*)&node_id, 4);
            infile.read((char*)&degree, 4);
            
            set<int> neighbour_list;
            for(int j = 0; j < degree; j++){
                int neighbour = 0;
                infile.read((char*)&neighbour, 4);
                neighbour_list.insert(neighbour);
            }
            process_graph[node_id] = neighbour_list;
        }
    }
    //triangle enumeration 
    //first store a list of all the edges owned by the process 
    int process_m = 0;
    for(auto it = process_graph.begin(); it != process_graph.end(); it++){
        int u = it->first;
        int degree_u = it->second.size();
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
            int v = *it2;
            //find the owner of edge <u,v>
            int offset = 0;
            headerfile.seekg(4*v, ios::beg);
            headerfile.read((char*)&offset, 4);

            int degree_v = 0;
            infile.seekg(offset+4, ios::beg);
            infile.read((char*)&degree_v, 4);

            //check if the edge belongs to the current process
            bool flag = cmp_for_vertexSort(make_pair(u, degree_u), make_pair(v, degree_v));
            if (!flag){
                continue;
            }
            process_m ++;
            support_map[make_pair(u,v)] = 0;
            triangle_map[make_pair(u,v)] = set<int>();
        }
    }
    vector<int> send_edge_counts(size, process_m);
    vector<int> recv_edge_counts(size, 0);
    MPI_Alltoall(&send_edge_counts[0], 1, MPI_INT, &recv_edge_counts[0], 1, MPI_INT, MPI_COMM_WORLD);
    int max_recv_edge_counts = 0;
    for(int i = 0; i < size; i++){
        // if(rank == 0)
            // cout<<recv_edge_counts[i]<< " ";
        if (recv_edge_counts[i] > max_recv_edge_counts){
            max_recv_edge_counts = recv_edge_counts[i];
        }
    }
    // if (rank == 0)
    //     cout<<endl;
    // if (rank == 0)
    //     cout<<"max_edges: "<<max_recv_edge_counts<<endl;
    int edge_count = 0;
    while(edge_count< max_recv_edge_counts){
        if(edge_count >= process_m){
            edge_count++;
            vector <int> ReqRecvCounts(size, 0);
            vector <int> ReqSendCounts(size, 0);
            MPI_Alltoall(&ReqSendCounts[0], 1, MPI_INT, &ReqRecvCounts[0], 1, MPI_INT, MPI_COMM_WORLD);

            vector <int> ReqSendDispls(size, 0);
            vector <int> ReqRecvDispls(size, 0);

            for(int i = 0; i<size; i++){
                if(i!=0){
                    ReqRecvDispls[i] = ReqRecvDispls[i-1] + ReqRecvCounts[i-1];
                    ReqSendDispls[i] = ReqSendDispls[i-1] + ReqSendCounts[i-1];
                }
            }
            int TotalReqRecv = ReqRecvDispls[size-1] + ReqRecvCounts[size-1];
            int TotalReqSend = ReqSendDispls[size-1] + ReqSendCounts[size-1];
            auto ReqRecvBuffer = new int[TotalReqRecv][4];
            auto ReqSendBuffer = new int[TotalReqSend][4];

            MPI_Alltoallv(ReqSendBuffer, &ReqSendCounts[0], &ReqSendDispls[0], tuple_type_4 , ReqRecvBuffer, &ReqRecvCounts[0], &ReqRecvDispls[0], tuple_type_4, MPI_COMM_WORLD);
            
            vector<vector<vector<int> > > ProcessResponseBuffer(size);
            vector <int> ResponseSendCounts(size, 0);
            vector <int> ResponseRecvCounts(size, 0);

            for(int i = 0; i < TotalReqRecv; i++){
                int u = ReqRecvBuffer[i][0];
                int v = ReqRecvBuffer[i][1];
                int w = ReqRecvBuffer[i][2];
                int id = ReqRecvBuffer[i][3];
                // if(rank==0)
                // cout<<u<<" "<< v<<" "<<w<<endl;
                if(id == 1){
                    if(process_graph[v].find(w) != process_graph[v].end()){
                        support_map[make_pair(v,w)]++;
                        triangle_map[make_pair(v,w)].insert(u);
                        int ResponseProcess = process_map[u];
                        vector<int> ResponseTuple = {u,v,w,2};
                        ProcessResponseBuffer[ResponseProcess].push_back(ResponseTuple);
                        ResponseSendCounts[ResponseProcess]++;
                    }
                }
                else{
                    if(process_graph[w].find(v) != process_graph[w].end()){
                        support_map[make_pair(w,v)]++;
                        triangle_map[make_pair(w,v)].insert(u);
                        int ResponseProcess = process_map[u];
                        vector <int> ResponseTuple = {u,v,w,2};
                        ProcessResponseBuffer[ResponseProcess].push_back(ResponseTuple);
                        ResponseSendCounts[ResponseProcess]++;
                    }
                }
            }
            // vector <int> ResponseRecvSize(size, 0);
            MPI_Alltoall(&ResponseSendCounts[0], 1, MPI_INT, &ResponseRecvCounts[0], 1, MPI_INT, MPI_COMM_WORLD);

            // vector<vector<int>> ResponseSendBuffer;
            vector <int> ResponseSendDispls(size, 0);
            vector <int> ResponseRecvDispls(size, 0);
            for(int i = 0; i< size; i++){
                if(i!=0){
                    ResponseSendDispls[i] = ResponseSendDispls[i-1] + ResponseSendCounts[i-1];
                    ResponseRecvDispls[i] = ResponseRecvDispls[i-1] + ResponseRecvCounts[i-1];
                }
            }
            int TotalResponseSend = ResponseSendDispls[size-1] + ResponseSendCounts[size-1];
            int TotalResponseReceived = ResponseRecvDispls[size-1] + ResponseRecvCounts[size-1];
            auto ResponseSendBuffer = new int[TotalResponseSend][4]; 
            auto ResponseRecvBuffer = new int[TotalResponseReceived][4];
            
            for(int i = 0; i < size; i++){
                for(int j = 0; j < ProcessResponseBuffer[i].size();j++){
                    int idx = ResponseSendDispls[i] + j;
                    for(int k = 0; k < 4; k++){
                        ResponseSendBuffer[idx][k] = ProcessResponseBuffer[i][j][k];
                    }
                }
            }
            MPI_Alltoallv(ResponseSendBuffer, &ResponseSendCounts[0], &ResponseSendDispls[0], tuple_type_4, ResponseRecvBuffer, &ResponseRecvCounts[0], &ResponseRecvDispls[0], tuple_type_4, MPI_COMM_WORLD);
            //as the zombie process is not sending any data, we need to check if the ReqRecvBuffer is empty or not
            if(TotalResponseReceived!=0){
                if(rank == 0)
                   cout<<"entered zombie zone"<<endl;
            }
        }
        //initialised the support map and triangle map
    
        else{
            for (auto it = process_graph.begin(); it != process_graph.end(); it++){
                int u = it->first;
                int degree_u = it->second.size();
                // cout<<"current vertex: "<<u<<endl;
                for ( auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                    int v = *it2;
                    // counter++;
                    //find the owner of edge <u,v>
                    int offset = 0;
                    headerfile.seekg(4*v, ios::beg);
                    headerfile.read((char*)&offset, 4);

                    int degree_v = 0;
                    infile.seekg(offset+4, ios::beg);
                    infile.read((char*)&degree_v, 4);

                    //check if the edge belongs to the current process 
                    bool flag = cmp_for_vertexSort(make_pair(u, degree_u), make_pair(v, degree_v));
                    if (!flag){
                        continue;
                    }
                    // if (rank ==0)
                    //     cout<< edge_count<<" <"<<u<<","<<v<<"> "<<endl; 
                    edge_count++;
                    //create a temporary send buffer and recv buffer for MPI_AlltoAllv
                    vector <int> ReqSendCounts(size, 0);
                    vector <int> ReqRecvCounts(size, 0);
                    vector <int> ReqSendDispls(size, 0);
                    vector <int> ReqRecvDispls(size, 0);
                    vector< vector< vector< int >>> ProcessReqBuffer(size);
                    for(auto it3 = it2; it3 != it->second.end(); it3++){
                        int w = *it3;
                        if (w == v){
                            continue;
                        }
                        int offset = 0;
                        headerfile.seekg(4*w, ios::beg);
                        headerfile.read((char*)&offset, 4);

                        int degree_w = 0;
                        infile.seekg(offset+4, ios::beg);
                        infile.read((char*)&degree_w, 4);

                        // check that <u,v> and <u,w> are monotonous wedge 
                        bool flag2 = cmp_for_vertexSort(make_pair(u, degree_u), make_pair(w, degree_w));
                        if (!flag2){
                            continue;
                        }
                        // if(rank == 0)
                        //     cout<< w<<endl;
                        //find the owner of edge <v,w>
                        bool flag3 = cmp_for_vertexSort(make_pair(v, degree_v), make_pair(w, degree_w));
                        int process_q;
                        vector<int> req_tuple;
                        if (!flag3){
                            process_q = process_map[w];
                            req_tuple = {u,v,w,2};
                        }
                        else{
                            process_q = process_map[v];
                            req_tuple = {u,v,w,1};
                        }
                        ProcessReqBuffer[process_q].push_back(req_tuple);
                        ReqSendCounts[process_q]++;
                    }
                    MPI_Alltoall(&ReqSendCounts[0], 1, MPI_INT, &ReqRecvCounts[0], 1 , MPI_INT, MPI_COMM_WORLD);
                    int TotalReqSend = 0;
                    int TotalReqReceived = 0;
                    for(int i = 0; i < size; i++){
                        TotalReqSend += ReqSendCounts[i];
                        TotalReqReceived += ReqRecvCounts[i]; 
                        if(i != 0){
                            ReqSendDispls[i] = ReqSendDispls[i-1] + ReqSendCounts[i-1];
                            ReqRecvDispls[i] = ReqRecvDispls[i-1] + ReqRecvCounts[i-1];
                        }
                    }
                    auto ReqRecvBuffer = new int[TotalReqReceived][4];
                    auto ReqSendBuffer = new int[TotalReqSend][4];
                    for(int i = 0; i < size; i++){
                        for(int j = 0; j < ProcessReqBuffer[i].size(); j++){
                            int idx = ReqSendDispls[i] + j;
                            for(int k = 0; k < 4; k++){
                                ReqSendBuffer[idx][k] = ProcessReqBuffer[i][j][k];
                            }   
                        }
                    }
                    // if(rank==0)
                        // cout<<endl;
                    MPI_Alltoallv(ReqSendBuffer, &ReqSendCounts[0], &ReqSendDispls[0], tuple_type_4 , ReqRecvBuffer, &ReqRecvCounts[0], &ReqRecvDispls[0], tuple_type_4, MPI_COMM_WORLD);

                    //now we have the recv buffer
                    //now we need to process the recv buffer
                    vector <int> ResponseSendCounts(size, 0);
                    vector <int> ResponseRecvCounts(size, 0);
                    vector <int> ResponseSendDispls(size, 0);
                    vector <int> ResponseRecvDispls(size, 0);
                    vector< vector< vector < int >>> ProcessResponseBuffer(size);
                    for(int i = 0; i < TotalReqReceived; i++){
                        // vector< int> ReqTuple = ReqRecvBuffer[i];
                        int Req_u = ReqRecvBuffer[i][0];
                        int Req_v = ReqRecvBuffer[i][1];
                        int Req_w = ReqRecvBuffer[i][2];
                        int id = ReqRecvBuffer[i][3];
                        // if(rank == 0)
                        //     cout<<"requested tuple" <<Req_u <<" "<<Req_v<<" "<<Req_w<< " "<< id<<endl;
                        if (id == 1){
                            if (process_graph[Req_v].find(Req_w)!=process_graph[Req_v].end()){
                                support_map[make_pair(Req_v, Req_w)]++;
                                triangle_map[make_pair(Req_v, Req_w)].insert(Req_u);
                                int ResponseProcess = process_map[Req_u];
                                vector<int> response_tuple= {Req_u, Req_v, Req_w, 1};
                                ProcessResponseBuffer[ResponseProcess].push_back(response_tuple);
                                ResponseSendCounts[ResponseProcess]++;
                            }
                        }
                        else{
                            if (process_graph[Req_w].find(Req_v)!=process_graph[Req_w].end()){
                                support_map[make_pair(Req_w, Req_v)]++;
                                triangle_map[make_pair(Req_w, Req_v)].insert(Req_u);
                                int ResponseProcess = process_map[Req_u];
                                vector<int> response_tuple = {Req_u, Req_v, Req_w, 2};
                                ProcessResponseBuffer[ResponseProcess].push_back(response_tuple);
                                ResponseSendCounts[ResponseProcess]++;
                            }
                        }
                    }
                    MPI_Alltoall(&ResponseSendCounts[0], 1, MPI_INT, &ResponseRecvCounts[0], 1, MPI_INT, MPI_COMM_WORLD);
                    // vector <vector<int> > ResponseSendBuffer;
                    for(int i = 0; i < size; i++){
                        if(i != 0){
                            ResponseSendDispls[i] = ResponseSendDispls[i-1] + ResponseSendCounts[i-1];
                            ResponseRecvDispls[i] = ResponseRecvDispls[i-1] + ResponseRecvCounts[i-1];
                        }
                    }
                    int TotalResponseReceived = ResponseRecvCounts[size-1] + ResponseRecvDispls[size-1];
                    int TotalResponseSend = ResponseSendCounts[size-1] + ResponseSendDispls[size-1];
                    auto ResponseRecvBuffer = new int[TotalResponseReceived][4];
                    auto ResponseSendBuffer = new int[TotalResponseSend][4];
                    for(int i = 0; i < size; i++){
                        for(int j = 0; j < ProcessResponseBuffer[i].size(); j++){
                            int idx = ResponseSendDispls[i] + j;
                            for(int k = 0; k < 4; k++){
                                ResponseSendBuffer[idx][k] = ProcessResponseBuffer[i][j][k];
                            }   
                        }
                    }
                    // vector <vector<int> > ResponseRecvBuffer(totalResponseReceived, vector<int>(4,0));
                    MPI_Alltoallv(ResponseSendBuffer, &ResponseSendCounts[0], &ResponseSendDispls[0], tuple_type_4 , ResponseRecvBuffer, &ResponseRecvCounts[0], &ResponseRecvDispls[0], tuple_type_4, MPI_COMM_WORLD);
                    // //now we have the result of the query in ResponseRecvBuffer
                    // //now we need to process the recv buffer
                    for(int i = 0; i < TotalResponseReceived; i++){
                        int u2 = ResponseRecvBuffer[i][0];
                        int v2 = ResponseRecvBuffer[i][1];
                        int w2 = ResponseRecvBuffer[i][2];
                        int id2 = ResponseRecvBuffer[i][3];
                        support_map[make_pair(u2, v2)]++;
                        triangle_map[make_pair(u2, v2)].insert(w2);
                        support_map[make_pair(u2, w2)]++;
                        triangle_map[make_pair(u2, w2)].insert(v2);
                    }
                }    
            }
        }
    }
    infile.close();
    headerfile.close();
    // for (auto it = triangle_map.begin(); it != triangle_map.end(); it++){
    //     cout << "Process " << rank << ": " << it->first.first << " " << it->first.second << ": " << it->second.size()<<" :";
    //     for(auto it1 = it->second.begin(); it1!= it->second.end(); it1++){
    //         cout<< *it1 << " ";
    //     }
    //     cout<<endl;
    // }

    
    int sizeOfMap = support_map.size();
    int sizeOfArray = sizeOfMap*3;

    // first send(gather) the size of the arraySupport
    int* recvdArraySizes;
    if (rank == 0)
    {
        recvdArraySizes = new int[size];
    }
    MPI_Gather(&sizeOfArray, 1, MPI_INT, recvdArraySizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* recvdArraySupport;
    int sizeOfRecvdArraySupport = 0;
    if(rank == 0){
        for (int it = 0; it < size; it++){
            sizeOfRecvdArraySupport+=recvdArraySizes[it];
        }
        recvdArraySupport = new int[sizeOfRecvdArraySupport];
    }
    
    int* recvDisplsG;
    if (rank==0)
    {
        recvDisplsG = new int[size];
        int offset = 0;
        for (int i = 0; i < size; i++)
        {
            recvDisplsG[i] = offset;
            offset+=recvdArraySizes[i];
        }
    }

    int* arrayOfSupport = new int[sizeOfMap*3];
    convertSupportMap(support_map, arrayOfSupport);

    MPI_Gatherv(arrayOfSupport, sizeOfArray, MPI_INT, recvdArraySupport, recvdArraySizes, recvDisplsG, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    // cout<<"hi"<<endl;

    //////////////////////  send triangle map to process 0   ///////////////////////////
    if (rank!=0)
    {
        int sizeOfTriangleMap = triangle_map.size();
        int totalSizeOfTriangleMap = 3*sizeOfTriangleMap;
        
        for (auto it = triangle_map.begin(); it != triangle_map.end(); it++){
            totalSizeOfTriangleMap += it->second.size();
        }

        int* triangleMapArray = new int[totalSizeOfTriangleMap];
        convertTriangleMap(triangle_map, triangleMapArray, totalSizeOfTriangleMap);

        // MPI_Send(&totalSizeOfTriangleMap, 1, MPI_INT, 0, 0, MPI_COMM_WORLD)

        MPI_Send(triangleMapArray ,totalSizeOfTriangleMap, MPI_INT, 0, rank, MPI_COMM_WORLD);
        // cout<<rank<<" send"<<endl;
        // MPI_Barrier(MPI_COMM_WORLD);
    }

    unordered_map<int,set<int>> graph;
    map <pair<int, int>, int> collectedSupport;
    map<pair<int,int>, set<int>> collectedTriangles;

    if (rank==0)
    {   

        int nn = 0, mm = 0;
        // unordered_map<int<set<int>> graph;
        readGraph(nn, mm, graph, iInputPath, iHeaderPath);


        // convert the recvd support array to map
        // map <pair<int, int>, int> collectedSupport;
        for (int i = 0; i < sizeOfRecvdArraySupport; i+=3)
        {
            int u = recvdArraySupport[i];
            int v = recvdArraySupport[i+1];
            int w = recvdArraySupport[i+2];
            pair <int, int> edge = make_pair(u, v);
            collectedSupport[edge] = w;
        }
        delete[] recvdArraySupport;

        // map<pair<int,int>, set<int>> collectedTriangles;

        for(int trank=1; trank<size; trank++)
        {
            int sizeToRecv;
            MPI_Status status;

            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &sizeToRecv);

            int* buffer = new int[sizeToRecv];
            int source = status.MPI_SOURCE;
            int tag = status.MPI_TAG;

            MPI_Recv(buffer, sizeToRecv, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // cout<<source<<" received"<<endl;
            // int offset = 0;
            int j = 0;
            while (j < sizeToRecv)
            {
                int u = buffer[j];      j++;
                int v = buffer[j];      j++;

                set<int> s;
                while (buffer[j] != -1)
                {
                    s.insert(buffer[j]);
                    j++;
                }
                
                pair<int,int> edge = make_pair(u, v);
                collectedTriangles[edge] = s;
                j++;
            }

            delete[] buffer;
        }
        
        // now we have the collected triangles in collectedTriangles
        // add the triangle support of this process (pid=0) to the collectedTriangles
        for (auto it = triangle_map.begin(); it != triangle_map.end(); it++)
        {
            pair<int,int> edge = it->first;
            set<int> s = it->second;
            collectedTriangles[edge] = s;
        }
        // MPI_Barrier(MPI_COMM_WORLD);

    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    if (iTaskId == 1){
        if(size == 1){
            if(rank == 0){
                int num_t = 1;
                vector<string> toPrint(iendK-istartK+1);
                ofstream outfile(iOutputPath);
                for(int k = istartK; k <= iendK; k++){
                    // cout<<"k = "<<k<<": ";
                    //calculate the truss value in rank 0
                    vector<pair<int,int>> deletable;
                    Initialize(k+2, deletable, collectedSupport);
            
                    FilterEdges(graph, k+2, deletable, collectedSupport, collectedTriangles, num_t);
            
                    // findKTruss(graph, k+2);
                    
                    modifyGraph(graph);

                    // printGraph(graph);
                    // cout<<graph.size()<<endl;
                    string s;
                    if (iVerbose==0){
                        s = "1";
                        if (graph.size() == 0){
                            s = "0";
                        }
                    }
                    else{
                        s = "1\n";
                        if (graph.size() == 0){
                            s = "0\n";
                        }
                        else{
                            vector<unordered_set<int>> connectedComponents = findConnectedComponents(graph);
                            
                            int noOfGroups = connectedComponents.size();
                            s = s + to_string(noOfGroups) + "\n";

                            sort(connectedComponents.begin(), connectedComponents.end(), [](unordered_set<int> a, unordered_set<int> b){return a.size() > b.size();});

                            for (int i=0; i<noOfGroups; i++){
                                vector<int> temp_vec(connectedComponents[i].begin(), connectedComponents[i].end());
                                sort(temp_vec.begin(), temp_vec.end());
                                for (auto it=temp_vec.begin(); it!=temp_vec.end(); it++){
                                    s = s + to_string(*it) + " ";
                                }
                                s = s + "\n";
                            }
                        }
                    }
                    outfile<< s;

                    // cout<< "hii"<<endl;
                    // MPI_Send(&k, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                    // int str_size = s.size();
                    // char send_str[str_size];
                    // strncpy(send_str, s.c_str(), str_size);
                    // MPI_Send(&str_size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                    // MPI_Send(send_str, str_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                    // cout<<"send for k "<<endl;
                } 
                outfile.close();
            }
            // cout<<graph.size()<<endl;
        }
        else{
            if(rank == 0){
                int num_t = 1;
                ofstream outfile(iOutputPath);
                for(int k = istartK; k <= iendK; k++){
                    // cout<<"k = "<<k<<": ";
                    //calculate the truss value in rank 0
                    vector<pair<int,int>> deletable;
                    Initialize(k+2, deletable, collectedSupport);
            
                    FilterEdges(graph, k+2, deletable, collectedSupport, collectedTriangles, num_t);
            
                    // findKTruss(graph, k+2);
                    
                    modifyGraph(graph);

                    // printGraph(graph);
                    // cout<<graph.size()<<endl;
                    string s;
                    if (iVerbose==0){
                        s = "1";
                        if (graph.size() == 0){
                            s = "0";
                        }
                    }
                    else{
                        s = "1\n";
                        if (graph.size() == 0){
                            s = "0\n";
                        }
                        else{
                            vector<unordered_set<int>> connectedComponents = findConnectedComponents(graph);
                            
                            int noOfGroups = connectedComponents.size();
                            s = s + to_string(noOfGroups) + "\n";

                            sort(connectedComponents.begin(), connectedComponents.end(), [](unordered_set<int> a, unordered_set<int> b){return a.size() > b.size();});

                            for (int i=0; i<noOfGroups; i++){
                                vector<int> temp_vec(connectedComponents[i].begin(), connectedComponents[i].end());
                                sort(temp_vec.begin(), temp_vec.end());
                                for (auto it=temp_vec.begin(); it!=temp_vec.end(); it++){
                                    s = s + to_string(*it) + " ";
                                }
                                s = s + "\n";
                            }
                        }
                    }
                    cout<< "hii"<<endl;
                    MPI_Send(&k, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                    int str_size = s.size();
                    char send_str[str_size];
                    strncpy(send_str, s.c_str(), str_size);
                    MPI_Send(&str_size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                    MPI_Send(send_str, str_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                    cout<<"send for k "<<endl;
                } 
            }
            else if(rank == 1){ //if rank != 0 and task == 1
                vector<string> toPrint(iendK-istartK+1);
                for(int k = istartK; k <= iendK; k++){
                    int k_recv, str_size;
                    MPI_Status status;

                    MPI_Recv(&k_recv, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&str_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

                    char recv_str[str_size];
                    MPI_Recv(recv_str, str_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
                    string s(recv_str, str_size);
                    toPrint[k-istartK] = s;
                }
                writeGraph(toPrint, iOutputPath, iVerbose);
            }
        } 
    }
    else{
        unordered_map <int, int> vertexTrussMap;
        int* recvInfluencerArray;
        int recvSizeForInfluencerSum;
        int* recvDisplsForInfluencer; 
        vector<unordered_set<int>> connectedComponents;
        if(rank == 0)
        {   
            int num_t = 1;
            int k = iendK;
            vector<pair<int,int>> deletable;
            Initialize(k+2, deletable, collectedSupport);
    
            FilterEdges(graph, k+2, deletable, collectedSupport, collectedTriangles, num_t);
    
            modifyGraph(graph);

            connectedComponents = findConnectedComponents(graph);
            int noOfGroups = connectedComponents.size();

            //broadcast the connected components to all the processes
            MPI_Bcast(&noOfGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            int* sendSizesArray = new int[noOfGroups];
            vector<int> sendArray;
            int TotalVerticesSize = 0;
            for (int i=0; i<connectedComponents.size(); i++){
                sendSizesArray[i] = connectedComponents[i].size();
                for (auto it2=connectedComponents[i].begin(); it2!=connectedComponents[i].end(); it2++){
                    // if(i == 0){
                    //     cout<< *it2<<" ";
                    // }
                    sendArray.push_back(*it2);
                    vertexTrussMap[*it2] = i;
                    TotalVerticesSize++;
                }
            }
            // cout<<endl;
            for (int i=1; i<size; i++){
                MPI_Send(sendSizesArray, noOfGroups, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&sendArray[0], TotalVerticesSize, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            // cout<<"send sizes array"<<endl;
            // for(int i = 0; i< noOfGroups; i++){
            //     cout<< sendSizesArray[i]<<" ";
            // }
            // cout<<endl;
            // unordered_map <int, int> vertexTrussMap;
            // int idx = 0;
            // for(int trussValue = 0; trussValue < noOfGroups; trussValue++)
            // {
            //     for(int i=0; i<sendSizesArray[trussValue]; i++){
            //         vertexTrussMap[sendArray[idx]] = trussValue;
            //         idx++;
            //     }
            // }
        }
        if (rank!=0)
        {
            vector<string> toPrint;
            int noOfGroups;
            MPI_Bcast(&noOfGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            int* recvSizesArray = new int[noOfGroups];
            MPI_Status status;
            MPI_Recv(recvSizesArray, noOfGroups, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            
            int TotalVerticesSize = 0;
            for (int i=0; i<noOfGroups; i++){
                TotalVerticesSize += recvSizesArray[i];
            }
            // if(rank==1){
            //     cout<<"recv sizes array"<<endl;
            //     for(int i = 0; i< noOfGroups; i++){
            //         cout<< recvSizesArray[i]<<" ";
            //     }
            //     cout<<endl;
            // }
            int* recvArray = new int[TotalVerticesSize];
            MPI_Recv(recvArray, TotalVerticesSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            
            // unordered_map <int, int> vertexTrussMap;
            int idx = 0;
            for(int trussValue = 0; trussValue < noOfGroups; trussValue++){
                for(int i=0; i<recvSizesArray[trussValue]; i++){
                    vertexTrussMap[recvArray[idx]] = trussValue;
                    idx++;
                }
            }
            delete[] recvSizesArray;
            delete[] recvArray;
        }
        MPI_Barrier(MPI_COMM_WORLD);


        // compuation done 
        unordered_map <int, set<int>> influencerMap;
        for(auto it = process_graph.begin(); it != process_graph.end(); it++){
            int vertex = it->first;
            influencerMap[vertex] = set<int>();
            if(vertexTrussMap.find(vertex)!=vertexTrussMap.end()){
                influencerMap[vertex].insert(vertexTrussMap[vertex]);
            }
            for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                int neighbour = *it2;
                if(vertexTrussMap.find(neighbour)!=vertexTrussMap.end())
                    influencerMap[vertex].insert(vertexTrussMap[neighbour]);
            }
        }

        int numVertices = influencerMap.size();
        int sendSizeForInfluencer = 0;
        vector <int> sendInfluencerArray;
        for(auto it = influencerMap.begin(); it != influencerMap.end(); it++)
        {
            int influencerVertex = it->first;
            if(it->second.size() < iP){
                numVertices--;
                continue;
            }
            sendInfluencerArray.push_back(influencerVertex);
            sendSizeForInfluencer++;
            sendInfluencerArray.push_back(it->second.size());
            sendSizeForInfluencer++;
            for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                sendInfluencerArray.push_back(*it2);
                sendSizeForInfluencer++;
            }
        }
        // cout<<rank<<" send influencer vertices "<<numVertices<<endl;

        int* recvSizeForInfluencer;
        if (rank==0)
        {
            recvSizeForInfluencer =  new int[size];
        } 
        
        
        MPI_Gather(&sendSizeForInfluencer, 1, MPI_INT, recvSizeForInfluencer, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);


        if (rank==0)
        {
            
            recvSizeForInfluencerSum = 0;
            int offset = 0;
            recvDisplsForInfluencer = new int[size];
            for(int i=0; i<size; i++){
                recvDisplsForInfluencer[i] = offset;
                offset += recvSizeForInfluencer[i];
                recvSizeForInfluencerSum += recvSizeForInfluencer[i];
            }

            recvInfluencerArray = new int[recvSizeForInfluencerSum];
        }
        int * sendInfArray = new int[sendSizeForInfluencer];
        for(int i = 0; i < sendInfluencerArray.size(); i++){
            sendInfArray[i] = sendInfluencerArray[i];
        }

        MPI_Gatherv(sendInfArray, sendSizeForInfluencer, MPI_INT, recvInfluencerArray, recvSizeForInfluencer, recvDisplsForInfluencer, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);


        if (rank==0)
        {
            unordered_map <int, set<int>> influencerMapRecv;
            convertInfluencerArrayToMap(recvInfluencerArray, recvSizeForInfluencerSum, influencerMapRecv);
            
            delete[] recvSizeForInfluencer;
            delete[] recvDisplsForInfluencer;
            delete[] recvInfluencerArray;


            for(auto it = influencerMap.begin(); it != influencerMap.end(); it++)
            {
                int vertex = it->first;
                if(it->second.size() < iP){
                    continue;
                }
                influencerMapRecv[vertex] = set<int>();
                for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                    int neighbour = *it2;
                    influencerMapRecv[vertex].insert(neighbour);
                }
            }

            string s;
            if(iVerbose == 0)
            {
                s = "0 \n";
                if(influencerMapRecv.size() > 0){
                    s = to_string(influencerMapRecv.size()) + "\n";
                    for(auto it = influencerMapRecv.begin(); it != influencerMapRecv.end(); it++){
                        s += to_string(it->first) + " ";
                    }
                    s += "\n";
                }
            }
            else
            {
                s = "0 \n";
                if(influencerMapRecv.size() > 0){
                    s = to_string(influencerMapRecv.size()) + "\n";
                    for(auto it = influencerMapRecv.begin(); it != influencerMapRecv.end(); it++){
                        s += to_string(it->first) + "\n";
                        for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++){
                            int component_no = *it2;
                            // if(it->first == 68752)
                            //     cout<< *it2<<": ";
                            for(auto it3 = connectedComponents[component_no].begin(); it3!= connectedComponents[component_no].end(); it3++){
                                // if(it->first == 68752)
                                //     cout<< *it3 << " ";
                                s += to_string(*it3) + " ";
                            }
                            // cout<<endl;
                        }
                        s += "\n";
                    }
                }
            }

            ofstream outfile(iOutputPath);
            outfile<<s;
            outfile.close();
             
        }
    }

    MPI_Finalize();
    return 0;
}