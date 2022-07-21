/*
This is code associated with the paper "Graph coloring for parallelization of quantum simulation of fermionic systems." (Bringewatt and Davoudi).

Author: Jacob Bringewatt
Affiliations: Joint Center for Quantum Information and Computer Science (QuICS); Joint Quantum Institute (JQI); University of Maryland, College Park; 
Date: 2022
*/
/*----------------------------------------------------------------------------------------------------------------*/

/* Include statements*/
#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
//Boost libraries
#include <boost/config.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/edge_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/utility.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/iterator/advance.hpp>

/*----------------------------------------------------------------------------------------------------------------*/

/*Structure of global parameters*/
struct Parameters{
    bool strong_coloring=false;        //determines if we're doing strong coloring or not
    bool verbose =false;                //determines if verbose output or not  
    int interactions;            //gives interaction type (0=all2all, 1=custom, reads in from file)
    bool output_conflict_graph;  //output conflict graph as graphviz file
    int samples =10;                 //number of samples
};
//initialize global parameters with default values
struct Parameters params={
    .strong_coloring=false, 
    .verbose=false,
    .interactions=0,
    .output_conflict_graph=false
};

/*----------------------------------------------------------------------------------------------------------------*/

/*Define system graph*/
struct SystemGraphVertexProperties{
    int node_id; //vertex id 
    std::string label; //vertex name
    std::string type; //keeps track if vertex is virtual or physical
    std::vector<int> paths_contained_in; //lists paths the vertex is contained in
    std::vector<int> enumeration; //keeps track of enumeration of edges incident on vertex
    int max_enum; //keeps track of maximum enumerated vertex for strong coloring
    int min_enum; //keeps track of minimum enumerated vertex for strong coloring
};
struct SystemGraphEdgeProperties{
    int edge_id;
    int weight=1; //weight for edge 
    int path_cnt=0; //keeps track of how many times an edge is used in a path
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, SystemGraphVertexProperties, SystemGraphEdgeProperties> SystemGraph;

/*----------------------------------------------------------------------------------------------------------------*/
/*Define conflict graph*/
struct ConflictGraphVertexProperties{
    int node_id; //vertex id
    std::string label;
    std::string interaction_type;
    std::string fillcolor; //color of conflict graph vertex
    std::string style="filled"; //fill vertices of conflict graph if strong coloring
    std::string colorscheme="paired12";
    std::vector<int> path; //gives path corresponding to vertex as vector of system graph vertices
};
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, ConflictGraphVertexProperties> ConflictGraph;
/*Coloring typedefs*/
typedef boost::graph_traits<ConflictGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<ConflictGraph>::vertices_size_type vertices_size_type;
typedef boost::property_map<ConflictGraph, boost::vertex_index_t>::const_type vertex_index_map;
typedef boost::iterator_property_map<vertices_size_type*, vertex_index_map>ColorMap;
typedef boost::iterator_property_map<int*, vertex_index_map> DegreeMap;

/*----------------------------------------------------------------------------------------------------------------*/
/*Define interaction struct*/
struct Interaction{
    std::vector<int> interacting_vertices; //list of vertices involved as determined by node_id
    std::string interaction_type; //interaction type, 'AB', 'B'
};

/*----------------------------------------------------------------------------------------------------------------*/
/*Define interaction graph*/
struct InteractionGraphVertexProperties{
    int node_id; //vertex id 
};
struct InteractionGraphEdgeProperties{
    int edge_id;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, InteractionGraphVertexProperties, InteractionGraphEdgeProperties> InteractionGraph;




/*Main, input to program*/
int main(int argc, char* argv[]){

    //initialize random engine
    std::random_device rd;
    std::default_random_engine rng(rd());


    //define char arrays for files
    char *interactions_filename = NULL;
    char *system_graph_filename = NULL;
    char *output_data_filename = NULL;

    /*Parse input*/
    int c;
    while ((c = getopt (argc, argv, "hvrgi:s:o:n:")) != -1)
    switch (c)
      {
      case 'h':
        std::cerr << "Usage: ./main <options>\n"
            << "<options> \n"
            << "-s <filename>\t filename for system graph (default is K_5)\n"  
            << "-o <filename>\t filename for outputfile (default is output)\n"  
            << "-i <filename>\t filename for interactions (default is all-to-all)\n" 
            << "-v\t Gives extra output updating on progress of program.\n"
            << "-r\t Sets to strong coloring problem (default is weak coloring).\n"
            << "-g\t Outputs conflict graph as .dot file.\n"
            << "-n\t Number of samples. (default is 10)"
            << std::endl;
            return 1;
        break;
      case 'v':
        params.verbose = true;
        break;
      case 'r':
        params.strong_coloring = true;
        break;
      case 'g':
        params.output_conflict_graph = true;
        break;
      case 'i':
        interactions_filename=optarg;
        break;
      case 's':
        system_graph_filename=optarg;
        break;
      case 'o':
        output_data_filename=optarg;
        break;
      case 'n':
        params.samples=std::stoi(optarg);
        break;
      case '?':
        if (optopt == 'i' || optopt == 's' || optopt=='o')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();      
    }
    
    /**************************************************************************************
    /*Create system graph*/
    if(params.verbose) std::cout<<"Initializing system graph..."<<std::endl;
    SystemGraph system_graph(0);
    //get dynamic properties
    boost::dynamic_properties system_graph_dynamic_properties(boost::ignore_other_properties);
    system_graph_dynamic_properties.property("node_id", get(&SystemGraphVertexProperties::node_id, system_graph));
    system_graph_dynamic_properties.property("label", get(&SystemGraphVertexProperties::label, system_graph));
    system_graph_dynamic_properties.property("type", get(&SystemGraphVertexProperties::type, system_graph));
    /*Read in system graph*/
    std::ifstream file_stream;
    int status;
    if(system_graph_filename==NULL){
        file_stream.open("src/K5.dot");      
    }else{
        file_stream.open(system_graph_filename);
    }
    std::string graph_string2((std::istreambuf_iterator<char>(file_stream)),(std::istreambuf_iterator<char>()));
    try{
        status=boost::read_graphviz(graph_string2, system_graph, system_graph_dynamic_properties, "node_id"); 
    }catch(...){
        std::cout << "Reading system graph failed. Check file exists and is a graphviz file."<< std::endl; 
        return 0;
    }   
    if(!status){
        std::cout << "Reading system graph failed. Check file exists and is a graphviz file."<< std::endl; 
        return 0;
    }  
    
    /**************************************************************************************/
    //read in interactions
    if(params.verbose) std::cout<<"Initializing interactions..."<<std::endl;
    std::vector<Interaction> interaction_list;
    if(interactions_filename==NULL){
        //interactions default to all-to-all between physical vertices
        for (auto u : boost::make_iterator_range(boost::vertices(system_graph))) {
            if(system_graph[u].type=="virtual") continue;
            for (auto v : boost::make_iterator_range(boost::vertices(system_graph))) {
                if(system_graph[v].type=="virtual") continue;
                if(u==v){
                    std::vector<int> interactions{system_graph[u].node_id};
                    interaction_list.push_back(Interaction{interactions, "B"}); 
                }else{
                    std::vector<int> interactions{system_graph[u].node_id,system_graph[v].node_id};
                    interaction_list.push_back(Interaction{interactions, "AB"});
                }
            }
        }
    }else{
      /* Interaction graph includes same vertices as system graph but no edges should exist between virtual vertices*/
      InteractionGraph interaction_graph(0);
      //get dynamic properties
      boost::dynamic_properties interaction_graph_dynamic_properties(boost::ignore_other_properties);
      interaction_graph_dynamic_properties.property("node_id", get(&InteractionGraphVertexProperties::node_id, interaction_graph));
      std::ifstream file_stream2;
      std::cout << interactions_filename << std::endl;
      file_stream2.open(interactions_filename);      
      std::string graph_string3((std::istreambuf_iterator<char>(file_stream2)),(std::istreambuf_iterator<char>()));
      try{
          status=boost::read_graphviz(graph_string3, interaction_graph, interaction_graph_dynamic_properties, "node_id"); 
      }catch(...){
          std::cout << "Reading interaction graph failed. Check file exists and is a graphviz file."<< std::endl; 
          return 0;
      }   
      for (auto ep = boost::edges(interaction_graph); ep.first != ep.second; ++ep.first)
      {
          // Get the two vertices that are joined by this edge...
          int u=boost::source(*ep.first,interaction_graph);
          int v=boost::target(*ep.first,interaction_graph);
          if (u==v){
            std::vector<int> interactions{system_graph[u].node_id};
            interaction_list.push_back(Interaction{interactions, "B"}); 
          }else{
            std::vector<int> interactions{system_graph[u].node_id,system_graph[v].node_id};
            interaction_list.push_back(Interaction{interactions, "AB"});
          }
      }
    }

    //Uncomment to double check interaction list is read in properly
    // for(auto interaction : interaction_list){
    //   for(auto  v : interaction.interacting_vertices){
    //     std::cout<<v<< " ";
    //   }
    //   std::cout<< "| " << interaction.interaction_type<<std::endl;
    // }

    /**************************************************************************************/

    //Randomize interaction list order and try to find paths
    std::vector<int> chromatic_numbers;
    for(int i=0; i<params.samples; i++){  
      if(params.verbose) std::cout<<"Sample "<< i << std::endl;    
      //shuffle interactions
      std::shuffle(interaction_list.begin(), interaction_list.end(), rng);

      //reset system graph
      int cntr=0;
      for (auto e : boost::make_iterator_range(boost::edges(system_graph))){   
          system_graph[e].edge_id=cntr;       
          system_graph[e].weight=1;
          system_graph[e].path_cnt=0;
          cntr++;
      }
      for (auto v : boost::make_iterator_range(boost::vertices(system_graph))){  
          system_graph[v].paths_contained_in.clear();   
          system_graph[v].enumeration.clear();
          system_graph[v].enumeration.resize(boost::out_degree(v, system_graph));
          std::fill(system_graph[v].enumeration.begin(), system_graph[v].enumeration.end(), -1);
          system_graph[v].max_enum=boost::out_degree(v, system_graph)-1;
          system_graph[v].min_enum=0;     
          if(system_graph[v].type=="physical"){
               for (auto e : boost::make_iterator_range(boost::out_edges(v, system_graph))){
                   system_graph[e].weight+=10; //weight edges with physical vertices more to try and avoid them
               }
          }
      }
      
      //initialize conflict graph 
      ConflictGraph conflict_graph(interaction_list.size()); //initialize with interaction list size number of vertices

      //now loop through interactions
      int pathidx=0;
      for(auto interaction : interaction_list){
        if(params.verbose){
          std::cout << std::endl << "Interaction: " << pathidx+1 << " of " << interaction_list.size() << std::endl;
          std::cout << "Interacting vertices: ";
          for(auto v : interaction.interacting_vertices) std::cout << v <<" ";
          std::cout << std::endl;
        }
        //assign interaction type to conflict graph node
        conflict_graph[pathidx].interaction_type=interaction.interaction_type; 
        std::vector<int> mypath;
        if(interaction.interaction_type=="AB"){
            //find shortest path for interaction if of AB type
            std::vector<double> distances(boost::num_vertices(system_graph));
            std::vector<SystemGraph::vertex_descriptor> predecessors(boost::num_vertices(system_graph));
            SystemGraph::vertex_descriptor from = boost::vertex(interaction.interacting_vertices[0], system_graph); //this line assumes only single path interactions (i.e. 2 body terms)
            boost::dijkstra_shortest_paths(system_graph, from, boost::weight_map(get(&SystemGraphEdgeProperties::weight, system_graph)).distance_map(boost::make_iterator_property_map(distances.begin(), get(boost::vertex_index, system_graph))).predecessor_map(boost::make_iterator_property_map(predecessors.begin(), get(boost::vertex_index, system_graph))));
            // Extract the shortest path for this edge
            SystemGraph::vertex_descriptor  to = boost::vertex(interaction.interacting_vertices[1], system_graph);
            //push vertex back onto path 
            mypath.push_back( to );             
            //check if to has any other paths going through it, if yes add edge to conflict graph appropriately
            for(int other_path : system_graph[to].paths_contained_in){
              if(!params.strong_coloring){
                boost::add_edge(pathidx, other_path, conflict_graph);
              }
            }
            //mark path as going through vertex 
            system_graph[to].paths_contained_in.push_back(pathidx); 

            for(SystemGraph::vertex_descriptor u = predecessors[to]; u != to; to=u, u=predecessors[to])
            {
                //add vertex u to path
                mypath.push_back( u );

                //check if to has any other paths going through it, if yes add edge to conflict graph appropriately, for weak coloring
                for(int other_path : system_graph[u].paths_contained_in){
                  if(!params.strong_coloring){
                    boost::add_edge(pathidx, other_path, conflict_graph);
                  }
                }

                //mark path as going through vertex u
                system_graph[u].paths_contained_in.push_back(pathidx);
                std::pair<SystemGraph::edge_descriptor, bool> edge=boost::edge(u, to, system_graph);
                if(edge.second){
                    system_graph[edge.first].path_cnt++; //mark edge as used
                    system_graph[edge.first].weight+=3; //weight edge more for being used, this can be tuned
                    //if(system_graph[edge.first].enumeration==-1){
                    //   system_graph[edge.first].enumeration
                   // } //do strong coloring enumeration
                }
            }
            //reverse so path is in correct order
            std::reverse(mypath.begin(), mypath.end());
            conflict_graph[pathidx].label=system_graph[interaction.interacting_vertices[0]].label+","+system_graph[interaction.interacting_vertices[1]].label;
            
        }else if(interaction.interaction_type=="B"){           
            mypath.push_back(interaction.interacting_vertices[0]); //path is just a single vertex
            //check if to has any other paths going through it, if yes add edge to conflict graph appropriately
            for(int other_path : system_graph[interaction.interacting_vertices[0]].paths_contained_in){
              //if(!params.strong_coloring){ //this type of interaction is the same for both weak and strong coloring
                boost::add_edge(pathidx, other_path, conflict_graph);
              //}
            }
            //mark path as going through initial vertex
            system_graph[interaction.interacting_vertices[0]].paths_contained_in.push_back(pathidx);
            conflict_graph[pathidx].label=system_graph[interaction.interacting_vertices[0]].label;
        }
        //add edges to conflict graph for strong coloring. this is after the for loop building up the path because requires full path enumeration information
        if(params.strong_coloring && interaction.interaction_type=="AB"){
          //loop through vertices in my_path
          int cnt=0;
          for(int u : mypath){
            //get enumeration of edges in current path
            int incoming_edge_enum_curr;
            int outgoing_edge_enum_curr;
            if(cnt==0){//starting vertex
              incoming_edge_enum_curr=-99; //null=-99
              std::pair<SystemGraph::edge_descriptor, bool> edge_out=boost::edge(u, mypath[cnt+1], system_graph);
              SystemGraph::out_edge_iterator ei, ei_end;
              int idx=0;
              for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){
                if(system_graph[*ei].edge_id==system_graph[edge_out.first].edge_id){
                  if(system_graph[u].enumeration[idx]==-1){
                    system_graph[u].enumeration[idx]=system_graph[u].min_enum; //want to make enumeration as small as possible for outgoing vertex
                    system_graph[u].min_enum++;                    
                  }
                  outgoing_edge_enum_curr=system_graph[u].enumeration[idx];
                  break;
                }         
                idx++;       
              }
            }else if(cnt==mypath.size()-1){//ending vertex
              outgoing_edge_enum_curr=-99;
              std::pair<SystemGraph::edge_descriptor, bool> edge_in=boost::edge(mypath[cnt-1], u, system_graph);
              SystemGraph::out_edge_iterator ei, ei_end;
              int idx=0;
              for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){
                if(system_graph[*ei].edge_id==system_graph[edge_in.first].edge_id){
                  if(system_graph[u].enumeration[idx]==-1){
                    system_graph[u].enumeration[idx]=system_graph[u].max_enum; //want to make enumeration as large as possible for ingoing vertex
                    system_graph[u].max_enum--;                    
                  }
                  incoming_edge_enum_curr=system_graph[u].enumeration[idx];
                  break;
                }         
                idx++;       
              }
            }else{//middle vertex
              //we want to make enumerations as close as possible for ingoing/outgoing
              std::pair<SystemGraph::edge_descriptor, bool> edge_in=boost::edge(mypath[cnt-1], u, system_graph);
              std::pair<SystemGraph::edge_descriptor, bool> edge_out=boost::edge(u, mypath[cnt+1], system_graph);
              SystemGraph::out_edge_iterator ei, ei_end;
              int idx=0;
              int flag=0;
              for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){
                if(system_graph[*ei].edge_id==system_graph[edge_in.first].edge_id){
                  if(system_graph[u].enumeration[idx]==-1){
                    system_graph[u].enumeration[idx]=system_graph[u].min_enum; //want to make enumeration as small as possible for ingoing vertex
                    system_graph[u].min_enum++;                       
                  }
                  incoming_edge_enum_curr=system_graph[u].enumeration[idx];
                  flag++;
                }else if(system_graph[*ei].edge_id==system_graph[edge_out.first].edge_id){
                  if(system_graph[u].enumeration[idx]==-1){
                    system_graph[u].enumeration[idx]=system_graph[u].min_enum; //want to make enumeration as small as possible for outgoing vertex
                    system_graph[u].min_enum++;                    
                  }
                  outgoing_edge_enum_curr=system_graph[u].enumeration[idx];
                  flag++;
                }     
                if(flag==2) break; //enumeration of both edges found    
                idx++;       
              }
            }

            //search for possible intersections
            for(int other_path : system_graph[u].paths_contained_in){
              if(params.verbose) std::cout << "Comparing path (" << conflict_graph[pathidx].label << ") with other path (" << conflict_graph[other_path].label << ")" <<std::endl;
              if(other_path==pathidx){
                if(params.verbose) std::cout << "same path, continue" <<std::endl;
                continue;
              }
              if(conflict_graph[other_path].interaction_type=="B"){
                if(params.verbose) std::cout << "Other path is type B, adding edge" <<std::endl;
                boost::add_edge(pathidx, other_path, conflict_graph);
                continue;
              }
              //get full path associated with index other_path
              std::vector<int> path1=conflict_graph[other_path].path;
              //find this vertex of interest in the path
              int val;
              auto res = std::find (path1.begin(), path1.end(), u); 
              int index_in_path=res-path1.begin();       
              int incoming_edge_enum_other;
              int outgoing_edge_enum_other;
              //get enumerations for other path and check for conflicts
              if(res-path1.begin()==0){                
                //starting vertex
                if(params.verbose) std::cout << "Intersecting at start of other path, ";
                std::pair<SystemGraph::edge_descriptor, bool> edge_out=boost::edge(u, path1[index_in_path+1], system_graph);                
                //get enumeration of overlapping path
                incoming_edge_enum_other=-99;
                SystemGraph::out_edge_iterator ei, ei_end;
                int idx=0;
                for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){                     
                  if(system_graph[*ei].edge_id==system_graph[edge_out.first].edge_id){                    
                    outgoing_edge_enum_other=system_graph[u].enumeration[idx];
                    break;
                  }         
                  idx++;       
                }
                //std::cout << "Incoming edge enumeration: " << incoming_edge_enum_curr << ", outgoing edge enumeration: " << outgoing_edge_enum_other << std::endl;
                if((incoming_edge_enum_curr==-99) //curr is starting vertex
                || (outgoing_edge_enum_curr==-99  && std::floor(incoming_edge_enum_curr/2.0)<=std::floor(outgoing_edge_enum_other/2.0)) //curr is ending vertex
                || (std::floor(incoming_edge_enum_curr/2.0)<=std::floor(outgoing_edge_enum_other/2.0))){//curr is middle vertex
                  boost::add_edge(pathidx, other_path, conflict_graph);
                  if(params.verbose) std::cout << "adding edge" <<std::endl;
                }    
              }else if(res-path1.begin()==path1.size()-1){
                //ending vertex
                if(params.verbose) std::cout << "Intersecting at end of other path, ";
                std::pair<SystemGraph::edge_descriptor, bool> edge_in=boost::edge(path1[res-path1.begin()-1], u, system_graph); 
                //get enumeration of overlapping path
                SystemGraph::out_edge_iterator ei, ei_end;
                int idx=0;
                for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){
                  if(system_graph[*ei].edge_id==system_graph[edge_in.first].edge_id){                    
                    outgoing_edge_enum_other=system_graph[u].enumeration[idx];
                    break;
                  }         
                  idx++;       
                }
                outgoing_edge_enum_other=-99;   
                if((incoming_edge_enum_curr==-99 && std::floor(outgoing_edge_enum_curr/2.0)>=std::floor(incoming_edge_enum_other/2.0)) //curr is starting vertex
                || (outgoing_edge_enum_curr==-99) //curr is ending vertex
                || (std::floor(outgoing_edge_enum_curr/2.0)>=std::floor(incoming_edge_enum_other/2.0))){//curr is middle vertex
                  boost::add_edge(pathidx, other_path, conflict_graph);
                  if(params.verbose) std::cout << "adding edge" <<std::endl;
                }     
              }else{
                //middle vertex
                std::pair<SystemGraph::edge_descriptor, bool> edge_in=boost::edge(path1[res-path1.begin()-1], u, system_graph); 
                std::pair<SystemGraph::edge_descriptor, bool> edge_out=boost::edge(u, path1[res-path1.begin()+1], system_graph); 
                //get enumeration of overlapping path
                SystemGraph::out_edge_iterator ei, ei_end;
                int idx=0;
                int flag=0;
                for (std::tie(ei, ei_end) = boost::out_edges(u, system_graph); ei != ei_end; ++ei){
                  if(system_graph[*ei].edge_id==system_graph[edge_in.first].edge_id){                  
                    incoming_edge_enum_other=system_graph[u].enumeration[idx];
                    flag++;
                  }else if(system_graph[*ei].edge_id==system_graph[edge_out.first].edge_id){                  
                    outgoing_edge_enum_other=system_graph[u].enumeration[idx];
                    flag++;
                  }     
                  if(flag==2) break; //enumeration of both edges found    
                  idx++;       
                } 
                //add edges to conflict graph
                int mincurr=std::min(std::floor(outgoing_edge_enum_curr/2.0), std::floor(incoming_edge_enum_curr/2.0));
                int maxcurr=std::max(std::floor(outgoing_edge_enum_curr/2.0), std::floor(incoming_edge_enum_curr/2.0));
                int minother=std::min(std::floor(outgoing_edge_enum_other/2.0), std::floor(incoming_edge_enum_other/2.0));
                int maxother=std::max(std::floor(outgoing_edge_enum_other/2.0), std::floor(incoming_edge_enum_other/2.0));
                if((incoming_edge_enum_curr==-99 && std::floor(outgoing_edge_enum_curr/2.0)>=std::floor(incoming_edge_enum_other/2.0)) //curr is starting vertex
                || (outgoing_edge_enum_curr==-99 && std::floor(incoming_edge_enum_curr/2.0)<=std::floor(outgoing_edge_enum_other/2.0)) //curr is ending vertex
                || (maxcurr>maxother && mincurr<=maxother) || (maxcurr<=maxother && minother<=maxcurr) ){//curr is middle vertex
                  boost::add_edge(pathidx, other_path, conflict_graph);
                  if(params.verbose) std::cout << "adding edge" <<std::endl;
                }            
              }              
            }
            cnt++;
          }   
        }

        conflict_graph[pathidx].node_id=pathidx;        
        conflict_graph[pathidx].path=mypath;
        pathidx++;
      }

      //color conflict graph
      // use a greedy Welsh-Powell algorithm (color according to max degree, largest first)
      std::vector<vertices_size_type> color_vec(num_vertices(conflict_graph));
      ColorMap color(&color_vec.front(), get(boost::vertex_index, conflict_graph));
      std::vector<int> order_vec(num_vertices(conflict_graph));
      std::iota(order_vec.begin(), order_vec.end(), 0);
      std::vector<int> degrees(num_vertices(conflict_graph));
      int temp_ctr=0;
      for (auto v : boost::make_iterator_range(boost::vertices(conflict_graph))){
        degrees[temp_ctr]=boost::out_degree(v, conflict_graph);
        temp_ctr++;
      }
      std::stable_sort(order_vec.begin(), order_vec.end(), [&degrees](size_t i1, size_t i2) {return degrees[i1] > degrees[i2];});
      /*for(auto idx : order_vec) std::cout << idx << " ";
      std::cout << std::endl;
      for(auto idx : degrees) std::cout << idx << " ";
      std::cout << std::endl;*/
      DegreeMap order(&order_vec.front(), get(boost::vertex_index, conflict_graph));
      vertices_size_type num_colors = boost::sequential_vertex_coloring(conflict_graph, order, color);
      if(params.verbose) std::cout << "Chromatic #: " << num_colors << std::endl;
      chromatic_numbers.push_back(num_colors);
     


      //print conflict graph if desired
      if(params.output_conflict_graph){
        std::ofstream fs;
        if(output_data_filename==NULL){
            fs.open("./out/conflict.dot");      
        }else{
          std::string output=output_data_filename;       
          fs.open(output+"_"+std::to_string(i)+".dot");
        }
        try{
          boost::dynamic_properties conflict_graph_dynamic_properties(boost::ignore_other_properties);
          conflict_graph_dynamic_properties.property("node_id", get(&ConflictGraphVertexProperties::node_id, conflict_graph));
          conflict_graph_dynamic_properties.property("label", get(&ConflictGraphVertexProperties::label, conflict_graph));
          for (auto v : boost::make_iterator_range(boost::vertices(conflict_graph))){ 
              conflict_graph[v].fillcolor=std::to_string(color[v]+1); 
          }       
          conflict_graph_dynamic_properties.property("fillcolor", get(&ConflictGraphVertexProperties::fillcolor, conflict_graph));
          conflict_graph_dynamic_properties.property("style", get(&ConflictGraphVertexProperties::style, conflict_graph));
          conflict_graph_dynamic_properties.property("colorscheme", get(&ConflictGraphVertexProperties::colorscheme, conflict_graph));
          
          boost::write_graphviz_dp(fs, conflict_graph, conflict_graph_dynamic_properties, "node_id");
        }catch(...){
            std::cout <<  "Writing conflict graph failed. Check file location exists."<< std::endl; 
            return 0;
        }
      }
    }

    //get best chromatic number found
    double min = *min_element(chromatic_numbers.begin(), chromatic_numbers.end());
    std::cout<<"Best chromatic number: "<<min<<std::endl;
    const size_t sz = chromatic_numbers.size();
    double mean = std::accumulate(chromatic_numbers.begin(), chromatic_numbers.end(), 0.0) / sz;
    std::cout<<"Mean chromatic number: "<<mean<<std::endl;
    // Now calculate the variance
    auto variance_func = [&mean, &sz](double accumulator, const double& val) {
        return accumulator + ((val - mean)*(val - mean) / (sz - 1));
    };
    double variance = std::accumulate(chromatic_numbers.begin(), chromatic_numbers.end(), 0.0, variance_func);
    std::cout<<"Variance chromatic number: "<<variance<<std::endl;

    //output results to file if requested
    std::ofstream ofile_stream;
    if(!(output_data_filename==NULL)){
        std::string output=output_data_filename;       
        ofile_stream.open(output+".csv", std::ios_base::app); //open for appending
        ofile_stream << min << "," << mean << "," << variance << std::endl;
    }

  }

