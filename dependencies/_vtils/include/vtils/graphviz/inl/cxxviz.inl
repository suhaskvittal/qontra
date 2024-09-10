/*
 *  author: Suhas Vittal
 *  date:   2 January 2024
 * */

#include <string.h>

inline Agraph_t*
cxx_agopen(std::string graph_name, Agdesc_t edge_type, Agdisc_t* _unknown_param) {
    char c_name[1024];
    copy_cxx_string_to_buffer(graph_name, c_name);
    return agopen(c_name, edge_type, _unknown_param);
}

inline Agnode_t*
cxx_agnode(Agraph_t* graph, std::string node_name, int create_flag) {
    char c_name[1024];
    copy_cxx_string_to_buffer(node_name, c_name);
    return agnode(graph, c_name, create_flag);
}

inline Agedge_t*
cxx_agedge(Agraph_t* graph, Agnode_t* src, Agnode_t* dst, std::string edge_name, int create_flag) {
    char c_name[1024];
    copy_cxx_string_to_buffer(edge_name, c_name);
    return agedge(graph, src, dst, c_name, create_flag);
}

template <class AGPTR> inline void
cxx_agset(AGPTR ax, std::string attr_name, std::string attr_value) {
    char c_name[1024];
    char c_value[1024];
    copy_cxx_string_to_buffer(attr_name, c_name);
    copy_cxx_string_to_buffer(attr_value, c_value);
    agsafeset(ax, c_name, c_value, "");
}

inline void
copy_cxx_string_to_buffer(std::string s, char* buf) {
    strcpy(buf, s.c_str()); 
}

