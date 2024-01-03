/*
 *  author: Suhas Vittal
 *  date:   2 January 2024
 * */

#ifndef QONTRA_CXX_GRAPHVIZ_h
#define QONTRA_CXX_GRAPHVIZ_h

#include <graphviz/gvc.h>

#include <string>

// This is glue code for modern C++ with cgraph and graphviz.

Agraph_t*   cxx_agopen(std::string, Agdesc_t, Agdisc_t*);
Agnode_t*   cxx_agnode(Agraph_t*, std::string, int);
Agedge_t*   cxx_agedge(Agraph_t*, Agnode_t*, Agnode_t*, std::string, int);

template <class AGPTR>
void    cxx_agset(AGPTR, std::string, std::string);

void    copy_cxx_string_to_buffer(std::string, char*);

#include "cxxviz.inl"

#endif  // QONTRA_CXX_GRAPHVIZ_h
