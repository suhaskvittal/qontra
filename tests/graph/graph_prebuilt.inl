/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

inline GRAPH
make_complete_graph(size_t n) {
    GRAPH gr = fast_make_graph_with_k_vertices(n);
    std::vector<sptr<VERTEX>> vertices = gr.get_vertices();
    for (size_t i = 0; i < n; i++) {
        sptr<VERTEX> v = vertices[i];
        for (size_t j = i+1; j < n; j++) {
            sptr<VERTEX> w = vertices[j];
            gr.add_edge(fast_make_edge(v, w));
        }
    }
    return gr;
}

inline GRAPH
make_bipartite_complete_graph(size_t n1, size_t n2) {
    GRAPH gr = fast_make_graph_with_k_vertices(n1+n2);
    for (size_t i = 0; i < n1; i++) {
        sptr<VERTEX> v = gr.get_vertex(i);
        for (size_t j = n1; j < n1+n2; j++) {
            sptr<VERTEX> w = gr.get_vertex(j);
            gr.add_edge(fast_make_edge(v, w));
        }
    }
    return gr;
}

inline GRAPH
fast_make_graph_with_k_vertices(size_t k) {
    GRAPH gr;
    for (size_t i = 0; i < k; i++) {
        sptr<VERTEX> v = std::make_shared<VERTEX>();
        v->id = i;
        gr.add_vertex(v);
    }
    return gr;
}

inline sptr<EDGE>
fast_make_edge(sptr<VERTEX> v, sptr<VERTEX> w) {
    sptr<EDGE> e = std::make_shared<EDGE>();
    e->src = v;
    e->dst = w;
    e->is_undirected = true;
    return e;
}
