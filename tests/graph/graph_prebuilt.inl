/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

template <size_t N> inline GRAPH
make_complete_graph() {
    GRAPH gr = make_graph_with_k_vertices(N);
    std::vector<sptr<VERTEX>> vertices = gr.get_vertices();
    for (size_t i = 0; i < N; i++) {
        sptr<VERTEX> v = vertices[i];
        for (size_t j = i+1; j < N; j++) {
            sptr<VERTEX> w = vertices[j];
            gr.add_edge(fast_make_edge(v, w));
        }
    }
    return gr;
}

inline GRAPH
fast_make_graph_with_k_vertices(size_t k) {
    GRAPH gr;
    for (size_t i = 0; i < k; i++) {
        sptr<VERTEX> v = std::make_shared<>();
        v->id = i;
        gr.add_vertex(v);
    }
    return gr;
}

inline sptr<EDGE>
fast_make_edge(sptr<VERTEX> v, sptr<VERTEX> w) {
    sptr<EDGE> e = std::make_shared<>();
    e->src = v;
    e->dst = w;
    e->is_undirected = true;
    return e;
}
