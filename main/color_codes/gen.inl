/*
 *  author: Suhas Vittal
 *  date:   23 May 2024
 * */

DetailedStimCircuit
make_capacity(
        TannerGraph* tgr, 
        fp_t p,
        bool is_memory_x,
        const std::map<sptr<tanner::vertex_t>, int>& color_map)
{
    DetailedStimCircuit circuit;
    // Measure observables and stabilizers once to reach valid state.
    size_t nobs = 0;
    for (auto obs : tgr->get_obs(is_memory_x)) {
        std::vector<uint32_t> operands;
        uint32_t flags = is_memory_x ? stim::TARGET_PAULI_X_BIT : stim::TARGET_PAULI_Z_BIT;
        bool first = true;
        for (sptr<tanner::vertex_t> tv : obs) {
            if (!first) operands.push_back(stim::TARGET_COMBINER);
            operands.push_back( static_cast<uint32_t>( tv->id ) | flags);
            first = false;
        }
        circuit.safe_append_u("MPP", operands);
        circuit.safe_append_ua("OBSERVABLE_INCLUDE", { stim::TARGET_RECORD_BIT | 1 },
                static_cast<fp_t>(nobs));
        nobs++;
    }
    std::vector<uint32_t> ch_operands;
    for (sptr<tanner::vertex_t> tch : tgr->get_checks()) {
        std::vector<uint32_t> sub_operands;
        bool is_x = tch->qubit_type == tanner::vertex_t::type::xparity;
        uint32_t flags = is_x ? stim::TARGET_PAULI_X_BIT : stim::TARGET_PAULI_Z_BIT;
        for (sptr<tanner::vertex_t> tv : tgr->get_neighbors(tch)) {
            sub_operands.push_back( static_cast<uint32_t>( tv->id ) | flags);
        }
        for (size_t i = 0; i < sub_operands.size(); i++) { 
            if (i > 0) ch_operands.push_back(stim::TARGET_COMBINER);
            ch_operands.push_back(sub_operands.at(i));
        }
    }
    circuit.safe_append_u("MPP", ch_operands);
    // Inject error on data qubits.
    std::vector<uint32_t> data_operands;
    for (sptr<tanner::vertex_t> tv : tgr->get_vertices_by_type(tanner::vertex_t::type::data)) {
        data_operands.push_back( static_cast<uint32_t>( tv->id ));
    }
    if (is_memory_x) {
        circuit.safe_append_ua("Z_ERROR", data_operands, p);
    } else {
        circuit.safe_append_ua("X_ERROR", data_operands, p);
    }
    // Measure stabilizers.
    circuit.safe_append_u("MPP", ch_operands);
    // Make detectors.
    uint32_t ndet = 0;
    uint32_t nch = 0;
    uint32_t totch = tgr->get_checks().size();
    for (sptr<tanner::vertex_t> tch : tgr->get_checks()) {
        bool is_x = tch->qubit_type == tanner::vertex_t::type::xparity;
        if (is_x == is_memory_x) {
            std::vector<uint32_t> operands{
                (stim::TARGET_RECORD_BIT | (totch-nch)),
                (stim::TARGET_RECORD_BIT | (2*totch-nch))
            };
            circuit.safe_append_u("DETECTOR", operands,
                    {0, 0, 0, static_cast<fp_t>( color_map.at(tch) )});
            circuit.detector_color_map[ndet] = color_map.at(tch);
            circuit.detector_base_map[ndet] = ndet;
            ndet++;
        }
        nch++;
    }
    nobs = 0;
    for (auto obs : tgr->get_obs(is_memory_x)) {
        std::vector<uint32_t> operands;
        uint32_t flags = is_memory_x ? stim::TARGET_PAULI_X_BIT : stim::TARGET_PAULI_Z_BIT;
        bool first = true;
        for (sptr<tanner::vertex_t> tv : obs) {
            if (!first) operands.push_back(stim::TARGET_COMBINER);
            operands.push_back( static_cast<uint32_t>( tv->id ) | flags);
            first = false;
        }
        circuit.safe_append_u("MPP", operands);
        circuit.safe_append_ua("OBSERVABLE_INCLUDE", { stim::TARGET_RECORD_BIT | 1 },
                static_cast<fp_t>(nobs));
        nobs++;
    }
    circuit.number_of_colors_in_circuit = 3;
    return circuit;
}

