/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#include <vtils/filesystem.h>

namespace qontra {
namespace protean {

inline void
write_network_to_folder(std::string output_folder, PhysicalNetwork* network) {
    // Make folder if it does not exist.
    vtils::safe_create_directory(output_folder);

    const std::string stats_file = output_folder + "/stats.txt";
    const std::string x_schedule_file = output_folder + "/memory_x.qes";
    const std::string z_schedule_file = output_folder + "/memory_z.qes";
    const std::string coupling_file = output_folder + "/coupling_graph.txt";
    const std::string role_file = output_folder + "/roles.txt";
    const std::string tanner_graph_file = output_folder + "/tanner_graph.txt";
    const std::string flag_assign_file = output_folder + "/flag_assignment.txt";

    network->write_stats_file(stats_file);
    network->write_coupling_file(coupling_file);
    network->write_role_file(role_file);
    network->write_tanner_graph_file(tanner_graph_file);
    network->write_flag_assignment_file(flag_assign_file);
    network->write_schedule_file(x_schedule_file, true);
    network->write_schedule_file(z_schedule_file, false);
}

}   // protean
}   // qontra

