/*
 *  author: Suhas Vittal
 *  date:   19 May 2023
 * */

#ifndef SIM_h
#define SIM_h

#include "defs.h"
#include "tables.h"
#include "instruction.h"

    
const uint64_t DMEM_SIZE = 1 << 14;     // 16 KB

namespace contra {

enum class Code { surf, hhex };

// Flags for enabling certain policies:
// 
// Implemented so far:
//  (0x1)   ERASER
//
// Use accordingly for your own research.

const uint8_t CONTRASIM_EN_ERASER = 0x1;

#define PR_FLAG(f)      if (this->sim_flags & (f)) {
#define END_PR          }

class ContraSim {
public:
    ContraSim(uint distance, fp_t p);
    ContraSim(uint distance, fp_t p, Code);


    stim::Circuit   get_canonical_circuit(void);    // Gets circuit to be consumed
                                                    // by decoder given *standard*
                                                    // calibration data. Thus,
                                                    // the canonical circuit only
                                                    // includes information about
                                                    // 1q errors, 2q errors, readout
                                                    // errors, and coherence/idling
                                                    // errors.
    void            load_decoder(/* TODO */);       // Loads decoder to use during
                                                    // simulation.
    void            load_flags(uint64_t);           // Loads flags to be used during
                                                    // simulation.
    void            load_imem(std::string);         // Loads instruction memory
                                                    // (assumed to be all SRAM) with
                                                    // control instructions.
    void            load_imem_with_memory_x(void);
    void            load_imem_with_memory_y(void);  // Loads instruction memory
                                                    // with control instructions for
                                                    // memory_Y experiment. Final
                                                    // measurement is replaced with
                                                    // PEEKZ and PEEKX operations.
    void            load_imem_with_memory_z(void);

    void            run(uint rounds);

    std::string     dump_stats(void);
    void            write_stats_to_file(std::string);

    bool en_real_time_constraint = false;   // If enabled, then operations may be
                                            // stalled if decoder cannot commit to
                                            // Pauli frame in time.
    bool en_transmission_latency = false;   // If enabled, then communication between
                                            // controls and FTQC has latency according
                                            // to transmission bandwidth.
    // Parameters
    fp_t transmission_bandwidth = 100e6;  // Bytes per second (Bps)
    fp_t control_clk_freq = 250e6;  // in Hz
    
    bool en_mpi = false;
    // Statistics
    uint64_t n_decodes;
    uint64_t n_logical_errors_x;
    uint64_t n_logical_errors_z;
    uint64_t normalized_latency;    // Normalized to depth * (round latency).
protected:
    // Communication pattern every round:
    //  FTQC                CP
    //      | S --------> R |       # Receives syndrome
    //      |        Decode |       # Decodes the syndrome history
    //      |        Commit |       # Updates Pauli Frame
    //         ROUND START
    //      | R <-------- S |
    //      | R <-------- S |
    //      |      ...      |
    //      | R <-------- S |       # Measurement of stabilizers
    //          ROUND END
    //      | S --------> R |       # Receives syndrome
    //      |      ...      |
    // FTQC Blocks
    void            qc_blk_recv_cmd(void);      // Receives operation from controls and
                                                // executes it. This may modify the
                                                // quantum state.
    void            qc_blk_send_syndrome(void); // Sends the contents of the syndrome
                                                // buffer to the control processor.
    // Control Processor + Decoder Blocks
    void            cp_blk_recv_syndrome(void); // Receives syndrome from the FTQC
                                                // and updates the syndrome history.
    void            cp_blk_decode(void);        // Decodes the syndrome history and
                                                // updates the Pauli frames.
    void            cp_blk_send_cmd(void);      // Sends commands to the control
                                                // processor.

    const uint distance;
    const fp_t p;

    Code    qec_code;
private:
    struct Register {
        uint64_t data = 0;
        bool valid = false;
    };

    // All tables below are 2d tables such that row-major tables
    // have the first index as the shot/index number and column
    // major tables have the second index as the shot/index number.
    //
    // FTQC structures
    stim::simd_bit_table    qc_syndrome_buf;        // Row-major
    // Control Processor + Decoder structures
    stim::simd_bit_table    cp_syndrome_history;    // Row-major
    stim::simd_bit_table    cp_pauli_frames;        // Column-major

    stim::simd_bit_table    cp_pauli_frames_ideal;  // Column-major, if we had a
                                                    // perfect decoder that never
                                                    // failed. A logical error
                                                    // occurs when the Pauli frames
                                                    // do not match this ideal.

    uint64_t    pc;

    std::array<Register, 16>        register_file;
    std::vector<cp::Instruction>    imem;
    std::array<uint64_t, DMEM_SIZE> dmem;
    
    // Transmission structures
    //
    // Assume FTQCs can execute operations as if they were SIMD. For example,
    // we can specify an H gate on qubits 1, 2, 3, 5, 8, 10 in one instruction.
    std::vector<qc::Instruction> instruction_buf;

    uint64_t sim_flags;
};

} // contra

#endif  // SIM_h
